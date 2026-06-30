#include "ahm_common.h"
#include <cstdlib>
#include <fstream>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List AHM_fixedQ(arma::mat Y,
                          arma::mat Q,
                          int N1 = 128,
                          int chain_length = 20000,
                          int burn_in = 10000,
                          double a_s0 = 1.0,
                          double a_g0 = 1.0,
                          double b_s0 = 1.0,
                          double b_g0 = 1.0,
                          double p_add = 0.5,
                          int chain_id = 1,
                          bool progress = true,
                          int print_every = 1000)
{
  const int N = Y.n_rows;
  const int J = Y.n_cols;
  const int K = Q.n_cols;
  const int L = std::pow(2, K);
  const int D = L - 1;
  const int T = chain_length - burn_in;
  const int progress_every = std::max(1, print_every);
  const char* progress_dir_env = std::getenv("AHM_PROGRESS_DIR");
  std::string progress_file;
  if (progress_dir_env != nullptr && std::string(progress_dir_env).size() > 0) {
    progress_file = std::string(progress_dir_env) + "/chain_" + std::to_string(chain_id) + ".txt";
  }

  arma::cube GG(K, K, T);
  arma::cube AA(N, K, T);
  arma::mat SLIP(J, T);
  arma::mat GUESS(J, T);
  arma::mat PIS(L, T);

  arma::vec classvec = arma::linspace<arma::vec>(0, L - 1, L);
  arma::mat alpha_all = Trans_10to2_mat(K, classvec);
  arma::mat Q_all = alpha_all.rows(1, L - 1);
  arma::vec delta0 = arma::ones(L, 1);

  arma::vec s = runif(J) * 0.5;
  arma::vec g = runif(J) * 0.5;
  arma::vec pi = runif(L);
  pi = pi / sum(pi);
  arma::mat G = arma::zeros(K, K);
  arma::mat alpha(N, K, arma::fill::zeros);

  arma::mat eta(L, D);
  for (int l = 0; l < L; l++) {
    for (int d = 0; d < D; d++) {
      eta(l, d) = dina_eta_row(alpha_all.row(l), Q_all.row(d));
    }
  }

  for (int t = 0; t < chain_length; t++) {
    AHM_update_fixQ(
      Y, alpha, s, g, Q, G, pi, p_add, alpha_all, Q_all, N, J, K, L,
      N1, a_s0, a_g0, b_s0, b_g0, delta0, eta
    );

    if (t >= burn_in) {
      const int tt = t - burn_in;
      SLIP.col(tt) = s;
      GUESS.col(tt) = g;
      AA.slice(tt) = alpha;
      GG.slice(tt) = G;
      PIS.col(tt) = pi;
    }

    if (t % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (progress && ((t + 1) == 1 || ((t + 1) % progress_every == 0) || ((t + 1) == chain_length))) {
      if (progress_file.empty()) {
        Rcpp::Rcout << "AHM chain " << chain_id << " iteration "
                    << (t + 1) << "/" << chain_length << "\n";
      } else {
        std::ofstream out(progress_file.c_str(), std::ios::trunc);
        if (out.is_open()) {
          out << (t + 1) << std::endl;
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("alpha") = AA,
    Rcpp::Named("s") = SLIP,
    Rcpp::Named("g") = GUESS,
    Rcpp::Named("pi") = PIS,
    Rcpp::Named("G") = GG
  );
}

