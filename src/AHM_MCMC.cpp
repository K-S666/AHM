#include "ahm_common.h"

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
      ahm_write_progress_file(chain_id, t + 1);
      Rcpp::Rcout << "\rAHM chain " << chain_id << " iteration "
                  << (t + 1) << "/" << chain_length;
      if ((t + 1) == chain_length) Rcpp::Rcout << "\n";
      Rcpp::Rcout.flush();
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
