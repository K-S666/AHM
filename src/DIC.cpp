#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double dina_logL(const arma::mat y,
                 unsigned int N,
                 unsigned int J,
                 unsigned int K,
                 const arma::vec s,
                 const arma::vec g,
                 const arma::mat alpha,
                 const arma::mat Q)
{
  double logL = 0.0;
  arma::mat AQ = alpha * Q.t();

  for (unsigned int j = 0; j < J; j++) {
    arma::vec ETA(N, arma::fill::zeros);
    double qq = arma::as_scalar(Q.row(j) * Q.row(j).t());
    ETA.elem(arma::find(AQ.col(j) == qq)).fill(1.0);

    for (unsigned int i = 0; i < N; i++) {
      double eta_ij = ETA(i);
      double p1 = (1.0 - s(j)) * eta_ij + g(j) * (1.0 - eta_ij);
      double p0 = s(j) * eta_ij + (1.0 - g(j)) * (1.0 - eta_ij);

      // Guard against exact 0/1 draws creating log(0) during diagnostics.
      p1 = std::min(std::max(p1, 1e-12), 1.0 - 1e-12);
      p0 = std::min(std::max(p0, 1e-12), 1.0 - 1e-12);
      logL += std::log(p1) * y(i, j) + std::log(p0) * (1.0 - y(i, j));
    }
  }

  return logL;
}

// [[Rcpp::export]]
arma::vec dina_logL_vec(const arma::mat y,
                        unsigned int N,
                        unsigned int J,
                        unsigned int K,
                        int chain,
                        const arma::mat s,
                        const arma::mat g,
                        const arma::cube alpha,
                        const arma::cube Q)
{
  arma::vec logL(chain);
  for (int i = 0; i < chain; i++) {
    logL(i) = dina_logL(y, N, J, K, s.col(i), g.col(i),
                        alpha.slice(i), Q.slice(i));
  }
  return logL;
}

// [[Rcpp::export]]
double bar_Dev_theta(const arma::mat y,
                     unsigned int N,
                     unsigned int J,
                     unsigned int K,
                     int chain,
                     const arma::mat s,
                     const arma::mat g,
                     const arma::cube alpha,
                     const arma::cube Q)
{
  arma::vec logL = dina_logL_vec(y, N, J, K, chain, s, g, alpha, Q);
  return -2.0 * arma::mean(logL);
}

// [[Rcpp::export]]
double compute_DIC(const arma::mat y,
                   unsigned int N,
                   unsigned int J,
                   unsigned int K,
                   int chain,
                   const arma::mat s,
                   const arma::mat g,
                   const arma::cube alpha,
                   const arma::cube Q)
{
  arma::vec logL = dina_logL_vec(y, N, J, K, chain, s, g, alpha, Q);
  double mean_deviance = -2.0 * arma::mean(logL);
  double deviance_at_best_draw = arma::min(-2.0 * logL);
  double effective_parameters = mean_deviance - deviance_at_best_draw;
  return deviance_at_best_draw + 2.0 * effective_parameters;
}

// [[Rcpp::export]]
arma::vec dina_logL_vec_fixedQ(const arma::mat y,
                               unsigned int N,
                               unsigned int J,
                               unsigned int K,
                               int chain,
                               const arma::mat s,
                               const arma::mat g,
                               const arma::cube alpha,
                               const arma::mat Q)
{
  arma::vec logL(chain);
  for (int i = 0; i < chain; i++) {
    logL(i) = dina_logL(y, N, J, K, s.col(i), g.col(i),
                        alpha.slice(i), Q);
  }
  return logL;
}

// [[Rcpp::export]]
double compute_DIC_fixedQ(const arma::mat y,
                          unsigned int N,
                          unsigned int J,
                          unsigned int K,
                          int chain,
                          const arma::mat s,
                          const arma::mat g,
                          const arma::cube alpha,
                          const arma::mat Q)
{
  arma::vec logL = dina_logL_vec_fixedQ(y, N, J, K, chain, s, g, alpha, Q);
  double mean_deviance = -2.0 * arma::mean(logL);
  double deviance_at_best_draw = arma::min(-2.0 * logL);
  double effective_parameters = mean_deviance - deviance_at_best_draw;
  return deviance_at_best_draw + 2.0 * effective_parameters;
}
