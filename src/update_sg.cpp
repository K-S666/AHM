#include "ahm_common.h"

using namespace Rcpp;

Rcpp::List AHM_sg_update(const arma::mat &Y, const arma::mat &Q,
                  const arma::mat &ALPHAS, const arma::vec &ss_old,
                  double as0, double bs0, double ag0, double bg0)
{
  
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  
  arma::vec ETA;
  arma::vec ss_new(J);
  arma::vec gs_new(J);
  arma::mat AQ = ALPHAS * Q.t();
  double T, S, G, y_dot_eta, qq, ps, pg;
  double ug, us;
  
  for (unsigned int j = 0; j < J; j++) {
    us = R::runif(0, 1);
    ug = R::runif(0, 1);
    ETA = arma::zeros<arma::vec>(N);
    qq = arma::as_scalar(Q.row(j) * (Q.row(j)).t());
    ETA.elem(arma::find(AQ.col(j) == qq)).fill(1.0);
    
    y_dot_eta = arma::as_scalar((Y.col(j)).t() * ETA);
    T = sum(ETA);
    S = T - y_dot_eta;
    G = sum(Y.col(j)) - y_dot_eta;
    
    // sample s and g as linearly truncated bivariate beta
    
    // draw g conditoned upon s_t-1
    pg = R::pbeta(1.0 - ss_old(j), G + ag0, N - T - G + bg0, 1, 0);
    gs_new(j) = R::qbeta(ug * pg, G + ag0, N - T - G + bg0, 1, 0);
    // draw s conditoned upon g
    ps = R::pbeta(1.0 - gs_new(j), S + as0, T - S + bs0, 1, 0);
    ss_new(j) = R::qbeta(us * ps, S + as0, T - S + bs0, 1, 0);
  }
  return Rcpp::List::create(Rcpp::Named("ss_new") = ss_new,
                            Rcpp::Named("gs_new") = gs_new);
}

