#include "ahm_common.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat AHM_Q_update(arma::mat y,arma::vec s,arma::vec g,arma::mat alpha,arma::mat Q_all,
                arma::mat Q,int N,int J,int K,arma::vec alpha_current_possible_index,
                arma::mat eta,arma::vec alpha_code)
{
  int D=pow(2,K)-1;
  int CC=alpha_current_possible_index.n_elem;
  arma::vec yj;
  arma::vec yy;
  arma::vec Q_class(J);
  arma::mat Q_new(J,K);
  arma::vec xi_j(K);
  arma::vec log_prior(D);
  arma::vec log_post(D);
  arma::vec p_xi(D);
  for(int j=0;j<J;j++){
    yj=y.col(j);
    for(int k=0;k<K;k++){
      xi_j(k)=R::rbeta(1+Q(j,k),2-Q(j,k));
    }
    for(int c=0;c<D;c++){
      log_prior(c)=0.0;
      for(int k=0;k<K;k++){
        // Chung-style auxiliary prior for each candidate q-vector.  This is
        // a product of Bernoulli probabilities, so it is accumulated on the
        // log scale before posterior normalization.
        double q_ck = Q_all(c,k);
        double p = q_ck > 0.5 ? xi_j(k) : 1.0 - xi_j(k);
        log_prior(c) += std::log(std::max(p, 1e-12));
      }
      log_post(c)=log_prior(c);
      for(int cg=0;cg<CC;cg++){
        yy=yj.elem(arma::find(alpha_code==alpha_current_possible_index(cg)));
        double eta_cj = eta(alpha_current_possible_index(cg),c);
        double n_correct = sum(yy);
        double n_total = yy.n_elem;
        log_post(c) +=
          n_correct * dina_log_response_prob(1.0, eta_cj, s(j), g(j)) +
          (n_total - n_correct) * dina_log_response_prob(0.0, eta_cj, s(j), g(j));
      }
    }
    // Softmax(log likelihood + log prior) gives the same categorical
    // distribution as the original ratio formula, without underflow.
    p_xi = ahm_softmax_from_log(log_post);
    Q_class(j) = rgen::rmultinomial(p_xi);
    Q_new.row(j) = Q_all.row(Q_class(j));
  }
  return Q_new;
}


