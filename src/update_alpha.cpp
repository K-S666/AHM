#include "ahm_common.h"

using namespace Rcpp;

arma::mat AHM_alpha_update(arma::mat y,arma::vec s,arma:: vec g,arma::mat Q,arma::mat G,
                    arma::mat alpha_all,arma::vec pi0,int N,int J,int K,int L)
{
  arma::mat R;
  arma::vec alpha_current_possible_index;
  arma::vec alpha_current_possible_code;
  //arma::mat alpha_current_possible;
  arma::mat alpha_all_current;
  Rcpp::List rr;
  R=Reachability(G,K);
  rr = Reduced_alpha(alpha_all,R,K);
  alpha_current_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_current_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  //alpha_current_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  alpha_all_current = Trans_10to2_mat(K,alpha_current_possible_index);
  
  int Cg=alpha_current_possible_index.n_elem;
  arma::mat log_post(N,Cg,arma::fill::zeros);
  arma::vec pi0_new(Cg,arma::fill::zeros);
  arma::mat alpha_new(N,K);
  arma::ivec alpha_class(N);
  arma::vec aa;
  for(int c=0;c<Cg;c++){
    aa=pi0.elem(arma::find(alpha_current_possible_code==alpha_current_possible_index(c)));
    pi0_new(c)=sum(aa);
    for(int i=0;i<N;i++){
      log_post(i,c)=std::log(std::max(pi0_new(c), 1e-12));
      for(int j=0;j<J;j++){
        double eta_ij = dina_eta_row(alpha_all_current.row(c), Q.row(j));
        log_post(i,c) += dina_log_response_prob(y(i,j), eta_ij, s(j), g(j));
      }
    }
  }
  for(int i=0;i<N;i++){
    // Normalize on the log scale to avoid underflow when J is large.
    alpha_class(i) = rgen::rmultinomial(ahm_softmax_from_log((log_post.row(i)).t()));
    alpha_new.row(i) = alpha_all_current.row(alpha_class(i));
  }
  return alpha_new;
}

arma::vec rDirichlet(const arma::vec& deltas){
  unsigned int C = deltas.n_elem;
  arma::vec Xgamma(C);
  
  //generating gamma(deltac,1)
  for(unsigned int c=0;c<C;c++){
    Xgamma(c) = R::rgamma(deltas(c),1.0);
  }
  return Xgamma/sum(Xgamma);
}

Rcpp::List AHM_alpha_pi_update(arma::mat y,arma:: vec s,arma:: vec g,arma::mat Q,arma::mat G,
                        arma::mat alpha_all,arma::vec pi0,arma::vec delta0,int N,int J,int K,int L)
{
  arma::mat R;
  arma::vec alpha_current_possible_index;
  arma::vec alpha_current_possible_code;
  //arma::mat alpha_current_possible;
  arma::mat alpha_all_current;
  Rcpp::List rr;
  R=Reachability(G,K);
  rr = Reduced_alpha(alpha_all,R,K);
  alpha_current_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_current_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  //alpha_current_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  alpha_all_current = Trans_10to2_mat(K,alpha_current_possible_index);
  
  int Cg=alpha_current_possible_index.n_elem;
  arma::mat log_post(N,Cg,arma::fill::zeros);
  arma::vec pi0_new(Cg,arma::fill::zeros);
  arma::vec delta0_new(L,arma::fill::zeros);
  arma::mat alpha_new(N,K);
  arma::ivec alpha_class(N);
  arma::vec aa;
  arma::vec count(L,arma::fill::zeros);
  for(int c=0;c<Cg;c++){
    aa=pi0.elem(arma::find(alpha_current_possible_code==alpha_current_possible_index(c)));
    pi0_new(c)=sum(aa);
    aa=delta0.elem(arma::find(alpha_current_possible_code==alpha_current_possible_index(c)));
    delta0_new(alpha_current_possible_index(c))=sum(aa);
    for(int i=0;i<N;i++){
      log_post(i,c)=std::log(std::max(pi0_new(c), 1e-12));
      for(int j=0;j<J;j++){
        double eta_ij = dina_eta_row(alpha_all_current.row(c), Q.row(j));
        log_post(i,c) += dina_log_response_prob(y(i,j), eta_ij, s(j), g(j));
      }
    }
  }
  int kk;
  for(int i=0;i<N;i++){
    // Sample alpha_i from its categorical posterior after log-scale
    // normalization. This preserves the same Gibbs update but avoids
    // multiplying many small item probabilities.
    alpha_class(i) = rgen::rmultinomial(ahm_softmax_from_log((log_post.row(i)).t()));
    kk=alpha_class(i);
    alpha_new.row(i) = alpha_all_current.row(alpha_class(i));
    count(alpha_current_possible_index(kk))++;
  }
  arma::vec pi_new=rDirichlet(count+delta0_new);
  return Rcpp::List::create(Rcpp::Named("alpha_new") = alpha_new,
                            Rcpp::Named("pi_new") = pi_new,
                            Rcpp::Named("count") = count,
                            Rcpp::Named("delta0_new") = delta0_new,
                            Rcpp::Named("alpha_class") = alpha_class,
                            Rcpp::Named("alpha_current_possible_index")=alpha_current_possible_index);
}



