#include "ahm_common.h"

using namespace Rcpp;

Rcpp::List add_path_sample1(arma::mat R,arma::mat G,int K)// sample an added edge
{
  arma::mat G_new=G;
  arma::mat R_double = R+R.t();
  //R_double.elem(arma::find(R_double > 0)).fill(1.0);
  arma::uvec index=arma::find(R_double == 0);
  int n_edge=index.n_elem;
  int change=1;
  int sample_location;
  double prob=0;
  if(n_edge==0){
    change = 0;
  }else{
    prob=1.0/n_edge;
    sample_location=rgen::rmultinomial(arma::ones(n_edge,1)/n_edge);
    G_new(index(sample_location))=1;
  }
  return Rcpp::List::create(Rcpp::Named("G_new") = G_new,Rcpp::Named("change") = change,Rcpp::Named("prob") = prob);
}


Rcpp::List reduce_path_sample1(arma::mat G,int K)// sample an removed edge
{
  arma::uvec index=arma::find(G == 1);
  int n_edge=index.n_elem;
  arma::mat G_new = G;
  int sample_location;
  double prob=0.0;
  int change = 1;
  if(n_edge==0){
    change = 0;
  }else{
    sample_location=rgen::rmultinomial(arma::ones(n_edge,1)/n_edge);
    G_new(index(sample_location))=0;
    prob=1.0/n_edge;
  }
  return Rcpp::List::create(Rcpp::Named("G_new") = G_new,Rcpp::Named("change") = change,Rcpp::Named("prob") = prob);
}


double AHM_G_acceptance_prob(arma::mat y,arma::vec s,arma::vec g,arma::vec pi0,arma::vec delta0,arma::mat alpha,arma::mat alpha_new,
                     arma::mat Q,arma::mat Q_new,arma::mat G,arma::mat G_new,arma::vec alpha_possible_index,
                     arma::vec alpha_new_possible_index,arma::vec alpha_possible_code,arma::vec alpha_new_possible_code,
                     arma::mat alpha_possible,arma::mat alpha_new_possible,double pr_G,int action,double p_add)
{ 
  int N=y.n_rows;
  int J=y.n_cols;
  int K=alpha.n_cols;
  int Cg = alpha_possible_index.n_elem;
  int Cg_new = alpha_new_possible_index.n_elem;
  arma::mat R_new = Reachability(G_new,K);
  arma::vec pi = arma::zeros(Cg,1);
  arma::vec pi_new = arma::zeros(Cg_new,1);
  arma::mat eta(J,Cg);
  arma::mat eta_new(J,Cg_new);
  arma::mat alpha_possible_all;
  arma::mat alpha_new_possible_all;
  
  alpha_possible_all = Trans_10to2_mat(K,alpha_possible_index);
  alpha_new_possible_all = Trans_10to2_mat(K,alpha_new_possible_index);
  for(int c=0;c<Cg;c++){
    pi(c) = sum(pi0.elem(arma::find(alpha_possible_code==alpha_possible_index(c))));
    for(int j=0;j<J;j++){
      eta(j,c)=dina_eta_row(alpha_possible_all.row(c), Q.row(j));
    }
  }
  for(int c=0;c<Cg_new;c++){
    pi_new(c) = sum(pi0.elem(arma::find(alpha_new_possible_code==alpha_new_possible_index(c))));
    for(int j=0;j<J;j++){
      eta_new(j,c)=dina_eta_row(alpha_new_possible_all.row(c), Q_new.row(j));
    }
  }
  
  double log_ratio_y = 0.0;
  for(int i=0;i<N;i++){
    arma::vec log_old(Cg);
    for(int c=0;c<Cg;c++){
      log_old(c) = std::log(std::max(pi(c), 1e-12)) +
        dina_log_response_vector(eta.col(c), y.row(i), s, g);
    }

    arma::vec log_new(Cg_new);
    for(int c=0;c<Cg_new;c++){
      log_new(c) = std::log(std::max(pi_new(c), 1e-12)) +
        dina_log_response_vector(eta_new.col(c), y.row(i), s, g);
    }
    // Marginal likelihood under G is a mixture over G-permissible attribute
    // classes. log-sum-exp keeps this stable for long tests and large N.
    log_ratio_y += ahm_log_sum_exp(log_new) - ahm_log_sum_exp(log_old);
  }

  double pr_G_new;
  double log_ratio_proposal;
  const double log_p_add = std::log(p_add);
  const double log_p_remove = std::log(1.0 - p_add);
  if(action==1){
    arma::uvec index=arma::find(G_new==1);
    int n_edge=index.n_elem;
    pr_G_new = 1.0/n_edge;
    log_ratio_proposal = std::log(pr_G_new) - std::log(pr_G) + log_p_add - log_p_remove;
  }else{
    R_new = R_new+R_new.t();
    arma::uvec index=arma::find(R_new==0);
    int n_edge=index.n_elem;
    pr_G_new = 1.0/n_edge;
    log_ratio_proposal = std::log(pr_G_new) - std::log(pr_G) + log_p_remove - log_p_add;
  }

  double log_accept = log_ratio_y + log_ratio_proposal;
  if (log_accept >= 0.0) {
    return 1.0;
  }
  return std::exp(log_accept);
}

arma::mat update_G(arma::mat y,arma::vec s,arma::vec g,arma::mat alpha,arma::mat alpha_all,arma::mat Q,arma::mat G,
                      arma::mat R,int N,int J,int K,int L,double p_add,arma::vec pi,arma::vec delta0,arma::mat G_new,
                      double prob_edge,int action)
{
  
  
  Rcpp::List rr;
  arma::mat R_new;
  arma::vec alpha_code;
  arma::vec alpha_new_possible_index;
  arma::vec alpha_new_possible_code;
  arma::mat alpha_new_possible;
  arma::vec alpha_possible_index;
  arma::vec alpha_possible_code;
  arma::mat alpha_possible;
  arma::mat alpha_new(N,K);
  
  
  double accept=1.0;
  R_new = Reachability(G_new,K);
  alpha_code = Trans_2to10_mat(alpha, K);
  rr = Reduced_alpha(alpha_all,R_new,K);
  alpha_new_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_new_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  alpha_new_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  for(int i=0;i<N;i++){
    alpha_new.row(i) = alpha_new_possible.row(alpha_code(i));
  }
  
  rr = Reduced_alpha(alpha_all,R,K);
  alpha_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  alpha_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  
  double rand=R::runif(0,1);
  accept = AHM_G_acceptance_prob(y,s,g,pi,delta0,alpha,alpha_new,Q,Q,G,G_new,alpha_possible_index,alpha_new_possible_index,
                         alpha_possible_code,alpha_new_possible_code,alpha_possible,alpha_new_possible,prob_edge,action,p_add);
  if(rand > accept){
    G_new=G;
    action = 0;
  }
  
  
  
  return G_new;
}

arma::mat AHM_G_update(arma::mat y,arma::vec s,arma::vec g,arma::mat alpha,arma::mat alpha_all,arma::mat Q,arma::mat G,
                arma::mat R,int N,int J,int K,int L,double p_add,arma::vec pi0,arma::vec delta0)
{
  Rcpp::List operate;
  arma::mat G_new(K,K);
  int change;
  double prob_edge;
  double u=R::runif(0,1);
  G_new=G;
  if(u < p_add){
    operate = add_path_sample1(R,G,K);
    change = operate["change"];
    if(change == 1){
      G_new = Rcpp::as<arma::mat>(operate["G_new"]);
      prob_edge = operate["prob"];
      G_new=update_G(y,s,g,alpha,alpha_all,Q, G,R, N, J, K, L, p_add, pi0, delta0,G_new,prob_edge, 1);
    }
  } else {
    operate = reduce_path_sample1(G,K);
    change = operate["change"];
    if(change == 1){
      G_new = Rcpp::as<arma::mat>(operate["G_new"]);
      prob_edge = operate["prob"];
      G_new=update_G(y,s,g,alpha,alpha_all,Q, G,R, N, J, K, L, p_add, pi0, delta0,G_new,prob_edge, 2);
    }
  }

  return G_new;
}  

