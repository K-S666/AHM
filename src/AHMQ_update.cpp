#include "ahm_common.h"

using namespace Rcpp;

// [[Rcpp::export]]
void AHM_update(arma::mat Y,arma::mat& alpha,arma:: vec& s,arma:: vec& g,arma::mat& Q,arma::mat& G,arma::vec& pi,
                       double p_add,arma::mat alpha_all,arma::mat Q_all,int N,int J,int K,int L,int N1,
                       double a_s0,double a_g0,double b_s0,double b_g0,arma::vec delta0,arma::mat eta)
{
  arma::mat y(N1,J);
  //update alpha, pi
  Rcpp::List r_alpha;
  r_alpha=AHM_alpha_pi_update(Y,s,g,Q,G,alpha_all,pi,delta0,N,J,K,L);
  alpha=Rcpp::as<arma::mat>(r_alpha["alpha_new"]);
  pi=Rcpp::as<arma::vec>(r_alpha["pi_new"]);
  arma::vec alpha_current_possible_index=Rcpp::as<arma::vec>(r_alpha["alpha_current_possible_index"]);
  //sample mini-batch
  arma::vec x = generate_sequence(N)-1;
  arma::ivec sampled_values = sample_int(N,N1);
  arma::uvec sample_index = arma::conv_to<arma::uvec>::from(sampled_values);
  y=Y.rows(sample_index);
  arma::mat alpha_mini=alpha.rows(sample_index);
  //update s,g
  Rcpp::List rr;
  rr = AHM_sg_update(y, Q,alpha_mini, s,a_s0,b_s0,a_g0,b_g0);
  s = Rcpp::as<arma::vec>(rr["ss_new"]);
  g = Rcpp::as<arma::vec>(rr["gs_new"]);
  //update Q
  arma::vec alpha_code1=Trans_2to10_mat(alpha,K);
  arma::vec alpha_code=alpha_code1.elem(sample_index);
  arma::mat R = Reachability(G,K);
  Q = AHM_Q_update(y,s,g,alpha_mini,Q_all,Q,N1,J,K,alpha_current_possible_index,eta,alpha_code);
  //update G
  G = AHM_G_update(y,s,g,alpha_mini,alpha_all,Q,G,R,N1,J,K,L,p_add,pi,delta0);
  
}


// [[Rcpp::export]]
void AHM_update_fixQ(arma::mat Y,arma::mat& alpha,arma:: vec& s,arma:: vec& g,arma::mat& Q,arma::mat& G,arma::vec& pi,
                               double p_add,arma::mat alpha_all,arma::mat Q_all,int N,int J,int K,int L,int N1,
                               double a_s0,double a_g0,double b_s0,double b_g0,arma::vec delta0,arma::mat eta)
{
  arma::mat y(N1,J);
  //update alpha, pi
  Rcpp::List r_alpha;
  r_alpha=AHM_alpha_pi_update(Y,s,g,Q,G,alpha_all,pi,delta0,N,J,K,L);
  alpha=Rcpp::as<arma::mat>(r_alpha["alpha_new"]);
  pi=Rcpp::as<arma::vec>(r_alpha["pi_new"]);
  arma::vec alpha_current_possible_index=Rcpp::as<arma::vec>(r_alpha["alpha_current_possible_index"]);
  //sample mini-batch
  arma::vec x = generate_sequence(N)-1;
  arma::ivec sampled_values = sample_int(N,N1);
  arma::uvec sample_index = arma::conv_to<arma::uvec>::from(sampled_values);
  y=Y.rows(sample_index);
  arma::mat alpha_mini=alpha.rows(sample_index);
  //update s,g
  Rcpp::List rr;
  rr = AHM_sg_update(y, Q,alpha_mini, s,a_s0,b_s0,a_g0,b_g0);
  s = Rcpp::as<arma::vec>(rr["ss_new"]);
  g = Rcpp::as<arma::vec>(rr["gs_new"]);
  arma::mat R = Reachability(G,K);
  // Q is fixed in this variant; only G is sampled.
  G = AHM_G_update(y,s,g,alpha_mini,alpha_all,Q,G,R,N1,J,K,L,p_add,pi,delta0);
  
}

