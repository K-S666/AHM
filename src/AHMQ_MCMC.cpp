#include "ahm_common.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List AHM_Q(arma::mat Y,int K,int N1=128,int chain_length=20000,int burn_in=10000,
                                    double a_s0=1.0,double a_g0=1.0,double b_s0=1.0,double b_g0=1.0,double p_add=0.5,
                                    int chain_id=1,bool progress=true,int print_every=1000)
{
  int N=Y.n_rows;int J=Y.n_cols;int L=pow(2,K);int D=L-1;
  int T = chain_length-burn_in;
  int progress_every = std::max(1, print_every);
  arma::cube QQ(J,K,T);
  arma::cube GG(K,K,T);
  arma::cube AA(N,K,T);
  arma::mat SLIP(J,T);arma::mat GUESS(J,T);
  arma::mat PIS(L,T);
  arma::mat alpha_all(L,K);
  arma::mat Q_all(L-1,K);
  arma::vec classvec=arma::linspace<arma::vec>(0, L-1, L);
  alpha_all=Trans_10to2_mat(K,classvec);
  Q_all=alpha_all.rows(1,L-1);
  // prior
  arma::vec delta0=arma::ones(L,1);
  arma::vec pi0=1.0/L*arma::ones(L,1);
  //initial value
  arma::vec s=runif(J)*0.5;
  arma::vec g=runif(J)*0.5;
  arma::vec pi=runif(L);pi=pi/sum(pi);
  arma::mat G=arma::zeros(K,K);
  arma::mat Q=random_Q(J,K);
  arma::mat alpha(N,K,arma::fill::zeros);
  arma::mat eta(L,D);
  int tmburn;
  for(int l=0;l<L;l++){
    for(int d=0;d<D;d++){
      eta(l,d)=dina_eta_row(alpha_all.row(l), Q_all.row(d));
    }
  }
  for(int t = 0; t < chain_length; t++){
    //update
    AHM_update(Y,alpha,s,g,Q,G,pi,p_add,alpha_all,Q_all,N,J,K,L,N1,a_s0,a_g0,b_s0,b_g0,delta0,eta);    
    if(t>burn_in-1){
      tmburn = t-burn_in;
      SLIP.col(tmburn)  = s;
      GUESS.col(tmburn) = g;
      AA.slice(tmburn)  = alpha;
      GG.slice(tmburn)  = G;
      QQ.slice(tmburn)  = Q;
      PIS.col(tmburn)   = pi;
      
    }
    if (t % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (progress && ((t + 1) == 1 || ((t + 1) % progress_every == 0) || ((t + 1) == chain_length))) {
      Rcpp::Rcout << "\rAHMQ chain " << chain_id << " iteration "
                  << (t + 1) << "/" << chain_length;
      if ((t + 1) == chain_length) Rcpp::Rcout << "\n";
      Rcpp::Rcout.flush();
    }
    
  }
  
  
  
  return Rcpp::List::create(Rcpp::Named("alpha") = AA,Rcpp::Named("s") = SLIP,Rcpp::Named("g") = GUESS,
                            Rcpp::Named("pi") = PIS,Rcpp::Named("Q") = QQ,Rcpp::Named("G") = GG);
}



