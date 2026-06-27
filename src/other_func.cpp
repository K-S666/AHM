#include "ahm_common.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec inv_bijectionvector(unsigned int K,double CL){
  arma::vec alpha(K);
  for(unsigned int k=0;k<K;k++){
    double twopow = pow(2,K-k-1);
    alpha(k) = (twopow<=CL);
    CL = CL - twopow*alpha(k);
  }
  return alpha;
}

// [[Rcpp::export]]
double twoToten(arma::ivec x)
{
  int K=x.n_elem;
  double res=0;
  for(int i=0;i<K;i++){
    res+=x(K-i-1)*pow(2,i);
  }
  return res;
}
// [[Rcpp::export]]
arma::mat ff(arma::mat a){
  arma::mat temp;
  temp=a;temp.diag().zeros();
  return temp;
}
double ahm_clamp_prob(double p)
{
  const double eps = 1e-12;
  if (p < eps) return eps;
  if (p > 1.0 - eps) return 1.0 - eps;
  return p;
}

double ahm_log_sum_exp(const arma::vec& log_values)
{
  double m = log_values.max();
  if (!std::isfinite(m)) {
    return m;
  }
  return m + std::log(arma::sum(arma::exp(log_values - m)));
}

arma::vec ahm_softmax_from_log(const arma::vec& log_values)
{
  double log_total = ahm_log_sum_exp(log_values);
  if (!std::isfinite(log_total)) {
    return arma::ones<arma::vec>(log_values.n_elem) / log_values.n_elem;
  }
  arma::vec probs = arma::exp(log_values - log_total);
  return probs / arma::sum(probs);
}

double dina_eta_row(const arma::rowvec& alpha, const arma::rowvec& q)
{
  for (arma::uword k = 0; k < q.n_elem; ++k) {
    if (q(k) > 0.5 && alpha(k) < 0.5) {
      return 0.0;
    }
  }
  return 1.0;
}

double dina_log_response_prob(double y, double eta, double s, double g)
{
  // DINA response probability:
  // eta = 1 -> P(Y=1)=1-s, P(Y=0)=s;
  // eta = 0 -> P(Y=1)=g,   P(Y=0)=1-g.
  double p = (eta > 0.5)
    ? (y > 0.5 ? 1.0 - s : s)
    : (y > 0.5 ? g : 1.0 - g);
  return std::log(ahm_clamp_prob(p));
}

double dina_log_response_vector(const arma::vec& eta, const arma::rowvec& y,
                                const arma::vec& s, const arma::vec& g)
{
  double log_lik = 0.0;
  for (arma::uword j = 0; j < eta.n_elem; ++j) {
    log_lik += dina_log_response_prob(y(j), eta(j), s(j), g(j));
  }
  return log_lik;
}

// [[Rcpp::export]]
arma::mat Boolean(arma::mat A){
  arma::mat R;
  R=arma::conv_to<arma::mat>::from(A>0);
  return R;
}

// [[Rcpp::export]]
arma::vec Booleanvec(arma::vec A){
  arma::vec R;
  R=arma::conv_to<arma::vec>::from(A>0);
  return R;
}

// [[Rcpp::export]]
arma::mat Reachability(const arma::mat& StrucMat,unsigned int K){
  arma::mat R=StrucMat;
  R.diag().zeros();
  arma::mat Identity=arma::eye(K,K);
  arma::mat Rnext=R+Identity;
  while(accu(R==Rnext)<(K*K)){
    R=Rnext;
    Rnext=Boolean(R*(R+Identity));
  }
  R.diag().ones();
  return R;
}
// [[Rcpp::export]]
arma::mat ConnectMat(const arma::mat& R,unsigned int K){
  arma::mat Connect(K,K);Connect.eye();
  for(unsigned int i=0;i<(K-1);i++){
    for(unsigned int j=(i+1);j<K;j++){
      if((R(i,j)==1)||(R(j,i)==1)){
        Connect(i,j)=1;Connect(j,i)=1;
      }
    }
  }
  return Connect;
}

// [[Rcpp::export]]
arma::mat  Transitive(const arma::mat& G,unsigned int K){
  arma::mat G_red=G;
  arma::mat R_c=Reachability(G_red,K);
  arma::mat C_c=ConnectMat(R_c,K);
  for(unsigned int i=0;i<K;i++){
    for(unsigned int j=0;j<K;j++){
      if(G(i,j)==1){
        arma::mat Gtemp=G_red;Gtemp(i,j)=0;
        arma::mat Rtemp=Reachability(Gtemp,K);
        arma::mat Ctemp=ConnectMat(Rtemp,K);
        if(arma::all(arma::vectorise(Ctemp==C_c))){
          G_red=Gtemp;
        }
      }
    }
  }
  return(G_red);
}

// [[Rcpp::export]]
arma::vec Trans_10to2(unsigned int K,double CL){
  arma::vec alpha(K);
  for(unsigned int k=0;k<K;k++){
    double twopow = pow(2,K-k-1);
    alpha(k) = (twopow<=CL);
    CL = CL - twopow*alpha(k);
  }
  return alpha;
}

// [[Rcpp::export]]
arma::mat Trans_10to2_mat(unsigned int K,const arma::vec& CL) {
  unsigned int Col=CL.n_elem;
  arma::mat alpha(Col,K);
  for(unsigned int i=0;i<Col;i++){
    arma::colvec alphai(K);double cl=CL(i);
    for(unsigned int k=0;k<K;k++){
      double twopow = pow(2,K-k-1);
      alphai(k) = (twopow<=cl);
      cl = cl - twopow*alphai(k);
    }
    alpha.row(i)=alphai.t();
  }
  return alpha;
}


// [[Rcpp::export]]
double Trans_2to10(arma::vec x,int K){
  double r=0.0;
  for(int k=0;k<K;k++){
    r=r+pow(2,k)*x(K-k-1);
  }
  return r;
}

// [[Rcpp::export]]
arma::vec Trans_2to10_mat(arma::mat x,int K){
  int n=x.n_rows;
  arma::vec code(n);
  for(int i=0;i<n;i++){
    code(i)=Trans_2to10((x.row(i)).t(),K);
  }
  return code;
}
// [[Rcpp::export]]
Rcpp::List Reduced_alpha(arma::mat alpha_all,arma::mat R,int K)//all possible alpha
{
  int L= pow(2,K);
  arma::mat alpha_current = alpha_all;
  arma::vec delta_new = arma::zeros(L,1);
  arma::uvec prerequisite;
  arma::uvec if_prerequisite;
  arma::mat alpha_prerequisite;
  arma::vec aa=arma::ones(L,1);
  for(int k=0;k<K;k++){
    aa=arma::ones(L,1);
    prerequisite = arma::find(R.col(k)==1);
    alpha_prerequisite=alpha_current.cols(prerequisite);
    if((prerequisite).n_elem > 0){
      if_prerequisite = arma::find(arma::sum(alpha_prerequisite,1) < (prerequisite).n_elem);
      aa.elem(if_prerequisite) = arma::zeros(if_prerequisite.n_elem,1);
      alpha_current.col(k)=aa;
    }
  }
  arma::vec alpha_current_binary;
  arma::vec index;
  alpha_current_binary = Trans_2to10_mat(alpha_current,K);
  index=sort(unique(alpha_current_binary));
  
  return Rcpp::List::create(Rcpp::Named("alpha_current") = alpha_current,
                            Rcpp::Named("alpha_current_binary") = alpha_current_binary,
                            Rcpp::Named("index") = index);
}


// [[Rcpp::export]]
arma::vec generate_sequence(int N) {
  arma::vec sequence = arma::linspace<arma::vec>(1, N, N);
  return sequence;
}

// [[Rcpp::export]]
Rcpp::IntegerVector sample_int(int N, int N1) {
  Rcpp::IntegerVector pool = Rcpp::seq(0, N - 1);
  return Rcpp::sample(pool, N1, false);
}


// [[Rcpp::export]]
arma::mat random_Q(unsigned int J,unsigned int K){
  if (J < K) {
    Rcpp::stop("random_Q requires J >= K so every attribute can be represented.");
  }
  
  //Generate identity matrices
  arma::vec one_K = arma::ones<arma::vec>(K);
  arma::mat I_K = arma::diagmat(one_K);
  //arma::mat Two_I_K = arma::join_cols(I_K,I_K);
  
  //generate Q1
  unsigned int JmK = J-K;
  if (JmK == 0) {
    return I_K.rows(arma::shuffle(arma::regspace<arma::uvec>(0, J - 1)));
  }
  unsigned int J1max = K;
  if(K>JmK){
    J1max = JmK;
  }
  unsigned int J1 = arma::randi<arma::uvec>(1, arma::distr_param(1, J1max))(0);
  arma::mat U1 = arma::randu<arma::mat>(J1,K);
  arma::mat Q1 = arma::zeros<arma::mat>(J1,K);
  
  //fix elements so columns are nonzero
  arma::vec row_ks = arma::randi<arma::vec>(K,arma::distr_param(0,J1-1) );
  for(unsigned int k=0;k<K;k++){
    Q1(row_ks(k),k) = 1;
  }
  
  Q1.elem(arma::find(Q1 > .5) ).fill(1.0);
  
  arma::mat Q = arma::join_cols(I_K,Q1);
  
  //Generating the remaining elements of Q in Q2 
  unsigned int JmKmJ1 = JmK - J1;
  arma::mat Q2 = arma::zeros<arma::mat>(JmKmJ1,K);
  if(JmKmJ1>0){
    arma::mat U2 = arma::randu<arma::mat>(JmKmJ1,K);
    Q2.elem(arma::find(U2 > .5) ).fill(1.0);
    Q = arma::join_cols(Q,Q2);
  }
  
  //Q
  arma::uvec P = arma::uvec(J);
  for(unsigned int j=0;j<J;j++){
    P(j)=j;
  }
  P = arma::shuffle(P);
  return Q.rows(P);
}

