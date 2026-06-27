#ifndef AHM_COMMON_H
#define AHM_COMMON_H

#include <RcppArmadillo.h>
#include <rgen.h>
#include <algorithm>
#include <R_ext/Utils.h>
#include <iomanip>

arma::vec inv_bijectionvector(unsigned int K, double CL);
double twoToten(arma::ivec x);
arma::mat ff(arma::mat a);
double ahm_log_sum_exp(const arma::vec& log_values);
arma::vec ahm_softmax_from_log(const arma::vec& log_values);
double dina_eta_row(const arma::rowvec& alpha, const arma::rowvec& q);
double dina_log_response_prob(double y, double eta, double s, double g);
double dina_log_response_vector(const arma::vec& eta, const arma::rowvec& y,
                                const arma::vec& s, const arma::vec& g);
arma::mat Boolean(arma::mat A);
arma::vec Booleanvec(arma::vec A);
arma::mat Reachability(const arma::mat& StrucMat, unsigned int K);
arma::mat ConnectMat(const arma::mat& R, unsigned int K);
arma::mat Transitive(const arma::mat& G, unsigned int K);
arma::vec Trans_10to2(unsigned int K, double CL);
arma::mat Trans_10to2_mat(unsigned int K, const arma::vec& CL);
double Trans_2to10(arma::vec x, int K);
arma::vec Trans_2to10_mat(arma::mat x, int K);
Rcpp::List Reduced_alpha(arma::mat alpha_all, arma::mat R, int K);
arma::vec rDirichlet(const arma::vec& deltas);
arma::vec generate_sequence(int N);
Rcpp::IntegerVector sample_int(int N, int N1);
arma::mat random_Q(unsigned int J, unsigned int K);

Rcpp::List AHM_sg_update(const arma::mat& Y, const arma::mat& Q,
                         const arma::mat& ALPHAS, const arma::vec& ss_old,
                         double as0, double bs0, double ag0, double bg0);

arma::mat AHM_alpha_update(arma::mat y, arma::vec s, arma::vec g,
                           arma::mat Q, arma::mat G, arma::mat alpha_all,
                           arma::vec pi0, int N, int J, int K, int L);
Rcpp::List AHM_alpha_pi_update(arma::mat y, arma::vec s, arma::vec g,
                               arma::mat Q, arma::mat G,
                               arma::mat alpha_all, arma::vec pi0,
                               arma::vec delta0, int N, int J, int K,
                               int L);

arma::mat AHM_Q_update(arma::mat y, arma::vec s, arma::vec g,
                       arma::mat alpha, arma::mat Q_all, arma::mat Q,
                       int N, int J, int K,
                       arma::vec alpha_current_possible_index,
                       arma::mat eta, arma::vec alpha_code);

Rcpp::List add_path_sample1(arma::mat R, arma::mat G, int K);
Rcpp::List reduce_path_sample1(arma::mat G, int K);
arma::mat update_G(arma::mat y, arma::vec s, arma::vec g,
                   arma::mat alpha, arma::mat alpha_all,
                   arma::mat Q, arma::mat G, arma::mat R,
                   int N, int J, int K, int L,
                   double p_add, arma::vec pi,
                   arma::vec delta0, arma::mat G_new,
                   double prob_edge, int action);
arma::mat AHM_G_update(arma::mat y, arma::vec s, arma::vec g,
                       arma::mat alpha, arma::mat alpha_all,
                       arma::mat Q, arma::mat G, arma::mat R,
                       int N, int J, int K, int L,
                       double p_add, arma::vec pi0,
                       arma::vec delta0);

void AHM_update(arma::mat Y, arma::mat& alpha,
                arma::vec& s, arma::vec& g,
                arma::mat& Q, arma::mat& G, arma::vec& pi,
                double p_add, arma::mat alpha_all,
                arma::mat Q_all, int N, int J, int K, int L,
                int N1, double a_s0, double a_g0,
                double b_s0, double b_g0, arma::vec delta0,
                arma::mat eta);
void AHM_update_fixQ(arma::mat Y, arma::mat& alpha,
                     arma::vec& s, arma::vec& g,
                     arma::mat& Q, arma::mat& G,
                     arma::vec& pi, double p_add,
                     arma::mat alpha_all, arma::mat Q_all,
                     int N, int J, int K, int L, int N1,
                     double a_s0, double a_g0,
                     double b_s0, double b_g0,
                     arma::vec delta0, arma::mat eta);

#endif

