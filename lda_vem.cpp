/**
 * Fully Bayesian LDA, fit using variational inference.
 * 
 */
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
List lda_vem(NumericMatrix dtm, int K, double alpha, double gam_tol, double em_tol,
             double em_max_iter) {
  int D = dtm.nrow();
  int W = dtm.ncol();
  List phis(D);
  NumericMatrix gammas(D, K);
  NumericMatrix lambda(K, W);
  NumericMatrix e_log_theta(D, K);
  NumericMatrix e_log_beta(K, W);
  
  for (int k = 0; k < lambda.nrow(); k++)
    lambda(k, _) = runif(W);
  
  double lold = 0.0;
  double lnew = 10000.0;
  int iter = 0;
  
  while ((std::abs(lnew - lold) / lold) > em_tol || iter < em_max_iter) {
    lold = lnew;
    lnew = 0.0;
    std::fill(gammas.begin(), gammas.end(), 1);
    compute_e_log_beta(lambda, e_log_beta);
    compute_e_log_theta(gammas, e_log_theta);
    for (int d = 0; d < D; d++) {
      while () {
        // phis[d] = phi_update
        // gammas(d, _) = gamma_update
        lnew += log_lik(dtm(d, _), phis[d], e_log_theta(d, _), e_log_beta,
                        lambda);
      } 
    }
    
    // lambda = lambda_update()
    lnew = log_lik()
  }
}

void compute_e_log_beta(NumericMatrix lambda, NumericMatrix e_log_beta) {
  for (int k = 0; k < lambda.ncol(); k++) {
    double lambda_sum = sum(lambda(_, k));
    double dig_sum = boost::math::digamma(lambda_sum);
    for (int w = 0; w < lambda.ncol(); w++) {
      e_log_beta(w, k) = boost::math::digamma(lambda(w, k)) - dig_sum;
    }
  }
}

// can be shared with dtm
void compute_e_log_theta(NumericMatrix gammas, NumericMatrix e_log_theta) {
  for (int d = 0; d < gammas.nrow(); d++) {
    double gamma_sum = sum(gammas(d, _));
    double dig_sum = boost::math::digamma(gamma_sum);
    for (int k = 0; k < gammas.ncol(); k++) {
      e_log_theta(d, k) = boost::math::digamma(gammas(d, k)) - dig_sum;
    }
  }
}

double log_lik(NumericVector n, NumericMatrix phi, NumericVector e_log_theta, 
               NumericMatrix e_log_beta, NumericMatrix lambda) {
  
}

