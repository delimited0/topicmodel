/**
 * Fully Bayesian LDA, fit using variational inference.
 * 
 */
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
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
  
  // lambda = random()
  
  double lold = 0.0;
  double lnew = 10000.0;
  int iter = 0;
  
  while ((std::abs(lnew - lold) / lold) > em_tol || iter < em_max_iter) {
    std::fill(gammas.begin(), gammas.end(), 1);
    // e_log_beta = 
    for (int d = 0; d < D; d++) {
      while () {
        // e_log_theta(d, _) = 
        // phis[d] = phi_update
        // gammas(d, _) = gamma_update
      }
    }
    
    // lambda = lambda_update()
  }
}

double log_lik(NumericVector n, NumericMatrix phi, NumericVector e_log_theta, 
               NumericMatrix e_log_beta, NumericMatrix lambda) {
  
}

