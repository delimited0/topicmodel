//
// Online variational Bayes for LDA, as in Hoffman et al. 2010
// Not truly online
//

#include "lda_funcs.h"

// [[Rcpp::export]]
List lda_online(NumericMatrix dtm, int K, double alpha, double eta, double kappa, 
                double tau, double gam_tol, double doc_max_iter) {
  int D = dtm.nrow();
  int W = dtm.ncol();
  List phis(D);
  for (int d = 0; d < D; d++) {
    NumericMatrix mat(W, K);
    phis[d] = mat;
  }
  NumericMatrix gammas(D, K);
  NumericMatrix lambda(W, K);
  for (int k = 0; k < K; k++)
    lambda(_, k) = runif(W);
  NumericMatrix e_log_theta(D, K);
  NumericMatrix e_log_beta(W, K);
  
  e_log_beta = compute_e_log_beta(lambda, e_log_beta);
  
  for (int t = 0; t < D; t++) {
    double avg_gamma_change = 10.0;
    NumericVector gamma_row(gammas.ncol(), 1.0);
    e_log_theta(t, _) = compute_e_log_theta(gamma_row, e_log_theta(t, _));
    int doc_iter = 0;
    if (t % 100 == 0)
      std::cout << "Fitting Document: " << t << std::endl;
    
    while (avg_gamma_change > gam_tol && doc_iter < doc_max_iter) {
      phis[t] = phi_update(phis[t], e_log_theta(t, _), e_log_beta);
      List gamma_result = gamma_update(gamma_row, phis[t], dtm(t, _), alpha);
      NumericVector gamma_row = gamma_result["gamma"];
      gammas(t, _) = gamma_row; 
      e_log_theta(t, _) = compute_e_log_theta(gamma_row, e_log_theta(t, _));
      avg_gamma_change = gamma_result["avg_gamma_change"];
      doc_iter++;
    }
    
    NumericMatrix phi = as<NumericMatrix>(phis[t]);
    double rho = pow(tau + t, -kappa);
    for (int w = 0; w < W; w++) {
      for (int k = 0; k < K; k++) {
        lambda(w, k) = (1 - rho) * lambda(w, k) + rho * (eta + D * dtm(t, w) * phi(w, k));   
      }
    }
  }
    
  List result;
  result["gamma"] = gammas;
  result["phi"] = phis;
  result["lambda"] = lambda;
  result["e_log_theta"] = e_log_theta;
  result["e_log_beta"] = e_log_beta;
  return result;  
}
  
  