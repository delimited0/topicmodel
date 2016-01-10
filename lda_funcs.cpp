#include "lda_funcs.h"
#include "util.h"

NumericMatrix compute_e_log_beta(NumericMatrix lambda, NumericMatrix e_log_beta) {
  NumericVector lambda_dig_sums(lambda.ncol());
  for (int k = 0; k < lambda.ncol(); k++) {
    lambda_dig_sums[k] = sum(lambda(_, k));  
    lambda_dig_sums[k] = digamma(lambda_dig_sums[k]);
  }
  
  for (int w = 0; w < lambda.nrow(); w++) {
    for (int k = 0; k < lambda.ncol(); k++) {
      e_log_beta(w, k) = digamma(lambda(w, k)) - lambda_dig_sums[k];
    }
  }
  
  return e_log_beta;
}

NumericVector compute_e_log_theta(NumericVector gammas, NumericVector e_log_theta) {
  double gamma_sum = sum(gammas);
  double dig_sum = digamma(gamma_sum);
  for (int k = 0; k < gammas.size(); k++) {
    e_log_theta[k] = digamma(gammas[k]) - dig_sum;
  }
  
  return e_log_theta;
}

NumericMatrix phi_update(NumericMatrix phi, NumericVector e_log_theta, NumericMatrix e_log_beta) {
  int W = phi.nrow();
  int K = phi.ncol();
  for (int w = 0; w < W; w++) {
    double phi_topic_sum = 0.0;
    for (int k = 0; k < K; k++) {
      phi(w, k) = e_log_theta[k] + e_log_beta(w, k);  
      if (k > 0)
        phi_topic_sum = log_sum(phi_topic_sum, phi(w, k));
      else
        phi_topic_sum = phi(w, k);
    }
    for (int k = 0; k < K; k++) {
      phi(w, k) = phi(w, k) - phi_topic_sum;
    }
  }
  
  return phi;
}

List gamma_update(NumericVector gamma_row, NumericMatrix phi, NumericVector n, double alpha) {
  int W = phi.nrow();
  int K = phi.ncol();
  NumericVector gamma_change(gamma_row.size());
  NumericVector gamma_old(gamma_row.size());
  for (int k = 0; k < K; k++) {
    gamma_old[k] = gamma_row[k];
  }
  std::fill(gamma_row.begin(), gamma_row.end(), alpha);
  for (int w = 0; w < W; w++) {
    for (int k = 0; k < K; k++) {
      gamma_row[k] += exp(phi(w, k)) * n[w];
    }
  }
  for (int k = 0; k < K; k++) {
    gamma_change[k] = std::abs(gamma_row[k] - gamma_old[k]);
  }
  
  double avg_gamma_change = sum(gamma_change) / gamma_change.size();
  List result;
  result["avg_gamma_change"] = avg_gamma_change;
  result["gamma"] = gamma_row;
  return result;
}