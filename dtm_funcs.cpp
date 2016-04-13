#include "dtm_funcs.h"

List compute_e_log_beta(List m, List v, NumericMatrix zeta, List e_log_beta) {
  int T = e_log_beta.size();
  NumericMatrix temp_mat = as<NumericMatrix>(v[0]);
  int W = temp_mat.nrow();
  int K = temp_mat.ncol();
  
  for (int t = 0; t < T; t++) {
    NumericMatrix mat(W, K);
    NumericMatrix m_t = m[t];
    NumericMatrix v_t = v[t];
    NumericMatrix exp_sum(m_t.nrow(), m_t.ncol());
    for (int w = 0; w < W; w++) {
      for (int k = 0; k < K; k++) {
        exp_sum(w, k) = exp(m_t(w, k) + v_t(w, k) / 2.0);
      }
    }
    
    for (int k = 0; k < K; k++) {
      double exp_sum_k = sum(exp_sum(_, k));
      for (int w = 0; w < W; w++) {
        mat(w, k) = m_t(w, k) - (1.0 / zeta(t, k)) * exp_sum_k - log(zeta(t, k)) + 1.0;
      }
    }
    e_log_beta[t] = mat;
  }
  
  return e_log_beta;
}

NumericVector compute_e_log_theta(NumericVector gammas, NumericVector e_log_theta) {
  for (int k = 0; k < gammas.size(); k++) {
    e_log_theta[k] = digamma(gammas[k]); 
  }
  
  return e_log_theta;
}

void phi_update(NumericMatrix& phi, NumericVector e_log_theta, NumericMatrix e_log_beta, 
                NumericVector uniq_words) {
  // we do things in log space to deal with underflow
  int W = phi.nrow();
  int K = phi.ncol();
  for (int w = 0; w < W; w++) {
    double phi_topic_sum = 0.0;
    for (int k = 0; k < K; k++) {
      phi(w, k) = e_log_theta[k] + e_log_beta(uniq_words[w], k);  
      if (k > 0)
        phi_topic_sum = log_sum(phi_topic_sum, phi(w, k));
      else
        phi_topic_sum = phi(w, k);
    }
    for (int k = 0; k < K; k++) {
      phi(w, k) = phi(w, k) - phi_topic_sum;
    }
  }
}

double gamma_update(NumericVector& gamma_row, NumericMatrix phi, NumericVector n, double alpha,
                    NumericVector uniq_words) {
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
      gamma_row[k] += exp(phi(w, k)) * n[uniq_words[w]];
    }
  }
  for (int k = 0; k < K; k++) {
    gamma_change[k] = std::abs(gamma_row[k] - gamma_old[k]);
  }
  
  double avg_gamma_change = sum(gamma_change) / gamma_change.size();
  return avg_gamma_change;
}

void forward_backward(List& m, List& v, double sigma, double nu) {
  
}