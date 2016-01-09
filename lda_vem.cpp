//
// Fully Bayesian LDA, fit using variational inference.
//
//
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include "util.h"
#include <boost/math/special_functions/digamma.hpp>
using namespace Rcpp;

void print_vector(NumericVector x) {
  for (int i = 0; i < x.size(); i++)
    std::cout << x[i] << ", ";
  std::cout << std::endl;
}

// [[Rcpp::export]]
NumericMatrix compute_e_log_beta(NumericMatrix lambda, NumericMatrix e_log_beta) {
  NumericVector lambda_dig_sums(lambda.ncol());
  for (int k = 0; k < lambda.ncol(); k++) {
    lambda_dig_sums[k] = sum(lambda(_, k));  
    lambda_dig_sums[k] = digamma(lambda_dig_sums[k]);
  }
  
  for (int w = 0; w < lambda.nrow(); w++) {
    for (int k = 0; k < lambda.ncol(); k++) {
      //std::cout << "Topic: " << k << std::endl;
      e_log_beta(w, k) = digamma(lambda(w, k)) - lambda_dig_sums[k];
      //std::cout << "Word: " << w << std::endl;
    }
  }
  
  return e_log_beta;
}

// can be shared with dtm
// [[Rcpp::export]]
NumericVector compute_e_log_theta(NumericVector gammas, NumericVector e_log_theta) {
  double gamma_sum = sum(gammas);
  double dig_sum = digamma(gamma_sum);
  for (int k = 0; k < gammas.size(); k++) {
    e_log_theta[k] = digamma(gammas[k]) - dig_sum;
  }
  
  return e_log_theta;
}

// can be shared with dtm
// [[Rcpp::export]]
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

// can be shared with dtm
// [[Rcpp::export]]
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

// [[Rcpp::export]]
NumericMatrix lambda_update(NumericMatrix lambda, NumericMatrix dtm, List phis, double eta) {
  std::fill(lambda.begin(), lambda.end(), eta);
  for (int d = 0; d < phis.size(); d++) {
    NumericMatrix phi = as<NumericMatrix>(phis[d]);
    for (int w = 0; w < lambda.nrow(); w++) {
      for (int k = 0; k < lambda.ncol(); k++) {
        lambda(w, k) += dtm(d, w) * exp(phi(w, k));
      }
    }
  }
  
  return lambda;
}

double log_lik(NumericVector n, NumericMatrix phi, NumericVector gamma_row, NumericVector e_log_theta, 
               NumericMatrix e_log_beta, NumericMatrix lambda, double alpha, double eta, int D) {
  int K = e_log_beta.ncol();
  int W = e_log_beta.nrow();
  
  double term1 = 0.0;
  for (int w = 0; w < W; w++) {
    double acc = 0.0;
    for (int k = 0; k < K; k++) {
      acc += exp(phi(w, k)) * (e_log_theta[k] + e_log_beta(w, k) - phi(w, k));
    }
    term1 += n[w] * acc;
  }
  // std::cout << "term1: " << term1 << std::endl;
  
  double term2 = 0.0;
  double neglgam = 0.0;
  for (int k = 0; k < K; k++) {
    term2 += (alpha - gamma_row[k]) * e_log_theta[k] + lgamma(gamma_row[k]);
    neglgam += gamma_row[k];
  }
  neglgam = - lgamma(neglgam);
  term2 += neglgam;
  // std::cout << "term2: " << term2 << std::endl;
  
  double term3 = 0.0;
  for (int w = 0; w < W; w++) {
    for (int k = 0; k < K; k++) {
      term3 += (eta - lambda(w, k)) * e_log_beta(w, k) + lgamma(lambda(w, k));
    }
  }
  for (int k = 0; k < K; k++) {
    double lamsum = sum(lambda(_, k));
    term3 -= lgamma(lamsum);
  }
  term3 = term3 / D;
  // std::cout << "term3: " << term3 << std::endl;
  
  double term4 = lgamma(K * alpha) - K * lgamma(alpha) + 
    (lgamma(W * eta) - W * lgamma(eta)) / D;
  //std::cout << "term4: " << term4 << std::endl;
  
  return term1 + term2 + term3 + term4;
}

// [[Rcpp::export]]
List lda_vem(NumericMatrix dtm, int K, double alpha, double eta, double gam_tol, double em_tol,
             double em_max_iter, double doc_max_iter) {
  int D = dtm.nrow();
  int W = dtm.ncol();
  List phis(D);
  for (int d = 0; d < D; d++) {
    // std::cout << "Filling Document: " << d << std::endl;
    NumericMatrix mat(W, K);
    std::fill(mat.begin(), mat.end(), 1.0 / K);
    phis[d] = mat;
  }
  NumericMatrix gammas(D, K);
  std::fill(gammas.begin(), gammas.end(), 1.0);
  NumericMatrix lambda(W, K);
  NumericMatrix e_log_theta(D, K);
  NumericMatrix e_log_beta(W, K);
  NumericVector log_liks(em_max_iter);
  
  for (int k = 0; k < K; k++)
    lambda(_, k) = runif(W);
  
  double lold = 1.0;
  double lnew = 10.0;
  int iter = 0;
  e_log_beta = compute_e_log_beta(lambda, e_log_beta);
  for (int d = 0; d < D; d++) {
    e_log_theta(d, _) = compute_e_log_theta(gammas(d, _), e_log_theta(d, _));
  }
  
  while ((std::abs((lnew - lold) / lold)) > em_tol && iter < em_max_iter) {
    std::cout << "Fractional change in lhood: " << std::abs((lnew - lold) / lold) << std::endl;
    std::cout << "Iteration: " << iter << std::endl;
    lold = lnew;
    lnew = 0.0;
    
    for (int d = 0; d < D; d++) {
      double avg_gamma_change = 10.0;
      NumericMatrix phi = as<NumericMatrix>(phis[d]);
      NumericVector gamma_row(gammas.ncol(), 1.0);
      e_log_theta(d, _) = compute_e_log_theta(gamma_row, e_log_theta(d, _));
      int doc_iter = 0;
      if (d % 100 == 0)
        std::cout << "Fitting Document: " << d << std::endl;
      while (avg_gamma_change > gam_tol && doc_iter < doc_max_iter) {
        phis[d] = phi_update(phi, e_log_theta(d, _), e_log_beta);
        List gamma_result = gamma_update(gamma_row, phis[d], dtm(d, _), alpha);
        NumericVector gamma_row = gamma_result["gamma"];
        gammas(d, _) = gamma_row; 
        e_log_theta(d, _) = compute_e_log_theta(gamma_row, e_log_theta(d, _));
        avg_gamma_change = gamma_result["avg_gamma_change"];
        doc_iter++;
      }
      lnew += log_lik(dtm(d, _), phis[d], gammas(d, _), e_log_theta(d, _), e_log_beta,
                      lambda, alpha, eta, D);
    }
    
    lambda = lambda_update(lambda, dtm, phis, eta);
    e_log_beta = compute_e_log_beta(lambda, e_log_beta);
    
    log_liks[iter] = lnew;
    iter++;
    std::cout << "Log Likelihood: " << lnew << std::endl;
  }
  
  List result;
  result["gamma"] = gammas;
  result["phi"] = phis;
  result["lambda"] = lambda;
  result["log_liks"] = log_liks;
  result["e_log_theta"] = e_log_theta;
  result["e_log_beta"] = e_log_beta;
  result.attr("class") = "lda";
  return result;
}
