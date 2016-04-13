//
// Fully Bayesian LDA, fit using variational inference.
//
//

#include "lda_funcs.h"
#include "util.h"

NumericMatrix lambda_update(NumericMatrix lambda, NumericMatrix dtm, List phis, double eta,
                            List doc_words) {
  std::fill(lambda.begin(), lambda.end(), eta);
  for (int d = 0; d < phis.size(); d++) {
    NumericMatrix phi = as<NumericMatrix>(phis[d]);
    NumericVector uniq_words = as<NumericVector>(doc_words[d]);
    for (int w = 0; w < phi.nrow(); w++) {
      for (int k = 0; k < lambda.ncol(); k++) {
        lambda(uniq_words[w], k) += dtm(d, uniq_words[w]) * exp(phi(w, k));
      }
    }
  }
  
  return lambda;
}

double log_lik(NumericVector n, NumericMatrix phi, NumericVector gamma_row, NumericVector e_log_theta, 
               NumericMatrix e_log_beta, NumericMatrix lambda, double alpha, double eta, int D,
               NumericVector uniq_words) {
  int K = e_log_beta.ncol();
  int W = uniq_words.size();
  
  double term1 = 0.0;
  for (int w = 0; w < W; w++) {
    double acc = 0.0;
    for (int k = 0; k < K; k++) {
      acc += exp(phi(w, k)) * (e_log_theta[k] + e_log_beta(uniq_words[w], k) - phi(w, k));
    }
    term1 += n[uniq_words[w]] * acc;
  }
  
  double term2 = 0.0;
  double neglgam = 0.0;
  for (int k = 0; k < K; k++) {
    term2 += (alpha - gamma_row[k]) * e_log_theta[k] + lgamma(gamma_row[k]);
    neglgam += gamma_row[k];
  }
  neglgam = - lgamma(neglgam);
  term2 += neglgam;
  
  double term3 = 0.0;
  for (int w = 0; w < W; w++) {
    for (int k = 0; k < K; k++) {
      term3 += (eta - lambda(uniq_words[w], k)) * e_log_beta(uniq_words[w], k) + 
        lgamma(lambda(uniq_words[w], k));
    }
  }
  for (int k = 0; k < K; k++) {
    double lamsum = sum(lambda(_, k));
    term3 -= lgamma(lamsum);
  }
  term3 = term3 / D;
  
  double term4 = lgamma(K * alpha) - K * lgamma(alpha) + 
    (lgamma(W * eta) - W * lgamma(eta)) / D;
  
  return term1 + term2 + term3 + term4;
}

// [[Rcpp::export]]
List lda_vb(NumericMatrix dtm, int K, double alpha, double eta, double gam_tol, double em_tol,
             double em_max_iter, double doc_max_iter) {
  int D = dtm.nrow();
  int W = dtm.ncol();
  List phis(D);
  List doc_words(D);
  for (int d = 0; d < D; d++) {
    int nonzero = 0;
    for (int w = 0; w < W; w++) {
      if (dtm(d, w) > 0)
        nonzero++;
    }
    NumericMatrix mat(nonzero, K);
    IntegerVector uniq_words = whichPositive(dtm(d, _));
    std::fill(mat.begin(), mat.end(), 1.0 / K);
    phis[d] = mat;
    doc_words[d] = uniq_words;
  }
  NumericMatrix gammas(D, K);
  std::fill(gammas.begin(), gammas.end(), 1.0);
  NumericMatrix lambda(W, K);
  NumericMatrix e_log_theta(D, K);
  NumericMatrix e_log_beta(W, K);
  NumericVector log_liks(em_max_iter);
  
  for (int k = 0; k < K; k++)
    lambda(_, k) = runif(W);
  
  double lold = 0;
  double lnew = 0;
  int iter = 0;
  
  e_log_beta = compute_e_log_beta(lambda, e_log_beta);
  for (int d = 0; d < D; d++) {
    e_log_theta(d, _) = compute_e_log_theta(gammas(d, _), e_log_theta(d, _));
  }
  
  while (((std::abs((lnew - lold) / lold)) > em_tol && iter < em_max_iter) || iter == 0) {
    Rcout << "Fractional change in lhood: " << std::abs((lnew - lold) / lold) << std::endl;
    Rcout << "Iteration: " << iter << std::endl;
    lold = lnew;
    lnew = 0.0;
    
    for (int d = 0; d < D; d++) {
      double avg_gamma_change = 10.0;
      NumericMatrix phi = as<NumericMatrix>(phis[d]);
      NumericVector gamma_row(gammas.ncol(), 1.0);
      e_log_theta(d, _) = compute_e_log_theta(gamma_row, e_log_theta(d, _));
      int doc_iter = 0;
      if (d % 100 == 0)
        Rcout << "Fitting Document: " << d << std::endl;
      while (avg_gamma_change > gam_tol && doc_iter < doc_max_iter) {
        phi_update(phi, e_log_theta(d, _), e_log_beta, doc_words[d]);
        avg_gamma_change = gamma_update(gamma_row, phis[d], dtm(d, _), alpha, doc_words[d]);
        gammas(d, _) = gamma_row; 
        e_log_theta(d, _) = compute_e_log_theta(gamma_row, e_log_theta(d, _));
        doc_iter++;
      }
    }
    
    lambda = lambda_update(lambda, dtm, phis, eta, doc_words);
    e_log_beta = compute_e_log_beta(lambda, e_log_beta);
    
    for (int d = 0; d < D; d++) {
      lnew += log_lik(dtm(d, _), phis[d], gammas(d, _), e_log_theta(d, _), e_log_beta,
                      lambda, alpha, eta, D, doc_words[d]);
    }
    log_liks[iter] = lnew;
    iter++;
    Rcout << "Log Likelihood: " << lnew << std::endl;
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
