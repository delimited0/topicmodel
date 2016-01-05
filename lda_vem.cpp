//
// Fully Bayesian LDA, fit using variational inference.
//
//
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
using namespace Rcpp;

void print_vector(NumericVector x) {
  for (int i = 0; i < x.size(); i++)
    std::cout << x[i] << ", ";
  std::cout << std::endl;
}

double log_sum(double log_a, double log_b) {
  double v;
  
  if (log_a < log_b) {
    v = log_b+log(1 + exp(log_a-log_b));
  }
  else {
    v = log_a+log(1 + exp(log_b-log_a));
  }
  return(v);
}

// [[Rcpp::export]]
NumericMatrix compute_e_log_beta(NumericMatrix lambda, NumericMatrix e_log_beta) {
  for (int k = 0; k < lambda.nrow(); k++) {
    double lambda_sum = sum(lambda(k, _));
    double dig_sum = boost::math::digamma(lambda_sum);
    //std::cout << "Topic: " << k << std::endl;
    for (int w = 0; w < lambda.ncol(); w++) {
      e_log_beta(k, w) = boost::math::digamma(lambda(k, w)) - dig_sum;
      //std::cout << "Word: " << w << std::endl;
    }
  }
  
  return e_log_beta;
}

// can be shared with dtm
// [[Rcpp::export]]
NumericVector compute_e_log_theta(NumericVector gammas, NumericVector e_log_theta) {
  double gamma_sum = sum(gammas);
  //std::cout << "Document: " << d << " gamma sum: " << gamma_sum << std::endl;
  double dig_sum = boost::math::digamma(gamma_sum);
  for (int k = 0; k < gammas.size(); k++) {
    e_log_theta[k] = boost::math::digamma(gammas[k]) - dig_sum;
  }
  
  return e_log_theta;
}

// can be shared with dtm
// [[Rcpp::export]]
NumericMatrix phi_update(NumericMatrix phi, NumericVector e_log_theta, NumericMatrix e_log_beta) {
  for (int w = 0; w < phi.ncol(); w++) {
    double phi_topic_sum = 0.0;
    for (int k = 0; k < phi.nrow(); k++) {
      phi(k, w) = e_log_theta[k] + e_log_beta(k, w);  
      if (k > 0)
        phi_topic_sum = log_sum(phi_topic_sum, phi(k, w));
      else
        phi_topic_sum = phi(k, w);
    }
    for (int k = 0; k < phi.nrow(); k++) {
      phi(k, w) = phi(k, w) - phi_topic_sum;
    }
  }
  
  return phi;
}

// can be shared with dtm
// [[Rcpp::export]]
List gamma_update(NumericVector gamma_row, NumericMatrix phi, NumericVector n, double alpha) {
  NumericVector gamma_change(gamma_row.size());
  for (int k = 0; k < phi.nrow(); k++) {
    double acc = 0.0;
    for (int w = 0; w < phi.ncol(); w++) {
      acc += exp(phi(k, w)) * n[w];
    }
    //std::cout << acc << std::endl;
    gamma_change[k] = std::abs(gamma_row[k] - alpha - acc);
    gamma_row[k] = alpha + acc;
  }
  double avg_gamma_change = sum(gamma_change) / gamma_change.size();
  List result;
  result["avg_gamma_change"] = avg_gamma_change;
  result["gamma"] = gamma_row;
  return result;
}

NumericMatrix lambda_update(NumericMatrix lambda, NumericMatrix dtm, List phis, double eta) {
  for (int k = 0; k < lambda.nrow(); k++) {
    for (int w = 0; w < lambda.ncol(); w++) {
      lambda(k, w) = eta;
      for (int d = 0; d < phis.size(); d++) {
        NumericMatrix phi = as<NumericMatrix>(phis[d]);
        lambda(k, w) += dtm(d, w) * exp(phi(k, w));
      }
    }
  }
  
  return lambda;
}

double log_lik(NumericVector n, NumericMatrix phi, NumericVector gamma_row, NumericVector e_log_theta, 
               NumericMatrix e_log_beta, NumericMatrix lambda, double alpha, double eta, int D) {
  int K  = e_log_beta.nrow();
  int W  = e_log_beta.ncol();
  
  double term1 = 0.0;
  for (int w = 0; w < W; w++) {
    double acc = 0.0;
    for (int k = 0; k < K; k++) {
      acc += exp(phi(k, w)) * (e_log_theta[k] + e_log_beta(k, w) - phi(k, w));
//       if (phi(k, w) == 0.0)
//         std::cout << "k: " << k << " w: " << w << std::endl;
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
  for (int k = 0; k < K; k++) {
    double lamsum = 0.0;
    for (int w = 0; w < W; w++) {
      term3 += (eta - lambda(k, w)) * e_log_beta(k, w) + lgamma(lambda(k, w));
      lamsum += lambda(k, w);
    }
    term3 += - lgamma(lamsum);
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
    NumericMatrix mat(K, W);
    std::fill(mat.begin(), mat.end(), 1.0 / K);
    phis[d] = mat;
  }
  NumericMatrix gammas(D, K);
  std::fill(gammas.begin(), gammas.end(), 1.0);
  NumericMatrix lambda(K, W);
  NumericMatrix e_log_theta(D, K);
  NumericMatrix e_log_beta(K, W);
  NumericVector log_liks(em_max_iter);
  
  for (int k = 0; k < K; k++)
    lambda(k, _) = runif(W);
  
  double lold = 1.0;
  double lnew = 10.0;
  int iter = 0;
  e_log_beta = compute_e_log_beta(lambda, e_log_beta);
  for (int d = 0; d < D; d++) {
    e_log_theta(d, _) = compute_e_log_theta(gammas(d, _), e_log_theta(d, _));
  }
  
  while ((std::abs((lnew - lold) / lold)) > em_tol && iter < em_max_iter) {
    std::cout << "Fractional change in lhood: " << (lnew - lold) << std::endl;
    std::cout << "Iteration: " << iter << std::endl;
    lold = lnew;
    lnew = 0.0;
    
    for (int d = 0; d < D; d++) {
      double avg_gamma_change = 10.0;
      NumericMatrix phi = as<NumericMatrix>(phis[d]);
      NumericVector gamma_row(gammas.ncol(), 1.0);
      e_log_theta(d, _) = compute_e_log_theta(gamma_row, e_log_theta(d, _));
      int doc_iter = 0;
      while (avg_gamma_change > gam_tol && doc_iter < doc_max_iter) {
        phis[d] = phi_update(phi, e_log_theta(d, _), e_log_beta);
        List gamma_result = gamma_update(gamma_row, phis[d], dtm(d, _), alpha);
        NumericVector gamma_row = gamma_result["gamma"];
        gammas(d, _) = gamma_row; 
        e_log_theta(d, _) = compute_e_log_theta(gamma_row, e_log_theta(d, _));
        avg_gamma_change = gamma_result["avg_gamma_change"];
        doc_iter++;
        //print_vector(e_log_theta(d, _));
        //std::cout << avg_gamma_change << std::endl;
      } 
    }
    
    lambda = lambda_update(lambda, dtm, phis, eta);
    for (int d = 0; d < D; d++) {
      lnew += log_lik(dtm(d, _), phis[d], gammas(d, _), e_log_theta(d, _), e_log_beta,
                      lambda, alpha, eta, D);
    }
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
  return result;
}
