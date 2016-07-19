#include "util.h"

using namespace Rcpp;

double phi_update(int d, NumericMatrix& phi, NumericVector uniq_words, NumericVector dtm_d,
                NumericMatrix& word_topic_count, NumericMatrix& doc_topic_count,
                NumericVector& topic_count, double alpha, double eta) {
  int W = uniq_words.size(); // number of unique words in document d
  int K = phi.ncol();
  NumericMatrix phi_old(phi.rows(), phi.cols());
  NumericMatrix phi_change(phi.rows(), phi.cols());
  
  for (int w = 0; w < W; w++) {
    NumericVector phi_new(K);
    for (int k = 0; k < K; k++) {
      phi_old(w, k) = phi(w, k);
      double remove_count = dtm_d[uniq_words[w]] * phi(w, k);
      word_topic_count(uniq_words[w], k) -= remove_count;
      doc_topic_count(d, k) -= remove_count;
      topic_count[k] -= remove_count;
      phi_new[k] = (eta + word_topic_count(uniq_words[w], k)) * 
             (alpha + doc_topic_count(d, k)) / 
             (W * eta + topic_count[k]);
    }
    double denom = sum(phi_new);
    for (int k = 0; k < K; k++) {
      phi(w, k) = phi_new(k) / denom;
      double add_count = dtm_d[uniq_words[w]] * phi(w, k);
      word_topic_count(uniq_words[w], k) += add_count;
      doc_topic_count(d, k) += add_count;
      topic_count[k] += add_count;
    }
  }
  
  for (int w = 0; w < W; w++)
    for (int k = 0; k < K; k++) 
      phi_change(w, k) = std::abs(phi(w, k) - phi_old(w, k));
  
  double avg_phi_change = sum(phi_change) / (phi_change.rows() * phi_change.cols());
  return avg_phi_change;
}

double log_lik(List phis, List doc_words, NumericMatrix dtm, NumericMatrix word_topic_count, 
               NumericMatrix doc_topic_count, NumericVector topic_count, 
               double alpha, double eta) {
  int D = doc_topic_count.rows();
  int W = word_topic_count.rows();
  int K = topic_count.size();
  
  double accum = 0;
  for (int d = 0; d < D; d++) {
    for (int k = 0; k < K; k++) {
      accum += log_gamma(alpha + doc_topic_count(d, k));
    }
    accum -= log_gamma(K * alpha + sum(doc_topic_count(d, _)));
    
    NumericMatrix phi = as<NumericMatrix>(phis[d]);
    NumericVector uniq_words = as<NumericVector>(doc_words[d]);
    
    for (int j = 0; j < phi.rows(); j++) {
      for (int k = 0; k < K; k++) {
        accum -= dtm(d, uniq_words[j]) * (phi(j, k) * log(phi(j, k)));
      }
    }
  }
  
  for (int k = 0; k < K; k++) {
    for (int w = 0; w < W; w++) {
      accum += log_gamma(eta + word_topic_count(w, k));
    }
    accum -= log_gamma(W * eta + topic_count(k));
  }
  
  return accum;
}

// [[Rcpp::export]]
List lda_cvb0(NumericMatrix dtm, int K, double alpha, double eta, double phi_tol, 
              double em_tol, double em_max_iter, double doc_max_iter) {
  int D = dtm.nrow();
  int W = dtm.ncol();  
  
  NumericVector log_liks(em_max_iter);
  NumericMatrix word_topic_count(W, K);
  NumericMatrix doc_topic_count(D, K);
  NumericVector topic_count(K);
  List phis(D);
  List doc_words(D);
  for (int d = 0; d < D; d++) {
    IntegerVector uniq_words = whichPositive(dtm(d, _));
    int nonzero = uniq_words.size();
    NumericMatrix mat(nonzero, K);
    for (int w = 0; w < nonzero; w++) {
      mat(w, _) = runif(K);
      double denom = sum(mat(w, _));
      for (int k = 0; k < K; k++)
        mat(w, k) = mat(w, k) / denom;
    }
    phis[d] = mat;
    doc_words[d] = uniq_words;
    
    for (int w = 0; w < uniq_words.size(); w++) {
      for (int k = 0; k < K; k++) {
        word_topic_count(uniq_words[w], k) += dtm(d, uniq_words[w]) * mat(w, k);
        doc_topic_count(d, k) += dtm(d, uniq_words[w]) * mat(w, k);
        topic_count[k] += dtm(d, uniq_words[w]) * mat(w, k);
      }
    }
  }
  
  double lold = 0;
  double lnew = 0;
  int iter = 0;
  
  while (((std::abs((lnew - lold) / lold)) > em_tol && iter < em_max_iter) || iter == 0) {
    Rcout << "Iteration: " << iter << std::endl;
    lold = lnew;
    lnew = 0.0;
    
    for (int d = 0; d < D; d++) {
      NumericMatrix phi = as<NumericMatrix>(phis[d]);
      double avg_phi_change = phi_tol * 100.0;
      if (d % 100 == 0)
        Rcout << "Fitting Document: " << d << std::endl;
      // phi_update(phi, doc_words[d], dtm(d, _), word_topic_count, doc_topic_count_d,
      //            topic_count, alpha, eta);
      int doc_iter = 0;
      while (avg_phi_change > phi_tol && doc_iter < doc_max_iter) {
        avg_phi_change = phi_update(d, phi, doc_words[d], dtm(d, _), word_topic_count, 
                                    doc_topic_count, topic_count, alpha, eta);
        doc_iter++;
      }
    }
    
    lnew = log_lik(phis, doc_words, dtm, word_topic_count, 
                   doc_topic_count, topic_count, alpha, eta);
    log_liks[iter] = lnew;
    iter++;
    Rcout << "Log Likelihood: " << lnew << std::endl;
  }
  
  NumericMatrix beta(W, K);
  NumericMatrix theta(D, K);
  for (int k = 0; k < K; k++) {
    for (int w = 0; w < W; w++) {
      beta(w, k) = (word_topic_count(w, k) + eta) / (topic_count(k) + W * eta);  
    }
    for (int d = 0; d < D; d++) {
      theta(d, k) = (doc_topic_count(d, k) + alpha) / (sum(doc_topic_count(d, _)) + K * alpha);
    }
  }
  
  List result;
  result["phis"] = phis;
  result["doc_words"] = doc_words;
  result["beta"] = beta;
  result["theta"] = theta;
  result["log_liks"] = log_liks;
  result["word_topic_count"] = word_topic_count;
  result["doc_topic_count"] = doc_topic_count;
  result["topic_count"] = topic_count;
  CharacterVector class_names(2);
  class_names[0] = "lda_cvb0"; 
  class_names[1] = "lda";
  result.attr("class") = class_names;
  return result;
}
