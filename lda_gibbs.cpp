// LDA, collapsed Gibbs sampler

#include "util.h"

// [[Rcpp::export]]
List lda_gibbs(NumericMatrix dtm, int K, double alpha, double eta, int burnin, int iter) {
  int D = dtm.nrow();
  int W = dtm.ncol();  
  // maybe replace these matrices with sparse ones
  NumericMatrix word_topic_count(W, K);
  NumericMatrix doc_topic_count(D, K);
  NumericVector topic_count(K);
  
  int tot_words = 0;
  for (int d = 0; d < D; d++) {
    tot_words += sum(dtm(d, _));
  }
  
  NumericVector word_assign = floor(runif(tot_words) * K);
  NumericVector word_id(tot_words);
  NumericVector doc_id(tot_words);
  int i = 0;
  for (int d = 0; d < D; d++) {
    Rcout << "Initializing doc " << d << std::endl;
    for (int w = 0; w < W; w++) {
      int word_count = dtm(d, w);
      while (word_count > 0) {
        doc_id[i] = d;
        word_id[i] = w;
        word_topic_count(w, word_assign[i])++;
        doc_topic_count(d, word_assign[i])++;
        word_count--;
        i++;
      }
    }
  }
  
  for (int k = 0; k < K; k++) {
    topic_count[k] = sum(word_topic_count(_, k));
  }
  Rcout << "Initialized" << std::endl;
  
  NumericVector us = runif(tot_words);
  for (int j = 0; j < iter; j++) {
    //if (j+1 % 100 == 0)
      Rcout << "Iteration: " << j << std::endl;
    for (int i = 0; i < tot_words; i++) {
      NumericVector Z(K);
      for (int k = 0; k < K; k++) {
        int subtract = 0;
        if (word_assign[i] == k)
          subtract = 1;
        Z[k] = (eta + word_topic_count(word_id[i], k) - subtract) * 
          (alpha + doc_topic_count(doc_id[i], k) - subtract) /
          (W * eta + topic_count[k] - subtract);
        if (k != 0)
          Z[k] = Z[k] + Z[k-1];
      }
      for (int k = 0; k < K; k++) {
        if (us[i] < Z[k] / Z[K-1]) {
          word_topic_count(word_id[i], word_assign[i])--;
          doc_topic_count(doc_id[i], word_assign[i])--;
          word_assign[i] = k;
          word_topic_count(word_id[i], word_assign[i])++;
          doc_topic_count(doc_id[i], word_assign[i])++;
          break;
        }
      }
    }
  }
  
  NumericMatrix beta(W, K);
  NumericMatrix theta(D, K);
  for (int k = 0; k < K; k++) {
    for (int w = 0; w < W; w++) {
      beta(w, k) = (word_topic_count(w, k) + eta) / (sum(word_topic_count(_, k)) + W * eta);  
    } 
    for (int d = 0; d < D; d++) {
      theta(d, k) = (doc_topic_count(d, k) + alpha) / (sum(doc_topic_count(d, _)) + K * alpha);
    }
  }
  
  List result;
  result["beta"] = beta;
  result["theta"] = theta;
  CharacterVector class_names(2); 
  class_names[0] = "lda_gibbs"; 
  class_names[1] = "lda";
  result.attr("class") = class_names;
  
  return result;
}