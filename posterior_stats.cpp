#include "util.h"

// [[Rcpp::export]]
NumericMatrix term_score(List lda, int T) {
  // calculate term score (Blei and Lafferty 2009)
  NumericMatrix beta = lda["beta"];
  int W = beta.rows();
  int K = beta.cols();
  NumericMatrix term_score(W, K);
  
  for (int w = 0; w < W; w++) {
    double product = 0.0;
    for (int k = 0; k < K; k++) {
      product *= beta(w, k);
    }
    for (int k = 0; k < K; k++) {
      term_score(w, k) = beta(w, k) * log(beta(w, k) / pow(product, (1 / K)));
    }
  }
  
  return term_score;
}