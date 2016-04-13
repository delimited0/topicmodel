//
// dynamic topic model
//

#include "util.h"
#include "dtm_funcs.h"

List dyn_tm(NumericMatrix dtm, IntegerVector period, int K, double alpha, double sigma, 
            double nu, double gam_tol, double em_tol,
            double em_max_iter, double doc_max_iter) {
  int D = dtm.nrow();
  int W = dtm.ncol();
  IntegerVector period_uniq = unique(period);
  int T = period_uniq.size();
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
  NumericMatrix e_log_theta(D, K);
  List betas_hat(T);
  List e_log_beta(T);
  List m(T);
  List v(T);
  List doc_period_id(T);
  for (int t = 0; t < T; t++) {
    NumericMatrix mat1(W, K);
    for (int k = 0; k < K; k++)
      mat1(_, k) = runif(W);
    betas_hat[t] = mat1;
    NumericMatrix mat2(W, K);
    m[t] = mat2;
    NumericMatrix mat3(W, K);
    v[t] = mat3;
    NumericMatrix mat4(W, K);
    e_log_beta[t] = mat4;
    IntegerVector period_t = whichEqual(period, t);
    doc_period_id[t] = period_t;
  }
  NumericMatrix zeta(T, K);
  
  for (int p = 0; p < period_uniq.length(); p++) {
    IntegerVector period_p = whichEqual(period, p);
    NumericMatrix mat4()
  }
  
  double lold = 0;
  double lnew = 0;
  int iter = 0;
  
  while () {
    for (int t = 0; t < T; t++) {
      Numeric Vector doc_ids = doc_period_id[t];
      for (int d = 0; d < doc_ids.size(); d++) {
        
      }
    }
  }
}