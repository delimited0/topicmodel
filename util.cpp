#include <util.h>

IntegerVector whichPositive(NumericVector x) {
  // return indexes of vector elements > 0
  IntegerVector v = seq(0, x.size()-1);
  return v[x > 0];
}

// [[Rcpp::export]]
IntegerVector whichEqual(IntegerVector x, int y) {
  // return indexes of vector elements == y
  IntegerVector v = seq(0, x.size()-1);
  return v[x==y];
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

  
double digamma(double x) {
  double p;
  x=x+6;
  p=1/(x*x);
  p=(((0.004166666666667*p-0.003968253986254)*p+
        0.008333333333333)*p-0.083333333333333)*p;
  p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
  return p;
}


double log_gamma(double x) {
  double z=1/(x*x);
  
  x=x+6;
  z=(((-0.000595238095238*z+0.000793650793651)
      *z-0.002777777777778)*z+0.083333333333333)/x;
  z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
    log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
  return z;
}