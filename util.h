#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

IntegerVector whichPositive(NumericVector x);
IntegerVector whichEqual(IntegerVector x, int y);
double log_sum(double log_a, double log_b);
double digamma(double x);
double log_gamma(double x);