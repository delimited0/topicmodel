#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix compute_e_log_beta(NumericMatrix lambda, NumericMatrix e_log_beta);
NumericVector compute_e_log_theta(NumericVector gammas, NumericVector e_log_theta);
NumericMatrix phi_update(NumericMatrix phi, NumericVector e_log_theta, NumericMatrix e_log_beta);
List gamma_update(NumericVector gamma_row, NumericMatrix phi, NumericVector n, double alpha);