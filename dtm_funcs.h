#include <util.h>

List compute_e_log_beta();
NumericVector compute_e_log_theta(NumericVector gammas, NumericVector e_log_theta);
void phi_update(NumericMatrix& phi, NumericVector e_log_theta, NumericMatrix e_log_beta, 
                NumericVector uniq_words);
double gamma_update(NumericVector& gamma_row, NumericMatrix phi, NumericVector n, double alpha,
                    NumericVector uniq_words);