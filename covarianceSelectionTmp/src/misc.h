#ifndef _MISC_H
#define _MISC_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

double c_quantile(const arma::mat X, const double quantile = 1);

#endif