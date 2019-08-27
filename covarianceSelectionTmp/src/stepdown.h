#ifndef _STEPDOWN_H
#define _STEPDOWN_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "misc.h"

arma::mat c_compute_sigma(const arma::mat X);
arma::mat c_compute_variance(const arma::mat X, const arma::mat cov_mat);
arma::mat c_compute_bootSigma(const arma::mat X, const arma::vec noise_vec, 
                              const arma::mat cov_mat);
double c_compute_covStat(const arma::mat num_x, const arma::mat num_y,
                         const arma::mat denom_x, const arma::mat denom_y,
                         const double quantile = 1, const bool squared = TRUE);

#endif