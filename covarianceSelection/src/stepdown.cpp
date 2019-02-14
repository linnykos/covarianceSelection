// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

// All the following functions assumes the matrix X is already centered
// Q: let's keep it as a matrix for now to see later if vectorizing is faster

// [[Rcpp::export()]]
arma::mat c_compute_sigma(const arma::mat X) {
  double n = X.n_rows;
  return(arma::trans(X) * X / n);
}

// [[Rcpp::export()]]
arma::mat c_compute_variance(const arma::mat X, const arma::mat cov_mat) {
  double n = X.n_rows;
  arma::mat X2 = arma::square(X);
  
  return(arma::trans(X2)*X2/n - arma::square(cov_mat));
}

// [[Rcpp::export()]]
arma::mat c_compute_bootSigma(const arma::mat X, const arma::vec noise_vec, 
                              const arma::mat cov_mat) {
  double n = X.n_rows;
  
  return(arma::trans(X)*arma::diagmat(noise_vec)*X/n - arma::sum(noise_vec)/n*cov_mat);
}