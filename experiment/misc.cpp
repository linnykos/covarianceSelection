// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export()]]
double c_quantile(const arma::mat X, const double quantile = 1){
  arma::vec X_vec = arma::sort(arma::vectorise(X), "ascend");
  return(X_vec[(X.n_elem - 1) * quantile]);
}