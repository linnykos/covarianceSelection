// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "misc.h"

// All the following functions assumes the matrix X is already centered
// Q: let's keep it as a matrix for now to see later if vectorizing is faster

// [[Rcpp::export()]]
arma::mat c_compute_sigma(const arma::mat& X) {
  double n = X.n_rows;
  return(X.t() * X / n);
}
 
// [[Rcpp::export()]]
arma::mat c_compute_variance(const arma::mat& X, const arma::mat& cov_mat) {
  double n = X.n_rows;
  arma::mat X2 = arma::square(X);
  
  return(X2.t() * X2/n - arma::square(cov_mat));
}

// [[Rcpp::export()]]
arma::mat c_compute_bootSigma(const arma::mat& X, const arma::vec& noise_vec, 
                              const arma::mat& cov_mat) {
  double n = X.n_rows;
  
  return(X.t() * arma::diagmat(noise_vec)*X/n - arma::sum(noise_vec)/n*cov_mat);
}

// // [[Rcpp::export()]]
// arma::mat c_compute_bootSigma_tmp(const arma::mat& X, const arma::vec& noise_vec, 
//                               const arma::mat& cov_mat) {
//   double n = X.n_rows;
//   arma::mat X2(X.n_rows, X.n_cols);
//   for(int i = 0; i < X.n_cols; i++){
//     X2.col(i) = noise_vec[i]*X/n;
//   }
//   
//   return(X.t() * X2 - arma::sum(noise_vec)/n*cov_mat);
// }

// [[Rcpp::export()]]
double c_compute_covStat(const arma::mat& num_x, const arma::mat& num_y,
                            const arma::mat& denom_x, const arma::mat& denom_y,
                            const double& quantile = 1, const bool& squared = true){
  if(squared){
    arma::mat res = arma::square(num_x - num_y)/(denom_x + denom_y);
    return(c_quantile(res));
  } else {
    arma::mat res = arma::abs(num_x - num_y);
    return(c_quantile(res));
  }
}

///////

// [[Rcpp::export()]]
Rcpp::List c_compute_all_denom(const Rcpp::List& dat_list, const Rcpp::List& cov_list){
  int len = dat_list.size();
  Rcpp::List res(len);
  for(int i = 0; i < len; i++){
    res[i] = c_compute_variance(dat_list[i], cov_list[i]);
  }
  
  return(res);
}

// [[Rcpp::export()]]
arma::vec c_compute_all_test_stat(const Rcpp::List& num_list, const Rcpp::List& denom_list,
                                   const arma::umat& combn_mat, const bool& squared = true){
  
  int len = combn_mat.n_cols;
  arma::vec res = arma::vec(len);
  for(int i = 0; i < len; i++){
    int j = combn_mat(0,i) - 1;
    int k = combn_mat(1,i) - 1;
    res[i] = c_compute_covStat(Rcpp::as<arma::mat>(num_list[j]),
                               Rcpp::as<arma::mat>(num_list[k]),
                               Rcpp::as<arma::mat>(denom_list[j]),
                               Rcpp::as<arma::mat>(denom_list[k]),
                                                   squared);
  }
  
  return(res);
}

// [[Rcpp::export()]]
Rcpp::List c_compute_all_numerator_bootstrap(const Rcpp::List& dat_list, const Rcpp::List& noise_list,
                                             const Rcpp::List& cov_list, const arma::uvec& remaining_idx){
  int k = remaining_idx.size();
  Rcpp::List res(k);
  
  for(int i = 0; i < k; i++){
    int j = remaining_idx[i] - 1;
    res[i] = c_compute_bootSigma(Rcpp::as<arma::mat>(dat_list[j]), 
                                 Rcpp::as<arma::vec>(noise_list[j]),
                                 Rcpp::as<arma::mat>(cov_list[j]));
  }
  
  return(res);
}

