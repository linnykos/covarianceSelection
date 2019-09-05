

// from https://github.com/fditraglia/selectr/blob/master/R/combn.cpp
// [[Rcpp::export]]
arma::umat c_combn(double n) {
  int k = 2;
  double n_subsets = Rf_choose(n, k);
  arma::umat out = arma::zeros<arma::umat>(k, n_subsets);
  arma::uvec a = arma::linspace<arma::uvec>(1, k, k);  
  out.col(0) = a;
  int m = 0;
  int h = k;
  arma::uvec j;
  
  for(long long i = 1; i < n_subsets; i++){
    if(m < (n - h)){  
      h = 1;
      m = a(k - 1);
      j = arma::linspace<arma::uvec>(1, 1, 1);
    }
    else{
      m = a(k - h - 1);
      ++h;
      j = arma::linspace<arma::uvec>(1, h, h);
    }
    a.elem(k - h - 1 + j) = m + j;
    out.col(i) = a;
  }
  return(out);
}

// [[Rcpp::export()]]
void c_armadillo(){
  arma::arma_version ver;
  std::cout << "ARMA version: "<< ver.as_string() << std::endl;
}


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

// [[Rcpp::export()]]
SEXP c_compute_bootSigma_tmp(SEXP X_, SEXP noise_, SEXP cov_){
  using namespace Rcpp;
  
  arma::mat X = as<arma::mat>(X_);
  arma::rowvec noise = as<arma::rowvec>(noise_);
  arma::mat cov = as<arma::mat>(cov_);
  
  
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