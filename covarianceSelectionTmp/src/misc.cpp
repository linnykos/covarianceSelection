// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export()]]
double c_quantile(const arma::mat X, const double quantile = 1){
  arma::vec X_vec = arma::sort(arma::vectorise(X), "ascend");
  return(X_vec[(X.n_elem - 1) * quantile]);
}

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