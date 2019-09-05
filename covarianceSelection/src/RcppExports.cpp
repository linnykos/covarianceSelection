// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// c_quantile
double c_quantile(const arma::mat X, const double quantile);
RcppExport SEXP _covarianceSelection_c_quantile(SEXP XSEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(c_quantile(X, quantile));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_covarianceSelection_c_quantile", (DL_FUNC) &_covarianceSelection_c_quantile, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_covarianceSelection(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}