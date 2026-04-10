#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dist_to_cor_inplace(NumericVector x) {
  const R_xlen_t n = x.size();
  double *ptr = REAL(x); // 直接获取底层指针

  for (R_xlen_t i = 0; i < n; ++i) {
    const double v = ptr[i];
    ptr[i] = 1.0 - 2.0 * v * v;
  }

  return x;
}

// [[Rcpp::export]]
NumericVector cor_to_dist_inplace(NumericVector x) {
  const R_xlen_t n = x.size();
  double *ptr = REAL(x);

  for (R_xlen_t i = 0; i < n; ++i) {
    double val = 0.5 * (1.0 - ptr[i]);
    ptr[i] = (val < 0.0) ? 0.0 : std::sqrt(val); 
  }

  return x;
}