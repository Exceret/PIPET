// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat corCosine_cpp(const arma::mat &x, Nullable<arma::mat> y = R_NilValue) {

  // Helper lambda: compute L2 norms with zero-protection
  auto compute_norms = [](const arma::mat &m) -> arma::vec {
    arma::vec norms =
        arma::sqrt(arma::sum(arma::square(m), 0).t()); // Key fix: .t()
    for (auto &v : norms) {
      if (v == 0.0)
        v = 1.0; // Avoid division by zero
    }
    return norms;
  };

  const auto x_norm = compute_norms(x);

  if (y.isNull()) {
    // Self-correlation
    const arma::mat cross = x.t() * x;
    const arma::mat denom = x_norm * x_norm.t();
    return cross / denom;
  }

  // Properly unwrap Nullable<arma::mat>
  const arma::mat y_mat = Rcpp::as<arma::mat>(y);

  if (x.n_rows != y_mat.n_rows) {
    stop("x and y must have the same number of rows (got %d and %d)",
         static_cast<int>(x.n_rows), static_cast<int>(y_mat.n_rows));
  }

  const auto y_norm = compute_norms(y_mat);
  const arma::mat cross = x.t() * y_mat;
  const arma::mat denom = x_norm * y_norm.t();

  return cross / denom;
}