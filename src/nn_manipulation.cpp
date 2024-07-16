#include <Rcpp.h>
#include <RcppCommon.h>

// [[Rcpp::export]]
Rcpp::List rcpp_get_sparse_dist(
    Rcpp::IntegerMatrix & knn_index,
    Rcpp::NumericMatrix & knn_dist,
    int n_obs,
    int n_neighbors
) {
  Rcpp::IntegerVector rows(n_obs * n_neighbors);
  Rcpp::IntegerVector cols(n_obs * n_neighbors);
  Rcpp::NumericVector vals(n_obs * n_neighbors);

  int nrow = knn_index.nrow();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < n_neighbors; ++j) {
      if (knn_index(i, j) == 0) continue;
      if (knn_index(i, j) == i + 1) {
        vals[i * n_neighbors + j] = 0;
      } else {
        vals[i * n_neighbors + j] = knn_dist(i, j);
      }
      cols[i * n_neighbors + j] = i;
      rows[i * n_neighbors + j] = knn_index(i, j) - 1;
    }
  }
  Rcpp::List res = Rcpp::List::create(
    Rcpp::_["i"] = rows,
    Rcpp::_["j"] = cols,
    Rcpp::_["x"] = vals
  );
  return res;
}
