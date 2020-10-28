#include <Rcpp.h>
using namespace Rcpp;

//' Clamp vector to interval [-1, 1].
//'
//' @param x Numeric vector.
//' @return Numeric vector with values in [-1, 1].
// [[Rcpp::export]]
NumericVector clamp1(NumericVector x) {
  for(int i = 0; i < x.size(); ++i) {
    x[i] = std::max(-1.0, std::min(x[i], 1.0));
  }
  return x;
}
