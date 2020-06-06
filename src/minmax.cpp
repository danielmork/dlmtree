#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector minmax(NumericVector x) {
  NumericVector range(2);
  range[0] = 2e16;
  range[1] = -2e16;
  for (auto &i : x) {
    if (i < range[0])
      range[0] = i;
    if (i > range[1])
      range[1] = i;
  }
  return range;
}
