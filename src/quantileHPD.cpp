// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP quantileHPD(arma::dvec samp,
         double prob)
{
  samp = sort(samp);
  int n = samp.size();
  int upper = ceil(n * prob);
  int stop = n - upper;
  arma::dvec lims(2);
  double diff = 1e16;
  for (int i = 0; i < stop; i++) {
    if ((samp[upper + i] - samp[i]) < diff) {
      lims[0] = samp[i];
      lims[1] = samp[upper + i];
      diff = (samp[upper + i] - samp[i]);
    }
  }
  return wrap(lims);
}

