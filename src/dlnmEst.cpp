// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

double phi2(double x1, double x2)
{
  return (erf(x2/sqrt(2)) - erf(x1/sqrt(2)))/2;
}


// [[Rcpp::export]]
SEXP dlnmEst(arma::dmat dlnm,
                   arma::dvec xsplits,
                   int nlags,
                   int nsamp,
                   double center,
                   double se,
                   bool smooth,
                   bool dlm)
{
  int rows = dlnm.n_rows;
  int nsplits = xsplits.n_elem - 1;
  if (!smooth)
    center--;
  arma::dcube C(nlags, nsplits, nsamp); C.fill(0.0);
  arma::dmat centerMat(nlags, nsamp); centerMat.fill(0.0);

  // Fill in estimates
  for (int i = 0; i < rows; i++) {
    int iter = dlnm(i, 0) - 1;
    double xmin = dlnm(i, 2);
    double xmax = dlnm(i, 3);
    int tmin = dlnm(i, 4) - 1;
    int tmax = dlnm(i, 5);
    double est = dlnm(i, 6);
    for (int t = tmin; t < tmax; t++) {
      for (int x = 0; x < nsplits; x++) {
        if (dlm) {
          C(t, x, iter) += est * xsplits[x];
        } else {
          if (smooth) {
            C(t, x, iter) += phi2((xmin - xsplits[x]) / se,
              (xmax - xsplits[x]) / se) * est;
          } else if ((xsplits[x] >= xmin) & (xsplits[x] <= xmax)) {
            C(t, x, iter) += est;
          }
        }
      }
      // Center value if SE is defined
      if (smooth) {
        centerMat(t, iter) += phi2((xmin - center) / se,
                                   (xmax - center) / se) * est;
      }
    }
  }

  // Center
  if (!dlm) {
    double cen = 0;
    for (int i = 0; i < nsamp; i++) {
      for (int t = 0; t < nlags; t++) {
        if (!smooth) {
          cen = C(t, center, i);
        } else {
          cen = centerMat(t, i);
        }
        for (int x = 0; x < nsplits; x++) {
          C(t, x, i) -= cen;
        }
      }
    }
  }

  return wrap(C);
}
