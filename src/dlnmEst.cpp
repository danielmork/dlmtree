// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588


double phi2(double x1, double x2)
{
  return (erf(x2 * MATH_SQRT1_2) - erf(x1 * MATH_SQRT1_2)) * 0.5;
}

// [[Rcpp::export]]
SEXP dlnmEst(arma::dmat dlnm,
             arma::dvec predAt,
             int nlags,
             int nsamp,
             double center,
             double se)
{
  int rows = dlnm.n_rows;
  int nsplits;
  bool smooth = 0;
  nsplits = predAt.n_elem;
  arma::dcube C(nlags, nsplits, nsamp); C.fill(0.0);
  arma::dmat centerMat(nlags, nsamp);

  if (se > 0) {
    smooth = 1;
    centerMat.fill(0.0);
  } else {
    center--;
  }

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
        if (smooth) {
          C(t, x, iter) += phi2((xmin - predAt[x]) / se,
                                (xmax - predAt[x]) / se) * est;
        } else if ((xmin <= predAt[x]) && (xmax > predAt[x])) {
          C(t, x, iter) += est;
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
  double cen = 0;
  for (int i = 0; i < nsamp; i++) {
    for (int t = 0; t < nlags; t++) {
      if (smooth) {
        cen = centerMat(t, i);
      } else {
        cen = C(t, center, i);
      }
      for (int x = 0; x < nsplits; x++) {
        C(t, x, i) -= cen;
      }
    }
  }

  return wrap(C);
}

// [[Rcpp::export]]
SEXP splitPIP(arma::dmat dlnm,
              int nlags)
{
  int rows = dlnm.n_rows;
  int tree = 0;
  int iter = 0;
  arma::vec splitCount(nlags);
  arma::vec splitIter(nlags);

  for (int i = 0; i < rows; ++i) {
    if (dlnm(i, 0) > iter) {
      splitCount += splitIter;
      splitIter.zeros();
      iter = dlnm(i, 0);
    }
    for (int t = dlnm(i, 4) - 1; t < dlnm(i, 5); ++t) {
      if (splitIter(t) == 0)
        splitIter(t) = 1.0;
    }
  }
  
  return wrap(splitCount);
}

// [[Rcpp::export]]
SEXP dlnmPLEst(arma::dmat dlnm,
               arma::dvec predAt,
               int nlags,
               int nsamp,
               double center)
{
  int rows = dlnm.n_rows;
  int nsplits;
  bool smooth = 0;
  nsplits = predAt.n_elem;
  arma::dcube C(nlags, nsplits, nsamp); C.fill(0.0);
  arma::dmat centerMat(nlags, nsamp);
  center--;
  double prevEst = 0.0;
  int iter = dlnm(0, 0);
  int tree = dlnm(0, 0);

  // Fill in estimates
  for (int i = 0; i < rows; i++) {
    if ((int(dlnm(i, 0) - 1) != iter) || (int(dlnm(i, 1)) != tree)) {
      prevEst = 0.0;
    }
    iter = dlnm(i, 0) - 1;
    tree = dlnm(i, 1);
    double xmin = dlnm(i, 2);
    double xmax = dlnm(i, 3);
    double den = xmax - xmin;
    int tmin = dlnm(i, 4) - 1;
    int tmax = dlnm(i, 5);
    double est = dlnm(i, 6);

    for (int t = tmin; t < tmax; t++) {
      for (int x = 0; x < nsplits; x++) {
        if ((xmin <= predAt[x]) && (xmax > predAt[x])) {
          C(t, x, iter) += prevEst + (est - prevEst) * (predAt[x] - xmin) / den;
        }
      }
    }
    prevEst = est;
  } // end loop over dlnm tree output

  // Center
  double cen = 0;
  for (int i = 0; i < nsamp; i++) {
    for (int t = 0; t < nlags; t++) {
      cen = C(t, center, i);
      for (int x = 0; x < nsplits; x++) {
        C(t, x, i) -= cen;
      }
    }
  }

  return wrap(C);
}


// [[Rcpp::export]]
SEXP dlmEst(arma::dmat dlm,
            int nlags,
            int nsamp)
{
  int rows = dlm.n_rows;
  arma::dmat C(nlags, nsamp); C.fill(0.0);

  // Fill in estimates
  for (int i = 0; i < rows; i++) {
    int iter = dlm(i, 0) - 1;
    int tmin = dlm(i, 2) - 1;
    int tmax = dlm(i, 3);
    double est = dlm(i, 4);
    for (int t = tmin; t < tmax; t++) {
      C(t, iter) += est;
    }
  }

  return wrap(C);
}

// [[Rcpp::export]]
SEXP mixEst(arma::dmat dlm,
            int nlags,
            int nsamp)
{
  int rows = dlm.n_rows;
  arma::dcube C(nlags, nlags, nsamp); C.fill(0.0);

  int i, t1, t2, iter, tmin1, tmax1, tmin2, tmax2;
  // Fill in estimates
  for (i = 0; i < rows; i++) {
    iter = dlm(i, 0) - 1;
    tmin1 = dlm(i, 3) - 1;
    tmax1 = dlm(i, 4);
    tmin2 = dlm(i, 6) - 1;
    tmax2 = dlm(i, 7);
    double est = dlm(i, 8);
    for (t1 = tmin1; t1 < tmax1; t1++) {
      for (t2 = tmin2; t2 < tmax2; t2++) {
        C(t1, t2, iter) += est;
      }
    }
  }

  return wrap(C);
}
