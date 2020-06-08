// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;


double phi(double x1, double x2)
{
  return (erf(x2/sqrt(2)) - erf(x1/sqrt(2)))/2;
}

// [[Rcpp::export]]

List nodeCount(const Eigen::Map<Eigen::MatrixXd> &X,
               const Eigen::Map<Eigen::MatrixXd> &Z,
               const Eigen::Map<Eigen::MatrixXd> &Vg,
               NumericVector parCount,
               double xmin,
               double xmax,
               double tmin,
               double tmax)
{
  NumericVector nvec(X.rows());
  NumericVector svec(X.rows());
  double inc = 1/(sqrt(X.rows()) * X.cols());
  double checkZero = inc/2;
  bool empty = true;
  bool sibempty = true;
  for (int i = 0; i < X.rows(); i++) {
    double curCount = parCount(i);
    if (curCount > checkZero) {
      svec(i) = curCount;
      for (int j = tmin-1; j < tmax; j++) {
        if ((X(i, j) >= xmin) & (X(i, j) < xmax)) {
          nvec(i) += inc;
          svec(i) -= inc;
          empty = false;
        }
      }
      if (svec(i) > checkZero) {
        sibempty = false;
      }
    }
  }
  // if (empty) {
  if (empty | sibempty) {
    return List::create(Named("Empty") = wrap(1));
  } else {
    Eigen::MatrixXd ZtX = Z.transpose() * as<Eigen::VectorXd>(nvec);
    Eigen::MatrixXd sZtX = Z.transpose() * as<Eigen::VectorXd>(svec);
    return List::create(Named("Empty") = wrap(0),
                        Named("Count") = wrap(nvec),
                        Named("ZtX") = wrap(ZtX),
                        Named("VgZtX") = wrap(Vg * ZtX),
                        Named("SibCount") = wrap(svec),
                        Named("SibZtX") = wrap(sZtX),
                        Named("SibVgZtX") = wrap(Vg * sZtX));
  }
}


// [[Rcpp::export]]

List nodeCountSE(const Eigen::Map<Eigen::MatrixXd> &X,
                 const Eigen::Map<Eigen::MatrixXd> &SE,
                 const Eigen::Map<Eigen::MatrixXd> &Z,
                 const Eigen::Map<Eigen::MatrixXd> &Vg,
                 NumericVector parCount,
                 double xmin,
                 double xmax,
                 double tmin,
                 double tmax)
{
  NumericVector nvec(X.rows());
  NumericVector svec(X.rows());
  double checkZero = 1.11e-16;
  bool empty = true;
  bool sibempty = true;


  for (int i = 0; i < X.rows(); i++) {
    double curCount = parCount(i);
    if (curCount > checkZero) {
      svec(i) = curCount;
      for (int j = tmin-1; j < tmax; j++) {
        if ((X(i, j) >= xmin) & (X(i, j) < xmax)) {
          double inc = phi((xmin - X(i,j)) / SE(i,j),
                           (xmax - X(i,j)) / SE(i, j))/(sqrt(X.rows()) * X.cols());
          nvec(i) += inc;
          svec(i) -= inc;
          empty = false;
        }
      }
      if (svec(i) > checkZero) {
        sibempty = false;
      }
    }
  }


  if (empty | sibempty) {
    return List::create(Named("Empty") = wrap(1));
  } else {
    Eigen::MatrixXd ZtX = Z.transpose() * as<Eigen::VectorXd>(nvec);
    Eigen::MatrixXd sZtX = Z.transpose() * as<Eigen::VectorXd>(svec);
    return List::create(Named("Empty") = wrap(0),
                        Named("Count") = wrap(nvec),
                        Named("ZtX") = wrap(ZtX),
                        Named("VgZtX") = wrap(Vg * ZtX),
                        Named("SibCount") = wrap(svec),
                        Named("SibZtX") = wrap(sZtX),
                        Named("SibVgZtX") = wrap(Vg * sZtX));
  }
}
