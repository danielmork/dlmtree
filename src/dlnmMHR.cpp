// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
//' @param
//'
//' @return
// [[Rcpp::export]]

List dlnmMHR1(const Eigen::Map<Eigen::VectorXd> &Y,
              const Eigen::Map<Eigen::MatrixXd> &Z,
              double &XX,
              const Eigen::Map<Eigen::VectorXd> &ZtX,
              const Eigen::Map<Eigen::VectorXd> &VgZtX,
              const Eigen::Map<Eigen::MatrixXd> &Vg,
              const Eigen::Map<Eigen::VectorXd> &Xd,
              const Eigen::Map<Eigen::VectorXd> &ZY,
              double LInv,
              double sigma)
{
  // const Eigen::VectorXd ZtX = Z.transpose() * Xd;
  // const Eigen::VectorXd VgZtX = Vg * ZtX;
  double VTheta = 1/(XX - ZtX.dot(VgZtX) + LInv);
  Eigen::VectorXd ZtY(Y.cols());
  if (ZY.size() == 1) {
    ZtY = Z.transpose() * Y;
  } else {
    ZtY = ZY;
  }
  const double XtVzInvY = (Xd.dot(Y)) - (VgZtX.dot(ZtY));
  const double ThetaHat = VTheta * XtVzInvY;
  const double VThetaChol = sqrt(VTheta);
  const double ThetaDraw = VThetaChol * as<double>(rnorm(1, 0, sigma)) + ThetaHat;
  const Eigen::VectorXd Yhat = Xd * ThetaDraw;


  return List::create(Named("ThetaDraw") = wrap(ThetaDraw),
                      Named("VThetaLogDet") = wrap(log(VThetaChol)),
                      Named("Beta") = wrap((Y.dot(Y)) - (ZtY.dot(Vg * ZtY)) - (ThetaHat * XtVzInvY)),
                      Named("TLiT") = wrap(ThetaDraw * ThetaDraw),
                      Named("ZtY") = wrap(ZtY),
                      Named("Yhat") = wrap(Yhat));
}

//' @param
//'
//' @return
// [[Rcpp::export]]

List dlnmMHR(const Eigen::Map<Eigen::VectorXd> &Y,
             const Eigen::Map<Eigen::MatrixXd> &Z,
             const Eigen::Map<Eigen::MatrixXd> &ZtX,
             const Eigen::Map<Eigen::MatrixXd> &VgZtX,
             const Eigen::Map<Eigen::MatrixXd> &Vg,
             const Eigen::Map<Eigen::MatrixXd> &Xd,
             const Eigen::Map<Eigen::VectorXd> &ZY,
             double LInv,
             double sigma)
{
  int pX = Xd.cols();
  // const Eigen::MatrixXd ZtX = Z.transpose() * Xd;
  // const Eigen::MatrixXd VgZtX = Vg * ZtX;
  Eigen::MatrixXd tempV = Xd.transpose() * Xd - ZtX.transpose() * VgZtX;
  tempV.diagonal().array() += LInv;
  const Eigen::MatrixXd VTheta = tempV.inverse();
  Eigen::VectorXd ZtY(Y.cols());
  if (ZY.size() == 1) {
    ZtY = Z.transpose() * Y;
  } else {
    ZtY = ZY;
  }
  const Eigen::VectorXd XtVzInvY = ((Xd.transpose() * Y) - (VgZtX.transpose() * ZtY));
  const Eigen::VectorXd ThetaHat = VTheta * XtVzInvY;
  const Eigen::MatrixXd VThetaChol = VTheta.llt().matrixL();
  const Eigen::VectorXd ThetaDraw = VThetaChol * as<Eigen::VectorXd>(rnorm(pX, 0, sigma)) + ThetaHat;
  const Eigen::VectorXd Yhat = Xd * ThetaDraw;


  return List::create(Named("ThetaDraw") = wrap(ThetaDraw),
                      Named("VThetaLogDet") = wrap(VThetaChol.diagonal().array().log().sum()),
                      Named("Beta") = wrap(Y.dot(Y) - ZtY.dot(Vg * ZtY) - ThetaHat.dot(XtVzInvY)),
                      Named("TLiT") = wrap(ThetaDraw.dot(ThetaDraw)),
                      Named("ZtY") = wrap(ZtY),
                      Named("Yhat") = wrap(Yhat));
}


//' @param
//'
//' @return
// [[Rcpp::export]]

// List dlnmMHRCalc (int nOldNodes,
//                   int nNewNodes,
//                   double treeVar,
//                   double sigma,
//                   )
