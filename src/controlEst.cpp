// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
//' @param R Vector of responses
//' @param Z Model matrix of covariates
//' @param V Posterior covariance for gamma
//' @param Vchol chol(V)
//' @param sigma2 Current estimate of sigma2
//' @param sumTermXP Number of terminal nodes in model times number of DLM times
//' @param sumFexp Sum over trees of theta * lambda inverse * theta transpose / tau
//'
//' @return List of
//' * new sigma2 estimate
//' * new gamma estimate
//' * new xi estimate (sigma2 hyperparameter)
// [[Rcpp::export]]

List controlEst(const Eigen::Map<Eigen::VectorXd> R,
                   const Eigen::Map<Eigen::MatrixXd> Z,
                   const Eigen::Map<Eigen::MatrixXd> V,
                   const Eigen::Map<Eigen::MatrixXd> Vchol,
                   double sigma2,
                   double sumTermXP,
                   double sumFexp)
{
  int n = R.size();
  int p = Z.cols();
  const Eigen::VectorXd ZR = Z.transpose() * R;
  const Eigen::VectorXd gammaHat = V * ZR;
  double xiSigma2 = 1 / as<double>(rgamma(1, 1, sigma2/(sigma2 + 1)));
  double sigma2New = 1 / as<double>(rgamma(1, (n + sumTermXP + 1)/2, 1/((R.dot(R) - ZR.dot(gammaHat) + sumFexp) / 2 + 1 / xiSigma2)));
  const Eigen::VectorXd gamma = Vchol * as<Eigen::VectorXd>(rnorm(p, 0, sqrt(sigma2New))) + gammaHat;

  return List::create(Named("sigma2") = wrap(sigma2New),
                      Named("xi.sigma2") = wrap(xiSigma2),
                      Named("gamma") = wrap(gamma));
}
