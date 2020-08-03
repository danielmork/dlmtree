#include <RcppEigen.h>
#include "Fncs.h"
using namespace Rcpp;

int sampleInt(std::vector<double> p, double totP = 1)
{
  double u = R::runif(0, totP);
  double sum = p[0];
  int i = 0;
  while (sum < u) {
    i++;
    sum += p[i];
  }
  return(i);
}

int sampleInt(Eigen::VectorXd p)
{
  double totP = p.sum();
  double u = R::runif(0, totP);
  double sum = p[0];
  int i = 0;
  while (sum < u) {
    i++;
    sum += p[i];
  }
  return(i);
}

double logPSplit(double alpha, double beta, int depth, bool terminal)
{
  double p = alpha * pow(1.0 + (double)depth, -beta);
  if (terminal) {
    return(log1p(-p));
  } else {
    return(log(p));
  }
}

double logDirichletDensity(Eigen::VectorXd x, Eigen::VectorXd alpha)
{
  if (x.size() != alpha.size())
    stop("logDirichletDensity incorrect size");
  double sumAlpha = alpha.sum();
  double sumAlphaX = ((alpha.array() - 1) * x.array().log()).matrix().sum();
  double out = sumAlphaX + lgamma(sumAlpha);
  for (int i = 0; i < alpha.size(); i++)
    out -= lgamma(alpha(i));
  return(out);
}

Eigen::VectorXd rDirichlet(Eigen::VectorXd alpha)
{
  Eigen::VectorXd out(alpha.size());
  double norm = 0;
  for (int i = 0; i < alpha.size(); i++) {
    out(i) = R::rgamma(alpha(i), 1);
    norm += out(i);
  }
  out /= norm;
  return(out);
}
