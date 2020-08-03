#include <RcppEigen.h>

int sampleInt(std::vector<double>, double);

int sampleInt(Eigen::VectorXd);

double logPSplit(double, double, int, bool);

double logDirichletDensity(Eigen::VectorXd, Eigen::VectorXd);

Eigen::VectorXd rDirichlet(Eigen::VectorXd);
