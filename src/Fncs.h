#include <RcppEigen.h>
using namespace Rcpp;

// General function library:
// * sampling
// * densities
// * sets

int sampleInt(const std::vector<double> &probs, double totProb);
int sampleInt(const Eigen::VectorXd &probs);
double logPSplit(double alpha, double beta, int depth, bool terminal);
double logDirichletDensity(const Eigen::VectorXd &x, 
                           const Eigen::VectorXd &alpha);
void rHalfCauchyFC(double* x2, double a, double b, double* yInv = 0);
Eigen::VectorXd rDirichlet(const Eigen::VectorXd &alpha);
// std::vector<std::vector<int> > 
std::pair<std::vector<int>, std::vector<int> >
  intersectAndDiff(const std::vector<int> &origVec, 
                   const std::vector<int> &newVec);
std::vector<int> cppIntersection(const IntegerVector& A, 
                                 const IntegerVector& B);