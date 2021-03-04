#include <RcppEigen.h>
using namespace Rcpp;

// General function library:
// * sampling
// * densities
// * sets

int sampleInt(std::vector<double> probs, double totProb);
int sampleInt(Eigen::VectorXd probs);
double logPSplit(double alpha, double beta, int depth, bool terminal);
double logDirichletDensity(Eigen::VectorXd x, Eigen::VectorXd alpha);
void rHalfCauchyFC(double* x2, double a, double b, double* yInv = 0);
Eigen::VectorXd rDirichlet(Eigen::VectorXd alpha);
std::vector<std::vector<int> > 
  intersectAndDiff(std::vector<int> origVec, std::vector<int> newVec);
std::vector<int> cppIntersection(const IntegerVector& A, 
                                 const IntegerVector& B);