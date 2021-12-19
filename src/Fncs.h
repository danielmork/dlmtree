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
std::vector<double> rtruncnorm(int n, double mu, double sigma, double lower, double upper);
double rtruncnorm(double mu, double sigma, double lower, double upper);
double dtruncnorm(double x, double mu, double sigma,  double lower, double upper);
Eigen::VectorXd rDirichlet(Eigen::VectorXd alpha);
std::vector<std::vector<int> > 
  intersectAndDiff(std::vector<int> origVec, std::vector<int> newVec);
std::vector<int> cppIntersection(const IntegerVector& A, 
                                 const IntegerVector& B);
Eigen::VectorXd selectInd(Eigen::VectorXd original, std::vector<int> indices);
Eigen::MatrixXd selectIndM(Eigen::MatrixXd original, std::vector<int> indices);