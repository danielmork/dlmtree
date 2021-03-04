#include <RcppEigen.h>
#include "Fncs.h"
using namespace Rcpp;

/**
 * sampleInt
 * @brief Method to sample integer using probabilities p
 * @param probs vector of probabilities
 * @param totP sum of p
 * @return integer from 0 to length of p minus 1
 */
int sampleInt(std::vector<double> probs, double totP = 1)
{
  double u = R::runif(0, totP);
  double sum = probs[0];
  int i = 0;
  while (sum < u) {
    ++i;
    sum += probs[i];
  }
  return(i);
}

/**
 * sampleInt
 * @brief Method to sample integer using probabilities p
 * @param probs vector of probabilities 
 * @return integer from 0 to length of p minus 1
 */
int sampleInt(Eigen::VectorXd probs)
{
  double totP = probs.sum();
  double u = R::runif(0, totP);
  double sum = probs(0);
  int i = 0;
  while (sum < u) {
    ++i;
    sum += probs(i);
  }
  return(i);
}

/**
 * logPSplit
 * @brief log probability of tree split at `depth`: p_split(eta)=alpha/(1+d_eta)^beta
 * @param alpha parameter in range (0, 1)
 * @param beta parameter > 0
 * @param depth depth of split (begins at zero)
 * @param terminal if true returns log(1-p)
 * @return log probability
 */
double logPSplit(double alpha, double beta, int depth, bool terminal)
{
  double p = alpha * pow(1.0 + (double)depth, -beta);
  if (terminal) {
    return(log1p(-p));
  } else {
    return(log(p));
  }
}

/**
 * logDirichletDensity
 * @brief log probability of Dirichlet with values x and parameters alpha
 * @param x vector of values
 * @param alpha vector of parameters
 * @return log probability
 */
double logDirichletDensity(Eigen::VectorXd x, Eigen::VectorXd alpha)
{
  if (x.size() != alpha.size()) // ! incorrect sizes
    stop("logDirichletDensity incorrect size");
  double out = lgamma(alpha.sum());
  for (int i = 0; i < alpha.size(); ++i)
    out += ((alpha(i) - 1) * log(x(i))) - lgamma(alpha(i));
  return(out);
}

/**
 * rDirichlet
 * @brief random draw from Dirichlet distribution with parameters alpha
 * @param alpha parameters
 * @return vector containing draw from Dirichlet
 */
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

/**
 * @brief draw C^+(0, 1) full conditional using hierarchy: x^2|y~IG(1/2,1/y), y~IG(1/2,1). Full conditional: y|-~IG(1,x^2/(x^2+1)), x^2|-~IG((a+1)/2,1/(b/2+y))
 * 
 * @param x2 pointer to current value of parameter x^2
 * @param a additional IG component for x^2 full conditional
 * @param b additional IG component for x^2 full conditional
 * @param yInv pointer to update 1/y
 * @return double x^2 draw from full conditional
 */
void rHalfCauchyFC(double* x2, double a, double b, double* yInv)
{
  double yi = R::rgamma(1.0, *x2 / (*x2 + 1.0));
  if (yInv != 0)
    *yInv = yi;
  *x2 = 1.0 / R::rgamma(0.5 * (a + 1.0), 2.0 / (b + 2.0 * yi));
}

/**
 * intersectAndDiff
 * @brief Function to simultaneously calculate intersection and difference 
 of origVec to newVec
 * @param origVec starting vector of sorted integers
 * @param newVec vector of unsorted integers to be compared to origVec
 * @return std vector with 2 elements: intersection and difference
 */
std::vector<std::vector<int> > 
  intersectAndDiff(std::vector<int> origVec, std::vector<int> newVec)
{
  // Assume origVec is sorted, sort newVec
  std::sort(newVec.begin(), newVec.end());
  std::vector<int> intVec;
  std::vector<int> diffVec;
  std::vector<std::vector<int> > iD;
  
  if (origVec.size() == 0) {
    iD.push_back(origVec); // intersection
    iD.push_back(origVec); // difference
    return(iD);
  }
  if (newVec.size() == 0) {
    iD.push_back(newVec); // intersection
    iD.push_back(origVec);// difference
    return(iD);
  }
    
  intVec.reserve(origVec.size());
  diffVec.reserve(origVec.size());

  std::size_t i = 0;
  std::size_t j = 0;
  // iterate over origVec and newVec
  do {
    if (origVec[i] == newVec[j]) { // intersection
      intVec.push_back(origVec[i]);
      ++i;
      if (j < (newVec.size() - 1)) {
        ++j;
      }
      
    } else if (origVec[i] < newVec[j]) { // difference
      diffVec.push_back(origVec[i]);
      ++i;
      
    } else { // origVec[i] > newVec[i]
      if (j < (newVec.size() - 1)) {
        ++j;
      } else { // exceeded newVec, difference
        diffVec.push_back(origVec[i]);
        ++i;
      }
    }

  } while (i < origVec.size());

  iD.push_back(intVec);
  iD.push_back(diffVec);
  return(iD);
}

/**
 * @brief fast set intersection tool assumes sorted vectors A and B
 * 
 * @param A sorted integer vector A
 * @param B sorted integer vector B
 * @return vector of resulting intersection
 */
// [[Rcpp::export]]
std::vector<int> cppIntersection(const IntegerVector& A, 
                                 const IntegerVector& B) {
  std::vector<int> output;
  std::set_intersection(A.begin(), A.end(), B.begin(), B.end(),
                        std::back_inserter(output));
  return output;
}