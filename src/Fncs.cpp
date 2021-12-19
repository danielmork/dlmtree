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
    return(log1p(-p)); // log(1 + (-p))
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
 * @brief randomly draw a value from zero-truncated normal (Utilizes inverse transform method)
 * 
 * @param n desired number of draws
 * @param mu mean of zero-truncated normal
 * @param sigma standard deviation of zero-truncated normal
 * @param lower lower bound of zero-truncated normal 
 * @param upper upper bound of zero-truncated normal 
 * @return a double vector from from zero-truncated normal
 */
std::vector<double> rtruncnorm(int n, double mu, double sigma, double lower, double upper){
  if(n < 1){
    stop("sampling size must be greater than 0");
  }

  if(sigma < 0){
    stop("standard deviation must be positive");
  }

  if(lower > upper){
    stop("lower bound value must be smaller than the upper bound value");
  }

  std::vector<double> sample;
  for(int i = 0; i < n; i++){
    double u = R::runif(0, 1); // Sample from a standard uniform
    double transformed = (R::qnorm(R::pnorm(lower, 0, 1, 1, 0) + u*(R::pnorm(upper, 0, 1, 1, 0) - R::pnorm(lower, 0, 1, 1, 0)), 0, 1, 1, 0))*sigma + mu;
    sample.push_back(transformed);
  }

  return sample;
}

/**
 * @brief randomly draw a SINGLE value from zero-truncated normal (Utilizes inverse transform method)
 * 
 * @param mu mean of zero-truncated normal
 * @param sigma standard deviation of zero-truncated normal
 * @param lower lower bound of zero-truncated normal 
 * @param upper upper bound of zero-truncated normal 
 * @return a double vector from from zero-truncated normal
 */
double rtruncnorm(double mu, double sigma, double lower, double upper){
  if(sigma < 0){
    stop("standard deviation must be positive");
  }

  if(!(lower < upper)){
    stop("lower bound value must be smaller than the upper bound value");
  }

  double u = R::runif(0, 1); // Sample from a standard uniform
  double transformed = (R::qnorm(R::pnorm(lower, 0, 1, 1, 0) + u*(R::pnorm(upper, 0, 1, 1, 0) - R::pnorm(lower, 0, 1, 1, 0)), 0, 1, 1, 0))*sigma + mu;

  if(transformed < lower || upper < transformed){
    stop("The sampled value is not in the support");
  }

  return transformed;
}


/**
 * @brief evaluates zero-truncated normal density
 * 
 * @param x a value to be evaluated
 * @param mu mean of zero-truncated normal
 * @param sigma standard deviation of zero-truncated normal
 * @param lower lower bound of zero-truncated normal 
 * @param upper upper bound of zero-truncated normal 
 * @return a double value: probability
 */
double dtruncnorm(double x, double mu, double sigma,  double lower, double upper){
  if(x < lower || upper < x){
    stop("The value to be evaluated must be between the lower and the upper bound");
  }

  if(sigma < 0){
    stop("standard deviation must be positive");
  }

  if(!(lower < upper)){
    stop("lower bound value must be smaller than the upper bound value");
  }

  double xz = (x - mu)/sigma; // _z notation for "scaled"
  double lowerz = (lower - mu)/sigma; // _z notation for "scaled"
  double upperz = (upper - mu)/sigma; // _z notation for "scaled"

  double numer = R::dnorm(xz, 0, 1, 0)/sigma;
  double denom = R::pnorm(upperz, 0, 1, 1, 0) - R::pnorm(lowerz, 0, 1, 1, 0);

  double prob = numer / denom;

  return prob;
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

/**
 * @brief Subset a vector only with given indices
 * 
 * @param original A vector to be subset
 * @param indices A vector containing wanted indices
 * @return A vector with values of given indices
 */

Eigen::VectorXd selectInd(Eigen::VectorXd original, std::vector<int> indices) {

  int m = indices.size(); // Get the total number of index

  Eigen::VectorXd subset; // define a subset
  subset.resize(m);

  // For loop to collect values with a corresponding index
  for(int i = 0; i < m; i++){
    int index = indices[i]; // Get an index from indices vector
    double val = original(index); // Find the value corresponding to the index
    subset(i) = val; // Save the value
  }

  return subset;
}

/**
 * @brief Subset a matrix only with given indices (rows)
 * 
 * @param original A vector to be subset
 * @param indices A vector containing wanted indices of row
 * @return A vector with values of given indices
 */

Eigen::MatrixXd selectIndM(Eigen::MatrixXd original, std::vector<int> indices) {

  int rownum = indices.size(); // row# of submat
  int colnum = original.cols(); // col# of submat

  Eigen::MatrixXd submat; // define a submat
  submat.resize(rownum, colnum);

  // For loop to collect values with a corresponding index
  for(int i = 0; i < rownum; i++){
    int index = indices[i]; // Get an index from indices vector
    for(int j = 0; j < colnum; j++){
      double val = original(index, j); // Find the value corresponding to the index
      submat(i, j) = val;              // Save the value
    }
  }

  return submat;
}