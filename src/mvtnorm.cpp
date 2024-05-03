/**
 * @file mvtnorm.cpp
 * @author Daniel Mork (danielmork.github.io)
 * @brief Functions for truncated multivariate normal density and sampling
 * @version 1.0
 * @date 2021-03-05
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// make fortran function accessible
extern "C" {
  extern void mvtdst_(int* n,
                      int* nu,
                      double* lower,
                      double* upper,
                      int* infin,
                      double* correl,
                      double* delta,
                      int* maxpts,
                      double* abseps,
                      double* releps,
                      double* error,
                      double* value,
                      int* inform);
}


//' Integrates (0,inf) over multivariate normal 
//'
//' @param mu vector of mean parameters
//' @param sigma covariance matrix
//' @returns double 
//' @export
// [[Rcpp::export]]
double zeroToInfNormCDF(Eigen::VectorXd mu, Eigen::MatrixXd sigma) {
  // Returns log P(X>0) for X~MVN(mean,sigma)
  int n = mu.size();

  // univariate
  if (n == 1) {
    return(R::pnorm5(0, mu(0), sqrt(sigma(0,0)), false, false));

    // multivariate
  } else {
    // vars for Fortran code
    int nu_         = 0;
    int maxpts_     = n * 1000;
    double abseps_  = 0.0001;
    double releps_  = 0;

    double* lower   = new double[n];
    double* upper   = new double[n];
    int* infin      = new int[n];
    double* delta   = new double[n];
    double* corrTri = new double[n * (n-1) / 2];

    // fill bounds, calculate lower tri correlation matrix
    int i, j;
    int k = 0;
    for (i = 0; i < n; ++i) {
      lower[i] = -mu(i) / sqrt(sigma(i, i));
      upper[i] = 0.0;
      infin[i] = 1; // 1 indicates [lower, +inf)
      delta[i] = 0.0;

      if (i > 0) {
        for (j = 0; j < i; ++j) {
          corrTri[k] = sigma(j, i) / sqrt(sigma(i, i) * sigma(j, j));
          ++k;
        }
      }
    }

    double error_ = 0.0;
    double value_ = 0.0;
    int inform_   = 0;

    mvtdst_(&n, &nu_, lower, upper, infin, corrTri, delta,
            &maxpts_, &abseps_, &releps_, &error_, &value_, &inform_);
    //Rcout << value_ << "\n" << error_ << "\n" << inform_;
    return(value_);
  }
}


// Code for truncated multivariate normal sampler
// adapted from R packages tmvmixnorm and bnmr.
// Algorithm is described in detail in paper:
// "Efficient sampling methods for truncated multivariate
//    normal and student-t distributions subject to linear inequality
//    constraints", by Yifang Li and Sujit Ghosh
//    J. of Statistical Theory and Practice 9:4 (2015), pp. 712-732.

/**
 * @brief univariate sampling of truncated [a, infinity)
 * 
 * @param a lower
 * @returns double 
 */
double rtnorm1(double a) {
  double x = 0.0;
  double y;

  if (a < 0) { // normal rejection sampling
    while (x == 0) {
      y = R::rnorm(0, 1);
      if (y > a)
        return(y);
    }
  } else if (a < 0.25696) { // half-normal rejection sampling
    while (x == 0) {
      y = fabs(R::rnorm(0, 1));
      if (y > a)
        return(y);
    }
  } else { // one-sided translated-exponential rejection sampling
    while (x == 0) {
      double lambdastar = (a + sqrt(a * a + 4.0)) / 2.0;
      y = R::rexp(1)/lambdastar + a;
      if (R::runif(0, 1) < exp(-0.5 * pow(y, 2.0) + lambdastar * y -
          0.5 * lambdastar + log(lambdastar)))
        return(y);
    }
  }
  return(x);
}

/**
 * @brief univariate sampling of truncated [a, b]; a < 0, b > 0
 * 
 * @param a lower
 * @param b upper
 * @returns double 
 */
double rtnorm2(double a, double b) {
  double x = 0.0;
  double y;

  if (b > a + sqrt(2 * M_PI)) { // normal rejection sampling
    while (x==0) {
      y = R::rnorm(0, 1);
      if ((y > a) && (y < b))
        return(y);
    }
  } else { // uniform rejection sampling
    while (x == 0) {
      y = R::runif(a, b);
      if (R::runif(0, 1) < exp(-pow(y, 2) / 2))
        return(y);
    }
  }
  return(x);
}

/**
 * @brief univariate sampling of truncated [a, b]; a > 0
 * 
 * @param a lower
 * @param b upper
 * @returns double 
 */
double rtnorm3(double a, double b) {
  double x = 0.0;
  double y;
  double a2 = pow(a, 2);

  if (a < 0.25696) {

    if (b > (a + sqrt(M_PI / 2) * exp(a2 / 2))) {
      // half-normal rejection sampling
      while(x==0){
        y = fabs(R::rnorm(0, 1));
        if((y > a) && (y < b))
          return(y);
      }
    } else { // uniform rejection sampling
      while (x == 0) {
        y = R::runif(a, b);
        if (R::runif(0, 1) < exp((a2 - pow(y, 2)) / 2))
          return(y);
      }
    }

  } else { // if (a > 0.25696)

    if (b > (a + 2/(a + sqrt(a2 + 4)) * exp((a2 - a*sqrt(a2 + 4))/4 + 0.5))) {
      // two-sided translated-exponential rejection sampling
      while (x == 0) {
        double lambdastar = (a + sqrt(a2 + 4)) / 2;
        y =  a - log(R::runif(exp((a - b) * lambdastar), 1)) / lambdastar;
        if (R::runif(0,1) <
          (exp(-pow(y, 2)/2 + lambdastar*y - lambdastar/2 + log(lambdastar))))
          return(y);
      }
    } else { // uniform rejection sampling
      while (x == 0) {
        y = R::runif(a, b);
        if (R::runif(0, 1) < (exp((a2-pow(y, 2))/2)))
          return(y);
      }
    }
  } // end if (a > 0.25696)
  return(x);
}

/**
 * @brief truncated univarite normal sampler, mean 0, sd 1
 * 
 * @param a lower truncation
 * @param b upper trancation
 * @returns double 
 */
double rtuvnorm(double a, double b) {
  double x = 0.0;

  if (!std::isfinite(a)) { // (-Inf, ?]

    if (!std::isfinite(b)) { // No truncation
      return(R::rnorm(0, 1));

    } else { // Case 1: truncated (-Inf, b]
      return(-rtnorm1(-b));
    }

  } else { // [a, ?]

    if (!std::isfinite(b)) { // Case 2: truncated [a, Inf)
      return(rtnorm1(a));

    } else if (a < 0) { // [a, b]

      if (b > 0) { // Case 3: truncated [a, b]; a < 0; b > 0
        return(rtnorm2(a, b));

      } else { // Case 4: truncated [a, b]; b < 0
        return(-rtnorm3(-b, -a));
      }

    } else { // Case 5: truncated [a, b]; a >= 0
      return(rtnorm3(a, b));
    }
  }
  return(x);
}


//' Truncated multivariate normal sampler, mean mu, cov sigma, truncated (0, Inf)
//'
//' @param mu vector of mean parameters
//' @param sigma covariance matrix
//' @param iter number of iterations
//' @returns VectorXd
//' @export
// [[Rcpp::export]]
Eigen::VectorXd rtmvnorm(Eigen::VectorXd mu, Eigen::MatrixXd sigma, int iter) 
{
  int n = mu.size();
  VectorXd z(n); z.setZero();
  const VectorXd a = -mu;
  const MatrixXd R = sigma.llt().matrixL();
  
  if (n == 1) {
    z(0) = rtuvnorm(-mu(0) / R(0, 0), INFINITY);
    return(R * z + mu);
  }
  
  z = R.triangularView<Eigen::Lower>().solve((a.array() + 1).matrix());
  int i, j, k;
  double amin, amax;
  VectorXd atmp = a - R * z;
  
  for (i = 0; i < iter; ++i) {  // # of iterations
    for (j = 0; j < n; ++j) { // loop over multivariate normal
      amin = -INFINITY;
      amax = INFINITY;
      atmp += R.col(j) * z(j);
      for (k = j; k < n; ++k) {
        if (R(k, j) > 0)
          amin = std::max(amin, atmp(k) / R(k, j));
        else
          amax = std::min(amax, atmp(k) / R(k, j));
      }
      z(j) = rtuvnorm(amin, amax);
      atmp -= R.col(j) * z(j);
    }
  }
  return(R * z + mu);
}
