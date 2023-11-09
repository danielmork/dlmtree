/**
 * This file implements the Polya-gamma sampler PG(1,z).
 * This is a C++ implementation of Algorithm 6 in PhD thesis of Jesse
 * Bennett Windle, 2013
 * URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
 *
 * References:
 *
 *   Jesse Bennett Windle
 *   Forecasting High-Dimensional, Time-Varying Variance-Covariance Matrices
 *   with High-Frequency Data and Sampling Polya-Gamma Random Variates for
 *   Posterior Distributions Derived from Logistic Likelihoods
 *   PhD Thesis, 2013
 *
 *   Damien, P. & Walker, S. G. Sampling Truncated Normal, Beta, and Gamma Densities
 *   Journal of Computational and Graphical Statistics, 2001, 10, 206-215
 *
 *   Chung, Y.: Simulation of truncated gamma variables
 *   Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
 *
 * (c) Copyright Enes Makalic and Daniel F Schmidt, 2018
 *
 * ! Changes for current work: Modified to work with RcppEigen
 */

#include "RcppEigen.h"

using namespace Rcpp;

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// FCN prototypes
double samplepg(double);
double exprnd(double);
double tinvgauss(double, double);
double truncgamma();
double randinvg(double);
double aterm(int, double, double);

/**
 * @brief multiple draw polya gamma latent variable for var c[i] with size b[i]
 * 
 * @param b vector of binomial sizes
 * @param c vector of parameters
 * @return Eigen::VectorXd 
 */
// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd rcpp_pgdraw(Eigen::VectorXd b,
                            Eigen::VectorXd c)
{
  int m = b.size();
  int n = c.size();
  Eigen::VectorXd y(n); y.setZero();

  // Setup
  int i, j, bi = 1;
  if (m == 1)
  {
    bi = b[0];
  }

  #pragma omp parallel
  #pragma omp for
  
  // Sample
  // TODO: add code for parallel draws: #pragma omp parallel for
  for (i = 0; i < n; i++){
    // Rcout << i << "\n";
    // Rcout << "First parameter: " << b[i] << "\n";
    if (m > 1){
      bi = b[i];
    }

    // Sample
    y[i] = 0;
    for (j = 0; j < (int)bi; j++){
      // Rcout << j << "\n";
      // Rcout << "Second parameter: " << c[j] << "\n";
      y[i] += samplepg(c[i]);
    }
  }

  return y;
}


/**
 * @brief single draw polya gamma latent variable for var c with size b
 * 
 * @param b double value of binomial sizes
 * @param c double value of parameter
 * @return double
 */
// [[Rcpp::depends(RcppEigen)]]
double rcpp_pgdraw(double b, double c){

  double y = 0;

  // Sample
  // TODO: add code for parallel draws: #pragma omp parallel for
  #pragma omp parallel
  #pragma omp for
  
  for (int i = 0; i < int(b); i++){
    y += samplepg(c);
  }
  
  return y;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = fabs(z) * 0.5;

  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;

  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = log(4) - MATH_LOG_PI - z;
  double logK = log(K);
  double Kt = K * t;
  double w = sqrt(MATH_PI_2);

  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = exp(logf1) + exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q);

  double u, X;

  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1)
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss(z, t);
    }

    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;

    while(1)
    {
      Sn = Sn + asgn * aterm(i, X, t);

      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }

      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }

      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// Generate exponential distribution random variates
double exprnd(double mu)
{
  return -mu * log(1.0 - R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + log(n + 0.5) + 1.5*(MATH_LOG_2_PI-log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }
  return exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - sqrt(4*mu*V + mu*mu * V*V) );

  if(R::runif(0.0,1.0) > mu /(mu+out)) {
    out = mu*mu / out;
  }
  return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / sqrt(X);

    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }

  return X;
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
// Note that mu is arbitrary constant
double tinvgauss(double z, double t){
  double X, u;
  double mu = 1.0/z;

  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();

      if(log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }
  return X;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------

const int grid_size = 81;
const double ygrid[] = {
  0.0625,0.06698584,0.07179365,0.07694653,0.08246924,
  0.08838835,0.09473229,0.1015315,0.1088188,0.1166291,
  0.125,0.1339717,0.1435873,0.1538931,0.1649385,
  0.1767767,0.1894646,0.2030631,0.2176376,0.2332582,
  0.25,0.2679434,0.2871746,0.3077861,0.329877,
  0.3535534,0.3789291,0.4061262,0.4352753,0.4665165,
  0.5,0.5358867,0.5743492,0.6155722,0.659754,
  0.7071068,0.7578583,0.8122524,0.8705506,0.933033,
  1,1.071773,1.148698,1.231144,1.319508,
  1.414214,1.515717,1.624505,1.741101,1.866066,
  2,2.143547,2.297397,2.462289,2.639016,
  2.828427,3.031433,3.24901,3.482202,3.732132,
  4,4.287094,4.594793,4.924578,5.278032,
  5.656854,6.062866,6.498019,6.964405,7.464264,
  8,8.574188,9.189587,9.849155,10.55606,
  11.31371,12.12573,12.99604,13.92881,14.92853,
  16};

const double vgrid[] = {
  -256,-222.8609,-194.0117,-168.897,-147.0334,
  -128,-111.4305,-97.00586,-84.4485,-73.51668,
  -63.99997,-55.71516,-48.50276,-42.22387,-36.75755,
  -31.99844,-27.85472,-24.24634,-21.10349,-18.36524,
  -15.97843,-13.89663,-12.07937,-10.49137,-9.101928,
  -7.884369,-6.815582,-5.875571,-5.047078,-4.315237,
  -3.667256,-3.092143,-2.580459,-2.124095,-1.716085,
  -1.350442,-1.022007,-0.7263359,-0.4595871,-0.2184366,
  0,0.1982309,0.3784427,0.5425468,0.6922181,
  0.828928,0.953973,1.068498,1.173516,1.269928,
  1.358533,1.440046,1.515105,1.584282,1.64809,
  1.706991,1.761401,1.811697,1.858218,1.901274,
  1.941143,1.978081,2.012318,2.044068,2.073521,
  2.100856,2.126234,2.149802,2.171696,2.192042,
  2.210954,2.228537,2.244889,2.260099,2.274249,
  2.287418,2.299673,2.311082,2.321703,2.331593,
  2.340804};

const double tol  = 1e-8;

// This is the code for faster PG sampling when h is a big integer in PG(h, z)
// This code is specifically written for negative binomial / zero-inflated version
// For the binary component of zero-inflated model, the modelEst code uses PG sampling code above
// For the negative binomial component, PG uses the saddle-point(sp) approach suggested by Windle 2013 (Algorithm 9)
// The reference is the same as written on top of this code.

struct Line {
  public:
    double slope;
    double itcpt;
};

struct FD {
  public:
    double val;
    double der;
};

// Functions from rcpp_pgdraw.cpp
double phi_mode(double);         // mode of phi = tanh(z) / z
double v_eval(double);            // alpha_l, alpha_r computing function
double tangent_to_eta(double, double, double, Line&); // Takes values and returns slope and intercept
double left_tgamma(double, double, double);


// Computes mode for the left point of the envelope
double phi_mode(double x){
  double tol = 1e-6; // This is fixed
  double y   = 0.0;
  double r   = sqrt(fabs(x));
  if (x > tol){
    y = tan(r) / r;
  }
  else if (x < -1*tol){
    y = tanh(r) / r;
  } else {
    y = 1 + (1/3) * x + (2/15) * x * x + (17/315) * x * x * x;
  }

  return y;
}

// Function for fdf_eval
void ydy_eval(double v, double* yp, double* dyp)
{
  double y = phi_mode(v);
  *yp = y;

  if (fabs(v) >= tol)
    *dyp = 0.5 * (y*y + (1-y) / v);
  else
    *dyp = 0.5 * (y*y - 1/3 - (2/15) * v);

}

// Function for v_eval
void fdf_eval(double v, void* params, double* fp, double* dfp){
  double y = *((double*)params);
  ydy_eval(v, fp, dfp);
  *fp  -= y;
}

// Lemma 2.24
double v_eval(double y){
  double ylower = ygrid[0];
  double yupper = ygrid[grid_size-1];

  if (y < ylower) {
    return -1. / (y*y);
  } else if (y > yupper) {
    double v = atan(0.5 * y * MATH_PI);
    return v*v;
  }

  else if (y==1) return 0.0;

  double id = (log(y) / log(2.0) + 4.0) / 0.1;

  int idlow  = (int)id;
  int idhigh = (int)id + 1;
  double vl  = vgrid[idlow];  // lower bound
  double vh  = vgrid[idhigh]; // upper bound

  double diff = tol + 1.0;
  double vnew = vl;
  double vold = vl;
  double f0, f1;

  while (diff > tol) {
    vold = vnew;
    fdf_eval(vold, &y, &f0, &f1);
    vnew = vold - f0 / f1;
    vnew = vnew > vh ? vh : vnew;
    vnew = vnew < vl ? vl : vnew;
    diff = fabs(vnew - vold);
  }

  return vnew; 
}

double cos_rt(double v){
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v >= 0)
    y = cos(r);
  else
    y = cosh(r);
  return y;
}

void delta_func(double x, double mid, FD& delta){
  if (x >= mid) {
    delta.val = log(x) - log(mid);
    delta.der = 1.0 / x;
  }
  else {
    delta.val = 0.5 * (1 - 1.0 / x) - 0.5 * (1 - 1.0 / mid);
    delta.der = 0.5 / (x*x);
  }
}

// function for tangent_to_eta
double phi_func(double x, double z, FD& phi)
{
  // double v = yv.v_func(x);
  double v = v_eval(x);
  double u = 0.5 * v;
  double t = u + 0.5 * z*z;

  phi.val = log(cosh(fabs(z))) - log(cos_rt(v)) - t * x;
  phi.der = -1.0 * t;

  return v;
}

double tangent_to_eta(double x, double z, double mid, Line& tl){
  FD phi, delta, eta;
  double v;

  v = phi_func(x, z, phi);
  delta_func(x, mid, delta);

  eta.val = phi.val - delta.val;
  eta.der = phi.der - delta.der;

  // printf("v=%g\nphi=%g, phi.d=%g\ndelta=%g, delta.d=%g\neta=%g, eta.d=%g\n",
  // 	 v, phi.val, phi.der, delta.val, delta.der, eta.val, eta.der);

  tl.slope = eta.der;
  tl.itcpt = eta.val - eta.der * x;
  
  return v;
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
// Note that mu is pre-specified
double tinvgauss_sp(double mu, double z, double t){
  double X, u;

  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();

      if(log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }
  return X;
}

double p_igauss(double x, double mu, double lambda)
{
    // z = 1 / mean
    double z = 1 / mu;
    double b = sqrt(lambda / x) * (x * z - 1);
    double a = sqrt(lambda / x) * (x * z + 1) * -1.0;
    double y = R::pnorm(b, 0.0, 1.0, 1, 0) + exp(2 * lambda * z) * R::pnorm(a, 0.0, 1.0, 1, 0);
    return y;
}

double ltgamma(double shape, double rate, double trunc){
  double a = shape;
  double b = rate * trunc;

  if (trunc <=0) {
      fprintf(stderr, "ltgamma: trunc = %g < 0\n", trunc);
      return 0;
  }
  if (shape < 1) {
      fprintf(stderr, "ltgamma: shape = %g < 1\n", shape);
      return 0;
  }

  if (shape == 1) {
    return R::rexp(1) / rate + trunc;
  }

  double d1 = b-a;
  double d3 = a-1;
  double c0 = 0.5 * (d1 + sqrt(d1*d1 + 4 * b)) / b;

  double x = 0.0;
  bool accept = false;

  while (!accept) {
      x = b + R::rexp(1) / c0;
      double u = R::runif(0, 1);

      double l_rho = d3 * log(x) - x * (1-c0);
      double l_M   = d3 * log(d3 / (1-c0)) - d3;

      accept = log(u) <= (l_rho - l_M);
  }

  return trunc * (x/b);
}