/* $Id: C_FORTRAN_interface.c 357 2020-02-07 14:14:01Z thothorn $
  *
  *  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include "mvtnorm.h"

double F77_SUB(unifrnd)(void) { return unif_rand(); }
double F77_SUB(sqrtqchisqint)(int *n, double *p) {
    return(sqrt(qchisq(p[0], (double) n[0], 0, 0)));
}
double F77_SUB(phid)(double *x){ return pnorm(*x, 0.0, 1.0, 1, 0); }
double F77_SUB(studnt)(int *nu, double *x){ return pt(x[0], (double) nu[0], 1, 0); }

double F77_SUB(mvphi)(double const *z){
  return pnorm5(*z, 0., 1., 1L, 0L);
}

double F77_SUB(mvphnv)(double const *p){
  return qnorm5(*p, 0., 1., 1L, 0L);
}
