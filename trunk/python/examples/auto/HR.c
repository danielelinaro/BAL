#include "auto_f2c.h"

/* Common Block Declarations */
extern struct {
	  integer itwist, istart, iequib, nfixed, npsi, nunstab, nstab;
} blhom_1;

/* Hindmarsh-Rose AUTO file */

int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp) {
  doublereal x, y, z;
  doublereal b, I, mi, s;
  doublereal xrest = -1.6;
  integer dfdu_dim1, dfdp_dim1;

  /* Evaluates the algebraic equations or ODE right hand side */
  
  /* Input arguments : */
  /*      ndim   :   Dimension of the ODE system */
  /*      u      :   State variables */
  /*      icp    :   Array indicating the free parameter(s) */
  /*      par    :   Equation parameters */
  
  /* Values to be returned : */
  /*      f      :   ODE right hand side values */
  
  
  x = u[0]; y = u[1]; z = u[2];
  b = par[0]; I = par[1]; mi = par[2]; s = par[3];

  f[0] = y - x*x*x + b*x*x + I - z;
  f[1] = 1 - 5*x*x - y;
  f[2] = mi*(s*(x-xrest) - z);

  if(ijac == 0) {
    return 0;
  }

  dfdu_dim1 = ndim;
  dfdp_dim1 = ndim;

  ARRAY2D(dfdu,0,0) = -3*x*x + 2*b*x;
  ARRAY2D(dfdu,0,1) = (doublereal)  1.;
  ARRAY2D(dfdu,0,2) = (doublereal) -1.;
  ARRAY2D(dfdu,1,0) = -10*x;
  ARRAY2D(dfdu,1,1) = (doublereal) -1.;
  ARRAY2D(dfdu,1,2) = (doublereal)  0.;
  ARRAY2D(dfdu,2,0) = mi*s;
  ARRAY2D(dfdu,2,1) = (doublereal)  0.;
  ARRAY2D(dfdu,2,2) = -mi;

  if(ijac == 1) {
    return 0;
  }

  ARRAY2D(dfdp,0,0) = x*x;
  ARRAY2D(dfdp,0,1) = (doublereal)  1.;
  ARRAY2D(dfdp,0,2) = (doublereal)  0.;
  ARRAY2D(dfdp,0,3) = (doublereal)  0.;
  ARRAY2D(dfdp,1,0) = (doublereal)  0.;
  ARRAY2D(dfdp,1,1) = (doublereal)  0.;
  ARRAY2D(dfdp,1,2) = (doublereal)  0.;
  ARRAY2D(dfdp,1,3) = (doublereal)  0.;
  ARRAY2D(dfdp,2,0) = (doublereal)  0.;
  ARRAY2D(dfdp,2,1) = (doublereal)  0.;
  ARRAY2D(dfdp,2,2) = s*(x-xrest) - z;
  ARRAY2D(dfdp,2,3) = mi*(x-xrest);

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par) {
  /* Input arguments : */
  /*      ndim   :   Dimension of the ODE system */
  
  /* Values to be returned : */
  /*      u      :   A starting solution vector */
  /*      par    :   The corresponding equation-parameter values */
  
  
  /* Initialize the equation parameters */
  
  par[0] = (doublereal)3.;
  par[1] = (doublereal)2.8;
  par[2] = (doublereal)1.0e-2;
  par[3] = (doublereal)4.;

  return 0;
}


/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc) {
  return 0;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint) {
    return 0;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp) {
    return 0;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par) {
  return 0;
}





