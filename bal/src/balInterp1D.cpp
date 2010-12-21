/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp1D.cpp
 *
 *   Copyright (C) 2009,2010 Daniele Linaro
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *=========================================================================*/

/**
 *  \file balInterp1D.cpp
 *  \brief Implementation of classes balBaseInterp1D, balLinearInterp1D, balPolyInterp1D and balSplineInterp1D. 
 */

#include "balInterp1D.h"

/***** balBaseInterp1D *****/

const char * balBaseInterp1D::GetClassName() const {
  return "balBaseInterp1D";
}

void balBaseInterp1D::Destroy() {
  delete this;
}

int balBaseInterp1D::Locate(const double x) {
  int ju, jm, jl, ascnd;
  
  if (n < 2 || mm < 2 || mm > n) {
    throw("locate size error");
  }
  
  ascnd = (xx[n-1] >= xx[0]); // 1 if ascending order of table, 0 otherwise. 
  jl = 0;   // Initialize lower 
  ju = n-1; // and upper limits. 
  while (ju-jl > 1) {  // If we are not yet done, 
    jm = (ju+jl) >> 1; // compute a midpoint, 
    if ((x >= xx[jm]) == ascnd) 
      jl = jm; // and replace either the lower limit 
    else 
      ju = jm; // or the upper limit, as appropriate. 
  } // Repeat until the test condition is satisﬁed. 
  return MAX(0,MIN(n-mm,jl-((mm-2)>>1))); 
}

int balBaseInterp1D::Hunt(const double x) {
  int jl=jsav, jm, ju, inc=1;
  if (n < 2 || mm < 2 || mm > n) 
    throw("hunt size error");
  bool ascnd=(xx[n-1] >= xx[0]);
  if (jl < 0 || jl > n-1) {
    jl=0;
    ju=n-1; 
  }
  else {
    if ((x >= xx[jl]) == ascnd) { 
      for (;;) {
	ju = jl + inc; 
	if (ju >= n-1) { 
	  ju = n-1; break;
	}
	else if ((x < xx[ju]) == ascnd) 
	  break;
	else { 
	  jl = ju;
	  inc += inc;
	}
      } 
    } 
    else {
      ju = jl;
      for (;;) {
	jl = jl - inc;
	if (jl <= 0) {
	  jl = 0;
	  break;
	} 
	else if ((x >= xx[jl]) == ascnd)
	  break;
	else {	
	  ju = jl;
	  inc += inc;
	}
      }
    }
  }
  
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1;
    if ((x >= xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  cor = fabs(jl-jsav) > dj ? 0 : 1; jsav = jl;
  return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

bool balBaseInterp1D::nextHunt() {
  if (cor)
    return true;
  else
    return false;
}

balBaseInterp1D::balBaseInterp1D()
  : n(0), mm(0), jsav(0), cor(0) { 
  dj = MIN(1, (int)pow((double)n,0.25));
  xx = NULL;
  yy = NULL;
  nnd = 1;
}

void balBaseInterp1D::SetInterpolationPoints(double * xi, double **yi, int length, int nf) {
  n = length;
  nnf = nf;
  xx = xi; 
  yy = yi;  
}

balBaseInterp1D::balBaseInterp1D(const balBaseInterp1D & interp) : balInterpolator(interp) {
  n = interp.n;
  mm = interp.mm;
  xx = interp.xx;
  yy = interp.yy;
  jsav = 0;
  cor = 0;
  dj = MIN(1, (int)pow((double)n,0.25));
}

balBaseInterp1D::~balBaseInterp1D() {
}

/***** balLinearInterp1D *****/

const char * balLinearInterp1D::GetClassName() const {
  return "balLinearInterp1D";
}

void balLinearInterp1D::Destroy() {
  delete this;
}

balLinearInterp1D * balLinearInterp1D::Create() {
  return new balLinearInterp1D();
}

balLinearInterp1D * balLinearInterp1D::Copy(balLinearInterp1D *interp) {
  return new balLinearInterp1D(*interp);
}

balLinearInterp1D * balLinearInterp1D::Clone() const {
  return new balLinearInterp1D(*this);
}

balLinearInterp1D::balLinearInterp1D() :
  balBaseInterp1D() {
  mm = 2;
}

balLinearInterp1D::balLinearInterp1D(const balLinearInterp1D & interp) :
  balBaseInterp1D(interp) {}

balLinearInterp1D::~balLinearInterp1D() {}

int balLinearInterp1D::Evaluate(double *x, double *y) {
  
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balLinearInterp1D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  
  int i;
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);
 
  if (xx[jlo]==xx[jlo+1]) {
    //The xa’s must be distinct. 
    cerr<<"Bad input to balLinearInterp1D::Evaluate()\n";
    return -1;
  }
  for (i=0; i<nnf; i++)
    y[i] = yy[i][jlo] + ((x[0]-xx[jlo])/(xx[jlo+1]-xx[jlo]))*(yy[i][jlo+1]-yy[i][jlo]);
  return 0;
}
  
int balLinearInterp1D::EvaluateJacobian(double *x, double **y) {
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balLinearInterp1D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  int i;
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);
  if (xx[jlo]==xx[jlo+1])  {
    //The xa’s must be distinct. 
    cerr<<"Bad input to balLinearInterp1D::EvaluateJacobian()\n";
    return -1;
  }
  for (i=0; i<nnf; i++)
    y[i][0] = (yy[i][jlo+1]-yy[i][jlo])/(xx[jlo+1]-xx[jlo]);
  return 0;
}

int balLinearInterp1D::EvaluateDivergence(double *x, double *y) {
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balLinearInterp1D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 1) {
    cerr<<"balLinearInterp1D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 1\n";
    return -1;
  }
  int i;
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);
  if (xx[jlo]==xx[jlo+1])  {
    //The xa’s must be distinct. 
    cerr<<"Bad input to balLinearInterp1D::EvaluateDivergence()\n";
    return -1;
  }
  y[0] = (yy[0][jlo+1]-yy[0][jlo])/(xx[jlo+1]-xx[jlo]);
  return 0;

}
/***** balPolyInterp1D *****/

const char * balPolyInterp1D::GetClassName() const {
  return "balPolyInterp1D";
}

void balPolyInterp1D::Destroy() {
  delete this;
}

balPolyInterp1D * balPolyInterp1D::Create() {
  return new balPolyInterp1D();
}

balPolyInterp1D * balPolyInterp1D::Copy(balPolyInterp1D *interp) {
  return new balPolyInterp1D(*interp);
}

balPolyInterp1D * balPolyInterp1D::Clone() const {
  return new balPolyInterp1D(*this);
}

balPolyInterp1D::balPolyInterp1D() :
  balBaseInterp1D() {
  dy = 0.0;
  mm = 2;
}

void balPolyInterp1D::SetInterpolationOrder(int m) {
  mm = m;
}

balPolyInterp1D::balPolyInterp1D(const balPolyInterp1D & interp) :
  balBaseInterp1D(interp) {
  dy = interp.dy;
}

balPolyInterp1D::~balPolyInterp1D() {}


int balPolyInterp1D::Evaluate(double *x, double *y) {
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balPolyInterp1D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);
  int i, j, m, ns; 
  double *den, dif, dift, ho, hp, *w; 
  double *xa;
  double **c, **d;
  
  c = new double * [mm];
  for (i=0; i<mm; i++) 
    c[i] = new double[nnf];
  d = new double * [mm];
  for (i=0; i<mm; i++) 
    d[i] = new double[nnf];
  w = new double[nnf];
  den = new double[nnf];

  ns = 0;
  xa = &xx[jlo];
  
  dif = fabs(x[0]-xa[0]); 
  for (i=0;i<mm;i++) { // Here we ﬁnd the index ns of the closest table entry, 
    if ((dift=fabs(x[0]-xa[i])) < dif) { 
      ns=i; 
      dif=dift; 
    } 
    for (j=0; j<nnf; j++) {
      c[i][j] = yy[j][i+jlo]; // and initialize the tableau of c’s and d’s. 
      d[i][j] = yy[j][i+jlo]; 
    }
  } 
  for (j=0; j<nnf; j++)
    y[j] = yy[j][ns+jlo]; // This is the initial approximation to y. 
  ns--;
  for (m=1; m<mm; m++) { // For each column of the tableau,
    for (i=0;i<mm-m;i++) { // we loop over the current c’s and d’s and update them.
      ho = xa[i]-x[0]; 
      hp = xa[i+m]-x[0]; 
      for (j=0; j<nnf; j++) {
        w[j] = c[i+1][j]-d[i][j];
        den[j] = ho-hp;
      }       
      if (den[0] == 0.0) { // This error can occur only if two input xa’s are (to within roundoﬀ ) identical. 
	for (j=0; j<mm; j++) {
          delete [] c[j];
	  delete [] d[j];
	}
	delete [] c;
	delete [] d;
	delete [] w;
	delete [] den;
	throw("PolyInterp1D error");
      }
      for (j=0; j<nnf; j++) {
	den[j] = w[j]/den[j]; 
	d[i][j] = hp*den[j]; // Here the c’s and d’s are updated. 
	c[i][j] = ho*den[j]; 
      }
    } 
    if (2*(ns+1) < (mm-m))
      for (j=0; j<nnf; j++)
	y[j] += c[ns+1][j];
    else {
      for (j=0; j<nnf; j++) {
	y[j] += d[ns][j];
      }
      ns--;
    }
  }
  for (j=0; j<mm; j++) {
    delete [] c[j];
    delete [] d[j];
  }
  delete [] c;
  delete [] d;
  delete [] w;
  delete [] den;
  return 0;
}

int balPolyInterp1D::EvaluateJacobian(double *x, double **y) { 
  // Implemented with finite differences
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balPolyInterp1D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  double x2 = *x+FINITE_DIFFERENCES_STEP; 
  double y2[nnf], y1[nnf];
  Evaluate(x,y1);
  Evaluate(&x2,y2);
  for (int i=0; i<nnf; i++) {
      y[i][0] = (y2[i]-y1[i])/FINITE_DIFFERENCES_STEP;
  }
    
  return 0;
}

int balPolyInterp1D::EvaluateDivergence(double *x, double *y) { 
  // Implemented with finite differences
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balPolyInterp1D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 1) {
    cerr<<"balPolyInterp1D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 1\n";
    return -1;
  }
  double x2 = *x+FINITE_DIFFERENCES_STEP; 
  double y2[nnf], y1[nnf];
  Evaluate(x,y1);
  Evaluate(&x2,y2);
  y[0] = (y2[0]-y1[0])/FINITE_DIFFERENCES_STEP;
    
  return 0;
}

/***** balSplineInterp1D *****/

const char * balSplineInterp1D::GetClassName() const {
  return "balSplineInterp1D";
}

void balSplineInterp1D::Destroy() {
  delete this;
}

balSplineInterp1D * balSplineInterp1D::Create() {
  return new balSplineInterp1D();
}

balSplineInterp1D * balSplineInterp1D::Copy(balSplineInterp1D *interp) {
  return new balSplineInterp1D(*interp);
}
  
balSplineInterp1D * balSplineInterp1D::Clone() const {
  return new balSplineInterp1D(*this);
}

balSplineInterp1D::balSplineInterp1D() :
  balBaseInterp1D() {
  yyp1 = NATURAL_SPLINE;
  yypn = NATURAL_SPLINE;
  y2 = NULL;
  mm = 2;
}

balSplineInterp1D::balSplineInterp1D(const balSplineInterp1D & interp) :
  balBaseInterp1D(interp) {
  if (interp.y2!=NULL) {
    y2 = new double * [n];
//    memcpy(y2,interp.y2,n*sizeof(double));
    for (int i=0; i<n; i++) { 
      y2[i] = new double[nnf];
      memcpy(y2[i],interp.y2[i],nnf*sizeof(double));
    }
  }
  else
    y2 = NULL;
}

balSplineInterp1D::~balSplineInterp1D() { 
  if (y2!=NULL) {
    for (int i=0; i<n; i++) 
      delete [] y2[i];
    delete [] y2;
  }
}

int balSplineInterp1D::Init() {
  if ((xx == NULL)||(yy == NULL)) {
    cerr<<"balSplineInterp1D::Init() - Interpolation points not set\n";
    return -1;
  }
  if (y2!=NULL) {
    for (int i=0; i<n; i++) 
      delete [] y2[i];
    delete [] y2;
  }
  y2 = new double * [n];
  for (int i=0; i<n; i++) 
    y2[i] = new double[nnf];
  balSplineInterp1D::Sety2();
  return 0;
}

void balSplineInterp1D::SetBoundaryConditions(double yp1, double ypn) {
  yyp1 = yp1;
  yypn = ypn;
}

void balSplineInterp1D::Sety2() {
  int i, j, k; 
  double p, qn, sig, *un; 
  double **u; 
  
  u = new double * [n-1];
  for (j=0; j<n-1; j++)
    u[j] = new double[nnf];
  un = new double[nnf];
  
  if (yyp1 > 0.99e99)	// The lower boundary condition is set either to be ``natural'' 
    for (j=0; j<nnf; j++)
      y2[0][j] = u[0][j]= 0.0;
  else {							// or else to have a speciﬁed ﬁrst derivative. 
    for (j=0; j<nnf; j++) {
      y2[0][j] = -0.5; 
      u[0][j] = (3.0/(xx[1]-xx[0]))*((yy[j][1]-yy[j][0])/(xx[1]-xx[0])-yyp1); 
    }
	} 
  for (i=1;i<n-1;i++) { 
    /*
     * This is the decomposition loop of the tridiagonal algorithm. 
     * y2 and u are used for temporary storage of the decomposed 
     * factors. 
     */
    sig = (xx[i]-xx[i-1])/(xx[i+1]-xx[i-1]); 
    
    for (j=0; j<nnf; j++) {
      p = sig*y2[i-1][j]+2.0; 
      y2[i][j] = (sig-1.0)/p; 
      u[i][j] = (yy[j][i+1]-yy[j][i])/(xx[i+1]-xx[i]) - (yy[j][i]-yy[j][i-1])/(xx[i]-xx[i-1]); 
      u[i][j] = (6.0*u[i][j]/(xx[i+1]-xx[i-1])-sig*u[i-1][j])/p; 
    }
  } 
  if (yypn > 0.99e99) {	// The upper boundary condition is set either to be ``natural''
    qn = 0.0;
    for (j=0; j<nnf; j++)
      un[j] = 0.0;
  }
  else {							// or else to have a speciﬁed ﬁrst derivative. 
    qn = 0.5; 
    for (j=0; j<nnf; j++) 
      un[j] = (3.0/(xx[n-1]-xx[n-2]))*(yypn-(yy[j][n-1]-yy[j][n-2])/(xx[n-1]-xx[n-2])); 
  } 
  for (j=0; j<nnf; j++) {
    y2[n-1][j] = (un[j]-qn*u[n-2][j])/(qn*y2[n-2][j]+1.0); 
    // This is the backsubstitution loop of the tridiagonal algorithm. 
    for (k=n-2;k>=0;k--) 
      y2[k][j]=y2[k][j]*y2[k+1][j]+u[k][j];
  }
  for (i=0; i<n-1; i++)
    delete [] u[i];
  delete [] u;
  delete [] un;
}

int balSplineInterp1D::Evaluate(double *x, double *y) {
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balSplineInterp1D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  if (y2 == NULL) {
    cerr<<"balSplineInterp1D::Evaluate() - Spline not initialized. Call method Init()\n";
    return -1;
  }
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);
  int i;
  int klo, khi; 
  double h, b, a;
  
  klo = jlo;
  khi = jlo+1;
  h = xx[khi]-xx[klo];
  
  if (h == 0.0) {
    //The xa’s must be distinct. 
    cerr<<"Bad input to balSplineInterp1D::Evaluate()\n";
    return -1;
  }
  a = (xx[khi]-x[0])/h;
  b = (x[0]-xx[klo])/h; // Cubic spline polynomial is now evaluated. 
  for (i=0; i<nnf; i++)
    y[i] = a*yy[i][klo] + b*yy[i][khi] + ((a*a*a-a)*y2[klo][i] + (b*b*b-b)*y2[khi][i])*(h*h)/6.0; 
  return 0;
}

int balSplineInterp1D::EvaluateJacobian(double *x, double **y) { 
  
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balSplineInterp1D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  
  if (y2 == NULL) {
    cerr<<"balSplineInterp1D::EvaluateJacobian() - Spline not initialized. Call method Init()\n";
    return -1;
  }
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);
  int i;
  int klo,khi;
  double h, b, a, adot, bdot;

  klo = jlo;
  khi = jlo+1;
  h = xx[khi]-xx[klo];
  
  if (h == 0.0) {
    //The xa’s must be distinct. 
    cerr<<"Bad input to balSplineInterp1D::EvaluateJacobian()\n";
    return -1;
  }

  a = (xx[khi]-x[0])/h;
  b = (x[0]-xx[klo])/h;
  adot = -1/h;
  bdot = 1/h;
  for (i=0; i<nnf; i++)
    y[i][0] = adot*yy[i][klo]+bdot*yy[i][khi]+((3*a*a*adot-adot)*y2[klo][i]+(3*b*b*bdot-bdot)*y2[khi][i])*(h*h)/6.0; 
  return 0;
}

int balSplineInterp1D::EvaluateDivergence(double *x, double *y) { 
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balSplineInterp1D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 1) {
    cerr<<"balSplineInterp1D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 1\n";
    return -1;
  }
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);
  int i;
  int klo,khi;
  double h, b, a, adot, bdot;

  klo = jlo;
  khi = jlo+1;
  h = xx[khi]-xx[klo];
  
  if (h == 0.0) {
    //The xa’s must be distinct. 
    cerr<<"Bad input to balSplineInterp1D::EvaluateDivergence()\n";
    return -1;
  }

  a = (xx[khi]-x[0])/h;
  b = (x[0]-xx[klo])/h;
  adot = -1/h;
  bdot = 1/h;
  y[0] = adot*yy[0][klo]+bdot*yy[0][khi]+((3*a*a*adot-adot)*y2[klo][i]+(3*b*b*bdot-bdot)*y2[khi][0])*(h*h)/6.0; 
    
  return 0;
}
/***** balSmoothingSplineInterp1D *****/

const char * balSmoothingSplineInterp1D::GetClassName() const {
  return "balSmoothingSplineInterp1D";
}

void balSmoothingSplineInterp1D::Destroy() {
  delete this;
}

balSmoothingSplineInterp1D * balSmoothingSplineInterp1D::Create() {
  return new balSmoothingSplineInterp1D();
}

balSmoothingSplineInterp1D * balSmoothingSplineInterp1D::Copy(balSmoothingSplineInterp1D *interp) {
  return new balSmoothingSplineInterp1D(*interp);
}

balSmoothingSplineInterp1D * balSmoothingSplineInterp1D::Clone() const {
  return new balSmoothingSplineInterp1D(*this);
}
  
balSmoothingSplineInterp1D::balSmoothingSplineInterp1D() :
  balBaseInterp1D() {
  mm = 2;
  SS = 0;
  a = NULL;
  b = NULL;
  c = NULL;
  d = NULL;
  ddy = NULL;
  _DEALLOC_DDY = 0;
}

int balSmoothingSplineInterp1D::Init() {
 if ((xx==NULL)||(yy==NULL)) {
   cerr<<"balSmoothingSplineInterp1D::Init() - Interpolation points not set\n";
   return -1;
 }
  int i;
  if (a!=NULL)
    for (i=0; i<nnf; i++)
      delete [] a[i];
    delete [] a;
  if (b!=NULL)
    for (i=0; i<nnf; i++)
      delete [] b[i];
    delete [] b; 
  if (c!=NULL)
    for (i=0; i<nnf; i++)
      delete [] c[i];
    delete [] c; 
  if (d!=NULL)
    for (i=0; i<nnf; i++)
      delete [] d[i];
    delete [] d; 
  a = new double * [nnf];
  b = new double * [nnf];
  c = new double * [nnf];
  d = new double * [nnf];
  for (int i=0; i<nnf; i++) {
    a[i] = new double[n];
    b[i] = new double[n];
    c[i] = new double[n];
    d[i] = new double[n];
    balSmoothingSplineInterp1D::ComputeCoefficients(i);
  }
  return 0;
}
 
int balSmoothingSplineInterp1D::ComputeCoefficients(int index) {

  if (SS == 0) {
    ddy = new double [n];
    for (int i=0; i<n; i++) 
      ddy[i] = 1.0;
  }

  int i;
  double *r,*r1,*r2,*t,*t1,*u,*v,p,h,f,f2,g,e;

  r = new double[n];
  r1 = new double[n+1];
  r2 = new double[n+2];
  t = new double[n-2];
  t1 = new double[n-2];
  u = new double[n+2];
  v = new double[n+2];

  r[0] = 0.0; 
  r[1] = 0.0; 
  r1[n] = 0.0; 
  r2[n] = 0.0; 
  r2[n+1] = 0.0; 
  u[0] = 0.0; 
  u[1] = 0.0; 
  u[n] = 0.0; 
  u[n+1] = 0.0; 
  p = 0.0;

  h = xx[1]-xx[0];
  f = (yy[index][1]-yy[index][0])/h;
 
  for (i=2; i<=n-1; i++) {
    g = h;
    h = xx[i]-xx[i-1];
    e = f;
    f = (yy[index][i]-yy[index][i-1])/h;
    a[index][i-1] = f-e;
    t[i-2] = 2*(g+h)/3.0;
    t1[i-2] = h/3.0;
    r2[i] = ddy[i-2]/g;
    r[i] = ddy[i]/h;
    r1[i] = -ddy[i-1]/g-ddy[i-1]/h;
  }

  for (i=2; i<=n-1; i++) {
    b[index][i-1] = r[i]*r[i]+r1[i]*r1[i]+r2[i]*r2[i];
    c[index][i-1] = r[i]*r1[i+1]+r1[i]*r2[i+1];
    d[index][i-1] = r[i]*r2[i+2];
  }
  f2 = -SS;

  // Next iteration
  while(1) {
    for (i=2; i<=n-1; i++) {
      r1[i-1] = f*r[i-1];
      r2[i-2] = g*r[i-2];
      r[i] = 1/(p*b[index][i-1]+t[i-2]-f*r1[i-1]-g*r2[i-2]);
      u[i] = a[index][i-1]-r1[i-1]*u[i-1]-r2[i-2]*u[i-2];
      f = p*c[index][i-1]+t1[i-2]-h*r1[i-1];
      g = h;
      h = d[index][i-1]*p;
    }
  
    for (i=n-1; i>=2; i--) {
      u[i] = r[i]*u[i]-r1[i]*u[i+1]-r2[i]*u[i+2];
    }
    e = 0.0;
    h = 0.0;
    
    for (i=1;i<=n-1;i++) {
      g = h;
      h = (u[i+1]-u[i])/(xx[i]-xx[i-1]);
      v[i] = (h-g)*ddy[i-1]*ddy[i-1];
      e = e+v[i]*(h-g);
    }
  
    g = -h*ddy[n-1]*ddy[n-1];
    v[n] = -h*ddy[n-1]*ddy[n-1];
    e = e-g*h;
    g = f2;
    f2 = e*p*p;
  
    if ((f2>=f)||(f2<=g)) {
      for (i=1;i<=n;i++) {
	a[index][i-1] = yy[index][i-1]-p*v[i];
	c[index][i-1] = u[i];
      }
      for (i=1;i<=n-1;i++) {
	h = xx[i]-xx[i-1];
	d[index][i-1] = (c[index][i]-c[index][i-1])/(3.0*h);
	b[index][i-1] = (a[index][i]-a[index][i-1])/h-(h*d[index][i-1]+c[index][i-1])*h;
      }

      delete [] r;
      delete [] r1;
      delete [] r2;
      delete [] t;
      delete [] t1;
      delete [] u;
      delete [] v;
      if (SS == 0) {
        delete [] ddy;
	ddy = NULL;
      }
      return 0;
    }

    f = 0.0;
    h = (v[2]-v[1])/(xx[1]-xx[0]);
    for (i=2;i<=n-1;i++) {
      g = h;
      h = (v[i+1]-v[i])/(xx[i]-xx[i-1]);
      g = h-g-r1[i-1]*r[i-1]-r2[i-2]*r[i-2];
      f = f+g*r[i]*g;
      r[i] = g;
    }
    h = e-p*f;
    if (h<=0) {
      for (i=1;i<=n;i++) {
	a[index][i-1] = yy[index][i-1]-p*v[i];
	c[index][i-1] = u[i];
      }
      for (i=1;i<=n-1;i++) {
	h = xx[i]-xx[i-1];
	d[index][i-1] = (c[index][i]-c[index][i-1])/(3*h);
	b[index][i-1] = (a[index][i]-a[index][i-1])/h-(h*d[index][i-1]+c[index][i-1])*h;
      }
      delete [] r;
      delete [] r1;
      delete [] r2;
      delete [] t;
      delete [] t1;
      delete [] u;
      delete [] v;
      if (SS == 0) {
        delete [] ddy;
	ddy = NULL;
      }
      return 0;
    } 
     	p = p+(SS-f2)/((sqrt(SS/e)+p)*h);
  }
}

balSmoothingSplineInterp1D::balSmoothingSplineInterp1D(const balSmoothingSplineInterp1D & interp) :
  balBaseInterp1D(interp) {
  int i;
  if (interp.a!=NULL) {
    a = new double * [nnf];
    //memcpy(a,interp.a,nnf*sizeof(double));
    for (i=0; i<nnf; i++) {
      a[i] = new double[n];
      memcpy(a[i],interp.a[i],n*sizeof(double));
    } 
  }
  else
    a = NULL;
  if (interp.b!=NULL) {
    b = new double * [nnf];
    //memcpy(b,interp.b,nnf*sizeof(double));
    for (i=0; i<nnf; i++) {
      b[i] = new double[n];
      memcpy(b[i],interp.b[i],n*sizeof(double));
    } 
  }
  else
    b = NULL;
  if (interp.c!=NULL) {
    c = new double * [nnf];
    //memcpy(c,interp.c,nnf*sizeof(double));
    for (i=0; i<nnf; i++) {
      c[i] = new double[n];
      memcpy(c[i],interp.c[i],n*sizeof(double));
    } 
  }
  else
    c = NULL;
  if (interp.d!=NULL) {
    d = new double * [nnf];
    //memcpy(d,interp.d,nnf*sizeof(double));
    for (i=0; i<nnf; i++) {
      d[i] = new double[n];
      memcpy(d[i],interp.d[i],n*sizeof(double));
    } 
  }
  else
    d = NULL;
  ddy = interp.ddy;
}

balSmoothingSplineInterp1D::~balSmoothingSplineInterp1D() { 
  int i;
  if (a!=NULL)
    for (i=0; i<nnf; i++)
      delete [] a[i];
    delete [] a;
  if (b!=NULL)
    for (i=0; i<nnf; i++)
      delete [] b[i];
    delete [] b; 
  if (c!=NULL)
    for (i=0; i<nnf; i++)
      delete [] c[i];
    delete [] c; 
  if (d!=NULL)
    for (i=0; i<nnf; i++)
      delete [] d[i];
    delete [] d; 
  if (_DEALLOC_DDY)
    delete [] ddy;

}

void balSmoothingSplineInterp1D::SetSmoothingParameters(double *dy, double S) {
  SS = S;
  if (_DEALLOC_DDY)
    delete [] ddy;
  ddy = dy;
  _DEALLOC_DDY = 0;
}

void balSmoothingSplineInterp1D::SetSmoothingParameters(double smooth) {
  SS = 1;
  if (_DEALLOC_DDY)
    delete [] ddy;
  ddy = new double[n];
  for (int i=0; i<n; i++)
    ddy[i] = smooth;
  _DEALLOC_DDY = 1;
}

int balSmoothingSplineInterp1D::Evaluate(double *x, double *y) {
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balSmoothingSplineInterp1D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  if ((a == NULL) || (b == NULL) || (c == NULL) || (d == NULL)) {
    cerr<<"balSmoothingSplineInterp1D::Evaluate() - Spline not initilized. Call method Init(). \n";
    return -1;
  }
  int i;
  double h;
  
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);

  h = x[0]-xx[jlo];
  
  for (i=0; i<nnf; i++)
    y[i] = ((d[i][jlo]*h+c[i][jlo])*h+b[i][jlo])*h+a[i][jlo];
  return 0;
}

int balSmoothingSplineInterp1D::EvaluateJacobian(double *x, double **y) { 
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balSmoothingSplineInterp1D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  if ((a == NULL) || (b == NULL) || (c == NULL) || (d == NULL)) {
    cerr<<"balSmoothingSplineInterp1D::EvaluateJacobian() - Spline not initilized. Call method Init(). \n";
    return -1;
  }
  int i;
  double h;
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);

  h = x[0]-xx[jlo];
  
  for (i=0; i<nnf; i++)
    y[i][0] = (3*d[i][jlo]*h+2*c[i][jlo])*h+b[i][jlo];
  return 0;
}

int balSmoothingSplineInterp1D::EvaluateDivergence(double *x, double *y) { 
  if ((xx == NULL) || (yy == NULL)) {
    cerr<<"balSmoothingSplineInterp1D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 1) {
    cerr<<"balSmoothingSplineInterp1D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 1\n";
    return -1;
  }
  int i;
  double h;
  int jlo = cor ? Hunt(x[0]) : Locate(x[0]);

  h = x[0]-xx[jlo];
  
  y[0] = (3*d[0][jlo]*h+2*c[0][jlo])*h+b[0][jlo];
  return 0;
}
