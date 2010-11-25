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
  
double balBaseInterp1D::interp(double x) {
  int jlo = cor ? hunt(x) : locate(x);
  return rawinterp(jlo,x);
}

int balBaseInterp1D::locate(const double x) {
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

int balBaseInterp1D::hunt(const double x) {
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

balBaseInterp1D::balBaseInterp1D(double * x, const double *y, int length, int m)
  : n(length), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) { 
  dj = MIN(1, (int)pow((double)n,0.25));
}

balBaseInterp1D::balBaseInterp1D(const balBaseInterp1D & interp) {
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

balLinearInterp1D * balLinearInterp1D::Create(double * xv, double * yv, int length) {
  return new balLinearInterp1D(xv,yv,length);
}

balLinearInterp1D * balLinearInterp1D::Copy(balLinearInterp1D *interp) {
  return new balLinearInterp1D(*interp);
}

balLinearInterp1D::balLinearInterp1D(double * xv, double * yv, int length) :
  balBaseInterp1D(xv,yv,length,2) {}

balLinearInterp1D::balLinearInterp1D(const balLinearInterp1D & interp) :
  balBaseInterp1D(interp) {}

balLinearInterp1D::~balLinearInterp1D() {}
  
double balLinearInterp1D::rawinterp(int j, double x) {
  if (xx[j]==xx[j+1]) 
    return yy[j];
  return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
}

/***** balPolyInterp1D *****/

const char * balPolyInterp1D::GetClassName() const {
  return "balPolyInterp1D";
}

void balPolyInterp1D::Destroy() {
  delete this;
}

balPolyInterp1D * balPolyInterp1D::Create(double * xv, double * yv, int length, int m) {
  return new balPolyInterp1D(xv,yv,length,m);
}

balPolyInterp1D * balPolyInterp1D::Copy(balPolyInterp1D *interp) {
  return new balPolyInterp1D(*interp);
}

balPolyInterp1D::balPolyInterp1D(double * xv, double * yv, int length, int m) :
  balBaseInterp1D(xv,yv,length,m), dy(0.0) {}

balPolyInterp1D::balPolyInterp1D(const balPolyInterp1D & interp) :
  balBaseInterp1D(interp) {
  dy = interp.dy;
}

balPolyInterp1D::~balPolyInterp1D() {}

double balPolyInterp1D::rawinterp(int j, double x) {
  int i, m, ns; 
  double y, den, dif, dift, ho, hp, w; 
  const double *xa, *ya;
  double *c, *d;
  
  c = new double[mm];
  if(c == NULL)
    throw("Memory allocation error");
  d = new double[mm];
  if(d == NULL) {
    delete c;
    throw("Memory allocation error");
  }
  ns = 0;
  xa = &xx[j];
  ya = &yy[j];
  
  dif = fabs(x-xa[0]); 
  for (i=0;i<mm;i++) { // Here we ﬁnd the index ns of the closest table entry, 
    if ((dift=fabs(x-xa[i])) < dif) { 
      ns=i; 
      dif=dift; 
    } 
    c[i] = ya[i]; // and initialize the tableau of c’s and d’s. 
    d[i] = ya[i]; 
  } 
  y = ya[ns--]; // This is the initial approximation to y. 
  for (m=1; m<mm; m++) { // For each column of the tableau,
    for (i=0;i<mm-m;i++) { // we loop over the current c’s and d’s and update them.
      ho = xa[i]-x; 
      hp = xa[i+m]-x; 
      w = c[i+1]-d[i]; 
      if ((den=ho-hp) == 0.0) { // This error can occur only if two input xa’s are (to within roundoﬀ ) identical. 
	delete c;
	delete d;
	throw("PolyInterp1D error");
      }
      den = w/den; 
      d[i] = hp*den; // Here the c’s and d’s are updated. 
      c[i] = ho*den; 
    } 
    y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--])); 
  }
  delete c;
  delete d;
  return y;
}

/***** balSplineInterp1D *****/

const char * balSplineInterp1D::GetClassName() const {
  return "balSplineInterp1D";
}

void balSplineInterp1D::Destroy() {
  delete this;
}

balSplineInterp1D * balSplineInterp1D::Create(double * xv, double * yv, int length, double yp1, double ypn) {
  return new balSplineInterp1D(xv,yv,length,yp1,ypn);
}

balSplineInterp1D * balSplineInterp1D::Copy(balSplineInterp1D *interp) {
  return new balSplineInterp1D(*interp);
}
  
balSplineInterp1D::balSplineInterp1D(double * xv, double * yv, int length, double yp1, double ypn) :
  balBaseInterp1D(xv,yv,length,2) {
  y2 = new double[length];
  sety2(xv,yv,yp1,ypn);
}

balSplineInterp1D::balSplineInterp1D(const balSplineInterp1D & interp) :
  balBaseInterp1D(interp) {
  y2 = new double[n];
  //sety2(xv,yv,yp1,ypn);
  memcpy(y2,interp.y2,n*sizeof(double));
}

balSplineInterp1D::~balSplineInterp1D() { delete y2; }

void balSplineInterp1D::sety2(double *xv, double *yv, double yp1, double ypn) {
  int i, k; 
  double p, qn, sig, un; 
  double *u; 
  
  u = new double[n-1];
  if(u == NULL)
    throw("Memory allocation failure");
  
  if (yp1 > 0.99e99)	// The lower boundary condition is set either to be ``natural'' 
    y2[0] = u[0]= 0.0;
  else {							// or else to have a speciﬁed ﬁrst derivative. 
    y2[0] = -0.5; 
    u[0] = (3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1); 
	} 
  for (i=1;i<n-1;i++) { 
    /*
     * This is the decomposition loop of the tridiagonal algorithm. 
     * y2 and u are used for temporary storage of the decomposed 
     * factors. 
     */
    sig = (xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]); 
    p = sig*y2[i-1]+2.0; 
    y2[i] = (sig-1.0)/p; 
    u[i] = (yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]); 
    u[i] = (6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p; 
  } 
  if (ypn > 0.99e99)	// The upper boundary condition is set either to be ``natural''
    qn = un = 0.0; 
  else {							// or else to have a speciﬁed ﬁrst derivative. 
    qn = 0.5; 
    un = (3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2])); 
  } 
  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0); 
  // This is the backsubstitution loop of the tridiagonal algorithm. 
  for (k=n-2;k>=0;k--) 
    y2[k]=y2[k]*y2[k+1]+u[k];
  
  delete u;
}

double balSplineInterp1D::rawinterp(int j, double x) {
  int klo, khi; 
  double y, h, b, a;
  
  klo = j;
  khi = j+1;
  h = xx[khi]-xx[klo];
  
  if (h == 0.0) {
    //The xa’s must be distinct. 
    throw("Bad input to balSplineInterp1D::rawinterp");
  }
  a = (xx[khi]-x)/h;
  b = (x-xx[klo])/h; // Cubic spline polynomial is now evaluated. 
  y = a*yy[klo] + b*yy[khi] + ((a*a*a-a)*y2[klo] + (b*b*b-b)*y2[khi])*(h*h)/6.0; 
  return y;
}

