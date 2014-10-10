/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balPLL.cpp
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
 * \file balPLL.cpp
 * \brief Implementation of the class PLL.
 */

#include <cmath>
#include "balPLL.h"

#ifndef WITHPHIERR
#define PLL_NEV 3
#else
#define PLL_NEV 0
#endif
#ifdef WITHPHIERR
// in this case, we have an additional state variable, PhiErr
#ifdef EXTEND
// we save the value of icp (the current in the charge pump)
// and divout (the output of the frequency divider)
#define PLL_NDIM 7
#else
#define PLL_NDIM 5
#endif
#else
#ifdef EXTEND
// we save the value of icp (the current in the charge pump)
// and divout (the output of the frequency divider)
#define PLL_NDIM 6
#else
// the number of state variables of the VCO
#define PLL_NDIM 4
#endif
#endif
#define PLL_NPAR 14

bal::DynamicalSystem* PLLFactory() {
  return new bal::PLL;
}

namespace bal {

const int PLL::npar = PLL_NPAR;
const char * PLL::parname[PLL_NPAR] = {"fref","r1","fvco","vdd","rho0",
        "rhoap","k0","krho","kap","alpha","kvcoa","kvcob","kvcoc","tuning"};

bool PLL::HasJacobian() const {
  return false;
}

bool PLL::HasEvents() const {
  return true;
}

bool PLL::HasEventsConstraints() const {
  return true;
}

PLL::PLL() : DynamicalSystem(PLL_NDIM, PLL_NPAR, PLL_NEV, false), pi(3.141592653589793) {
  // fixed parameters
  dt = 1e-9;
  tau_d = 0.0;
  C0 = 0.85e-9;
  C1 = 6e-9;
  Aud = 0.8e-3;
#ifdef MISMATCH
  // increase the delay in the output of the AND gate
  tau_d = 100e-9;
  // 10 percent mismatch between the current generators in the charge pump
  Aud_mismatch = 0.1;
#endif

#ifndef FRACTIONAL
  N = 2400;
#else
  for(int i=0; i<ndiv/2; i++) {
    N[i] = 2400;
    N[i+ndiv/2] = 2401;
  }
  idx = 0;
#endif

#ifndef WITHPHIERR
  divout = 0;
  cnt = 0;
#endif
  
  // initialize variables that need to be reset at every new integration
  Reset();
}

PLL::~PLL() {
}

int PLL::RHS (realtype t, N_Vector X, N_Vector Xdot, void *sys) {
  realtype x, y, r, w;
#ifdef WITHPHIERR
  realtype phierr;
#endif
  realtype icp;
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  Parameters *parameters = ds->GetParameters();

  fREF = parameters->At(0);
  T = 1.0/fREF;
  R1 = parameters->At(1);
  omega0 = 2*pi*parameters->At(2);
  Vdd = parameters->At(3);
  rho0 = parameters->At(4);
  rhoap = parameters->At(5);
  k0 = parameters->At(6);
  Krho = parameters->At(7);
  Kap = parameters->At(8);
  alpha = parameters->At(9);
  KVCOa = parameters->At(10);
  KVCOb = parameters->At(11);
  KVCOc = parameters->At(12);
  tuning_coeff = parameters->At(13);
  
  x = Ith (X, 0);
  y = Ith (X, 1);
  r = Ith (X, 2);
  w = Ith (X, 3);
#ifdef WITHPHIERR
  phierr = Ith (X, 4);
#endif
  
  realtype gamma = sqrt(std::max(x*x+y*y,1e-12));
  // saturation
  if(w > 3.0) {
    w = 3.0;
    Ith(X,3) = 3.0;
  }
  if(w < -3.0) {
    w = -3.0;
    Ith(X,3) = -3.0;
  }
  
  realtype vtune = w/tuning_coeff;
  
  realtype A = ((rho0+Krho*vtune)/gamma - 1) * k0;
  realtype B = ((1-alpha)*Kap*(gamma-rhoap) + 1 + alpha*vtune*(KVCOa + KVCOb*vtune + KVCOc*vtune*vtune))*omega0;
  Ith (Xdot, 0) = A*x - B*y;
  Ith (Xdot, 1) = A*y + B*x;

  // ODE for r
  Ith (Xdot, 2) =  1.0 / (R1*C1) * (w-r);
  
  // computation of icp, the current in the charge pump
#ifdef WITHPHIERR

  realtype phihat, that, tmp;
  tmp = fabs(phierr/(2*pi));
  phihat = (tmp - floor(tmp));
  tmp = t/T;
  that = (tmp - floor(tmp));
  S = (phihat > that ? true : false);
  icp = S*Aud*(fabs(phierr)/phierr);

#else

#ifndef MISMATCH
  icp = zu*Aud - zd*Aud;
#else
  realtype Au, Ad;
  Au = Aud*(1.0+Aud_mismatch);
  Ad = Aud;
  icp = zu*Au - zd*Ad;
#endif

#endif
  
  // ODE for w
  Ith (Xdot, 3) =  1.0 / C0 * (icp + (r-w)/R1);

#ifdef WITHPHIERR
  // ODE for phierr
#ifndef FRACTIONAL
  Ith (Xdot, 4) = 2*pi*fREF - (x*Ith(Xdot,1) - y*Ith(Xdot,0))/(N*gamma*gamma);
#else
  Ith (Xdot, 4) = 2*pi*fREF - (x*Ith(Xdot,1) - y*Ith(Xdot,0))/(N[idx]*gamma*gamma);
#endif
#endif
  
#ifdef EXTEND
  Ith (X, GetDimension()-2) = (realtype) divout;
  Ith (X, GetDimension()-1) = icp;
  Ith (Xdot, GetDimension()-2) = 0.;
  Ith (Xdot, GetDimension()-1) = 0.;
#endif
  
  return CV_SUCCESS;
}

#ifdef CVODE25
int PLL::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		      void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int PLL::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
		      void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  return ! CV_SUCCESS;
}

void PLL::Reset() {
#ifndef WITHPHIERR
  wait_reset = false;
  n = 1;
  treset = 0;
  zu = zd = false;
  fprintf(stderr, "%c%s", ESC, RED);
  fprintf(stderr, "zu = %d zd = %d @ 0\n", zu, zd);
  fprintf(stderr, "%c%s", ESC, NORMAL);
#else
  S = true;
#endif
}

int PLL::Events (realtype t, N_Vector X, realtype * event, void * data) {
#ifndef WITHPHIERR
  realtype x, y, r, w;
  Parameters * parameters = (Parameters *) data;
  
  x = Ith (X, 0);
  y = Ith (X, 1);
  r = Ith (X, 2);
  w = Ith (X, 3);
  
  fREF = parameters->At(0);
  T = 1.0/fREF;
  R1 = parameters->At(1);
  omega0 = 2*pi*parameters->At(2);
  Vdd = parameters->At(3);
  rho0 = parameters->At(4);
  rhoap = parameters->At(5);
  k0 = parameters->At(6);
  Krho = parameters->At(7);
  Kap = parameters->At(8);
  alpha = parameters->At(9);
  KVCOa = parameters->At(10);
  KVCOb = parameters->At(11);
  KVCOc = parameters->At(12);
  tuning_coeff = parameters->At(13);
  
  // rising edge of the clock
  event[0] = t - n*T;
  // extrema of x
  realtype gamma = sqrt(x*x+y*y);
  realtype vtune = w/tuning_coeff;
  event[1] = ((rho0+Krho*vtune)/gamma - 1)*k0*x - ((1-alpha)*Kap*(gamma-rhoap) + 1 + alpha*vtune*(KVCOa + KVCOb*vtune + KVCOc*vtune*vtune))*omega0*y;
  // reset
  event[2] = t - (treset+tau_d+dt);
#endif
  return CV_SUCCESS;
}
 
void PLL::EventsConstraints (realtype t, N_Vector X, int * constraints, void * data) {
#ifndef WITHPHIERR
  constraints[0] = 1;
  // minimum of x [we know that x is a sinusoid with mean 0]
  constraints[1] = (Ith(X,0) < 0 ? 1 : 0);
  constraints[2] = 1;
#endif
}

void PLL::ManageEvents(realtype t, N_Vector X, int * events, int * constraints) {
#ifndef WITHPHIERR
  // everything is off for t < T
  if(t < T) {
    return;
  }
  
  int c[3] = {1,1,1};
  if(constraints != NULL) {
    for(int i=0; i<3; i++) c[i] = constraints[i];
  }
  if(events[0]) {
    n++;
    zu = true;
    fprintf(stderr, "%c%s", ESC, CYAN);
    fprintf(stderr, "zu <= 1 @ %e\n", t);
    fprintf(stderr, "%c%s", ESC, NORMAL);
  }
  if(events[1] && c[1]) {
    cnt++;
#ifndef FRACTIONAL
    if((divout && cnt == floor((double) N/2)) || (!divout && cnt == ceil((double) N/2))) {
#else
    if((divout && cnt == floor((double) N[idx]/2)) || (!divout && cnt == ceil((double) N[idx]/2))) {
#endif
      cnt = 0;
      divout = !divout;
      if(divout) { // rising edge of the output of the frequency divider
	zd = true;
#ifdef FRACTIONAL
	idx = (idx+1) % ndiv;
#endif
	fprintf(stderr, "%c%s", ESC, YELLOW);
	fprintf(stderr, "zd <= 1 @ %e\n", t);
	fprintf(stderr, "%c%s", ESC, NORMAL);
      }
    }
  }
  if(events[2]) {
    zu = zd = false;
    fprintf(stderr, "%c%s", ESC, MAGENTA);
    fprintf(stderr, "zu <= zd <= 0 @ %e\tw = %e\n", t, Ith(X,3));
    fprintf(stderr, "%c%s", ESC, NORMAL);
    wait_reset = false;
  }
  if(zu && zd && ! wait_reset) {
    treset = t;
    wait_reset = true;
  }
}

DynamicalSystem* PLL::Clone() const {
  return new PLL(*this);
}

} // namespace bal

#endif

