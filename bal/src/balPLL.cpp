/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balPLL.cpp
 *
 *   Copyright (C) 2009 Daniele Linaro
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
 * \brief Implementation of the class balPLL
 */

#include "balPLL.h"

const int balPLL::npar = 14;
const char * balPLL::parname[14] = {"fref","r1","fvco","vdd","rho0",
																	 "rhoap","k0","krho","kap","alpha",
																	 "kvcoa","kvcob","kvcoc","tuning"};

balDynamicalSystem* balPLLFactory() {
	return balPLL::Create();
}

balPLL::balPLL() : pi(3.141592653589793) {
#ifndef WITHPHIERR
	SetNumberOfEvents(3);
#else
	SetNumberOfEvents(0);
#endif
	
	// the number of state variables of the VCO
	int ndim = 4;

#ifdef WITHPHIERR
	// in this case, we have an additional state variable, PhiErr
	ndim++;
#endif

	SetDimension(ndim);
	SetNumberOfParameters(npar);

	// fixed parameters
	dt = 1e-9;
	C0 = 0.85e-9;
	C1 = 6e-9;
	Aud = 0.8e-3;
	N = 2400;
#ifndef WITHPHIERR
	divout = 0;
	cnt = 0;
#endif

	// initialize variables that need to be reset at every new integration
	Reset();
}

balPLL::~balPLL() {}

int balPLL::RHS (realtype t, N_Vector X, N_Vector Xdot, void * data) {
  realtype x, y, r, w;
#ifdef WITHPHIERR
	realtype phierr;
#endif
	realtype icp;
	balParameters * parameters = (balParameters *) data;

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

	realtype gamma = sqrt(max(x*x+y*y,1e-12));
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

#ifdef WITHPHIERR
	Ith (Xdot, 4) = 2*pi*fREF - (x*Ith(Xdot,1) - y*Ith(Xdot,0))/(N*gamma*gamma);
#endif

	// ODE for r
  Ith (Xdot, 2) =  1.0 / (R1*C1) * (w-r);

	// compute icp, the current in the charge pump
#ifdef WITHPHIERR
	realtype phihat, that, tmp;
	tmp = fabs(phierr/(2*pi));
	phihat = (tmp - floor(tmp));
	tmp = t/T;
	that = (tmp - floor(tmp));
	S = (phihat > that ? true : false);
	icp = S*Aud*(fabs(phierr)/phierr);
#else
	icp = zu*Aud - zd*Aud;
#endif

	// ODE for w
	Ith (Xdot, 3) =  1.0 / C0 * (icp + (r-w)/R1);

	return CV_SUCCESS;
}

#ifdef CVODE25
int balPLL::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int balPLL::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
				void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
	return ! CV_SUCCESS;
}

void balPLL::Reset() {
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

int balPLL::Events (realtype t, N_Vector X, realtype * event, void * data) {
#ifndef WITHPHIERR
	realtype x, y, r, w;
	balParameters * parameters = (balParameters *) data;

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
	event[2] = t - (treset+dt);
#endif
	return CV_SUCCESS;
}

void balPLL::EventsConstraints (realtype t, N_Vector X, int * constraints, void * data) {
#ifndef WITHPHIERR
  constraints[0] = 1;
	// minimum of x [we know that x is a sinusoid with mean 0]
  constraints[1] = (Ith(X,0) < 0 ? 1 : 0);
  constraints[2] = 1;
#endif
}

void balPLL::ManageEvents(realtype t, N_Vector X, int * events, int * constraints) {
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
		if(cnt >= N/2) {
			cnt = 0;
			divout = !divout;
			if(divout) {
				zd = true;
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
#endif
}


