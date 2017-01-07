/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balDynamicalSystem.cpp
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
 * \file balDynamicalSystem.cpp
 * \brief Implementation of the class DynamicalSystem
 */

#include "balDynamicalSystem.h"
#include <iostream>

namespace bal {

DynamicalSystem::DynamicalSystem(int ndim, int npar, int nevents, bool extended) {
#ifdef DEBUG
  std::cout << "DynamicalSystem constructor.\n";
#endif
  if (ndim <= 0)
    throw "DynamicalSystem::DynamicalSystem - number of dimensions must be positive";
  if (npar < 0)
    throw "DynamicalSystem::DynamicalSystem - number of parameters cannot be negative";
  if (nevents < 0)
    throw "DynamicalSystem::DynamicalSystem - number of events cannot be negative";
  n = ndim;
  p = npar;
  pars = NULL;
  nExt = n*(n+1);
  nev = nevents;
  ext = extended;
#ifdef CVODE25
  jac = newDenseMat(n,n);
#endif
#ifdef CVODE26
  jac = NewDenseMat(n,n);
#endif
}

DynamicalSystem::DynamicalSystem(int ndim, Parameters* params, int nevents, bool extended) : 
        DynamicalSystem(ndim, params->GetNumber(), nevents, extended) {
  pars = params;
}

DynamicalSystem::DynamicalSystem(const DynamicalSystem& system) : pars(system.pars) {
#ifdef DEBUG
  std::cout << "DynamicalSystem copy constructor.\n";
#endif
  n = system.n;
  p = system.p;
  nev = system.nev;
  nExt = system.nExt;
  ext = system.ext;
#ifdef CVODE25
  jac = newDenseMat(n,n);
#endif
#ifdef CVODE26
  jac = NewDenseMat(n,n);
#endif
}

DynamicalSystem::~DynamicalSystem() {
#ifdef DEBUG
  std::cout << "DynamicalSystem destructor.\n";
#endif
#ifdef CVODE25
  destroyMat(jac);
#endif
#ifdef CVODE26
  DestroyMat(jac);
#endif
}

int DynamicalSystem::RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void *sys) {
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);

  // the first n components of the vector field
  int flag = ds->RHS(t, x, xdot, sys);

  if(! ds->IsExtended() || flag != CV_SUCCESS) {
    return flag;
  }
  
  // the Jacobian matrix
  if(ds->HasJacobian()) {
#ifdef CVODE25
    DynamicalSystem::JacobianWrapper(ds->n, ds->jac, t, x, NULL, sys, NULL, NULL, NULL);
#endif
#ifdef CVODE26
    DynamicalSystem::JacobianWrapper(ds->n, t, x, NULL, ds->jac, sys, NULL, NULL, NULL);
#endif
  }
  else {
    DynamicalSystem::JacobianFiniteDifferences(ds->n, t, x, ds->jac, sys);
  }

  realtype Y[ds->n][ds->n], F[ds->n][ds->n];
  int i, j, k;
  
  // extend
  for(i=0, k=0; i<ds->n; i++) {
    for(j=0; j<ds->n; j++, k++)
      Y[j][i] = Ith(x,3+k);
  }
  // multiply the matrices
  for(i=0; i<ds->n; i++) {
    for(j=0; j<ds->n; j++) {
      F[i][j] = 0.;
      for(k=0; k<ds->n; k++)
	F[i][j] += IJth(ds->jac,i,k)*Y[k][j];
    }
  }
  // the last n*n components of the vector field
  for(i=0, k=0; i<ds->n; i++) {
    for(j=0; j<ds->n; j++,k++)
      Ith(xdot,3+k) = F[j][i];
  }
  
  return CV_SUCCESS;
}

#ifdef CVODE25
int DynamicalSystem::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, 
			       N_Vector fy, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int DynamicalSystem::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, 
				  DlsMat J, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  return ! CV_SUCCESS;
}

#ifdef CVODE25
int DynamicalSystem::JacobianWrapper (long int N, DenseMat J, realtype t, N_Vector x, 
					 N_Vector fy, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int DynamicalSystem::JacobianWrapper (int N, realtype t, N_Vector x, N_Vector fy, 
					 DlsMat J, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
#ifdef CVODE25
  return ds->Jacobian(N, J, t, x, fy, sys, tmp1, tmp2, tmp3);
#endif
#ifdef CVODE26
  return ds->Jacobian(N, t, x, fy, J, sys, tmp1, tmp2, tmp3);
#endif
}

#ifdef CVODE25
int DynamicalSystem::JacobianFiniteDifferences (long int N, realtype t, N_Vector x, DenseMat J, void *sys) {
#endif
#ifdef CVODE26
int DynamicalSystem::JacobianFiniteDifferences (long int N, realtype t, N_Vector x, DlsMat J, void *sys) {
#endif
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  N_Vector ref, perturb;
  double eps = 1e-8;
  int i, j;

  ref = N_VNew_Serial(N);
  perturb = N_VNew_Serial(N);

  ds->RHS(t, x, ref, sys);

  for(j=0; j<N; j++) {
    Ith(x,j) = Ith(x,j) + eps;
    ds->RHS(t, x, perturb, sys);
    for(i=0; i<N; i++)
      IJth(J,i,j) = (Ith(perturb,i)-Ith(ref,i)) / eps;
    Ith(x,j) = Ith(x,j) - eps;
  }

  N_VDestroy_Serial(perturb);
  N_VDestroy_Serial(ref);

  return CV_SUCCESS;
}
  

int DynamicalSystem::Events (realtype t, N_Vector x, realtype *event, void *sys) {
  return ! CV_SUCCESS;
}

int DynamicalSystem::EventsWrapper (realtype t, N_Vector x, realtype *event, void *sys) {
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  return ds->Events(t, x, event, sys);
}

void DynamicalSystem::EventsConstraints (realtype t, N_Vector x, int *constraints, void *sys) {
}

Parameters* DynamicalSystem::GetParameters () const {
  return pars;
}

void DynamicalSystem::SetParameters(Parameters* params) {
  if (params->GetNumber() < p) {
    throw "The number of parameters passed is smaller than that of this system.";
  }
  else if (params->GetNumber() != p) {
    fprintf(stderr, "The number of parameters passed is greater than that of this system.\n");
  }
  pars = params;
}

int DynamicalSystem::GetDimension() const {
  return (ext ? nExt : n);
}

int DynamicalSystem::GetOriginalDimension() const {
  return n;
}

void DynamicalSystem::Extend(bool extend) {
  ext = extend;
}

bool DynamicalSystem::IsExtended() const {
  return ext;
}

bool DynamicalSystem::SpecialOptions(const void *opt) {
  std::cout << "DynamicalSystem::SpecialOptions\n";
  return false;
}

bool DynamicalSystem::HasJacobian() const {
  return false;
}
  
bool DynamicalSystem::HasEvents() const {
  return false;
}

bool DynamicalSystem::HasEventsConstraints() const {
  return false;
}
  
void DynamicalSystem::Reset() {
}

void DynamicalSystem::ManageEvents(realtype t, N_Vector X, int *events, int *constraints) {
}
  
int DynamicalSystem::GetNumberOfEvents() const {
  return nev;
}

int DynamicalSystem::GetNumberOfParameters() const {
  return p;
}
  
void DynamicalSystem::SetNumberOfEvents(int nev_) {
  if(nev_ >= 0)
    nev = nev_;
}

} // namespace bal

