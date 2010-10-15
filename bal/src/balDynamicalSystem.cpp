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
 * \brief Implementation of the class balDynamicalSystem
 */

#include "balDynamicalSystem.h"

balDynamicalSystem::balDynamicalSystem() {
  n = 0;
  p = 0;
  nev = 0;
  ext = false;
  nExt = 0;
  pars = NULL;
  jac = NULL;
  _dealloc = false;
}

balDynamicalSystem::balDynamicalSystem(const balDynamicalSystem& system) {
  n = system.n;
  nev = system.nev;
  pars = system.pars; // non rialloco spazio parametri
  p = (system.GetParameters())->GetNumber();
  nExt = system.nExt;
  ext = system.ext;
#ifdef CVODE25
  jac = newDenseMat(n,n);
#endif
#ifdef CVODE26
  jac = NewDenseMat(n,n);
#endif
  _dealloc = true;
}

balDynamicalSystem::~balDynamicalSystem() {
  if(_dealloc)
#ifdef CVODE25
    destroyMat(jac);
#endif
#ifdef CVODE26
  DestroyMat(jac);
#endif
}

const char * balDynamicalSystem::GetClassName () const {
  return "balDynamicalSystem";
}

balDynamicalSystem * balDynamicalSystem::Create () { 
  return new balDynamicalSystem;
}

balDynamicalSystem * balDynamicalSystem::Copy () {
  return new balDynamicalSystem(*this);
}

void balDynamicalSystem::Destroy () {
  delete this;
}

int balDynamicalSystem::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) {
  return ! CV_SUCCESS;
}

int balDynamicalSystem::RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void * sys) {
  balDynamicalSystem * bds = (balDynamicalSystem *) sys;

  // the first n components of the vector field
  int flag = bds->RHS(t,x,xdot,(void *)bds->GetParameters());

  if(! bds->IsExtended() || flag != CV_SUCCESS) {
    return flag;
  }
  
  
  // the Jacobian matrix
  if(bds->HasJacobian()) {
#ifdef CVODE25
    balDynamicalSystem::JacobianWrapper(bds->n,bds->jac,t,x,NULL,(void *)bds,NULL,NULL,NULL);
#endif
#ifdef CVODE26
    balDynamicalSystem::JacobianWrapper(bds->n,t,x,NULL,bds->jac,(void *)bds,NULL,NULL,NULL);
#endif
  }
  else {
    balDynamicalSystem::JacobianFiniteDifferences(bds->n,t,x,bds->jac,(void *)bds);
  }

  realtype Y[bds->n][bds->n], F[bds->n][bds->n];
  int i, j, k;
  
  // extend
  for(i=0, k=0; i<bds->n; i++) {
    for(j=0; j<bds->n; j++, k++)
      Y[j][i] = Ith(x,3+k);
  }
  // multiply the matrices
  for(i=0; i<bds->n; i++) {
    for(j=0; j<bds->n; j++) {
      F[i][j] = 0.;
      for(k=0; k<bds->n; k++)
				F[i][j] += IJth(bds->jac,i,k)*Y[k][j];
    }
  }
  // the last n*n components of the vector field
  for(i=0, k=0; i<bds->n; i++) {
    for(j=0; j<bds->n; j++,k++)
      Ith(xdot,3+k) = F[j][i];
  }
  
  return CV_SUCCESS;
}

#ifdef CVODE25
int balDynamicalSystem::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, 
				  N_Vector fy, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int balDynamicalSystem::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, 
				  DlsMat J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  return ! CV_SUCCESS;
}

#ifdef CVODE25
int balDynamicalSystem::JacobianWrapper (long int N, DenseMat J, realtype t, N_Vector x, 
					 N_Vector fy, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int balDynamicalSystem::JacobianWrapper (int N, realtype t, N_Vector x, N_Vector fy, 
					 DlsMat J, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  balDynamicalSystem * bds = (balDynamicalSystem *) sys;
#ifdef CVODE25
  return bds->Jacobian(N,J,t,x,fy,(void *)bds->GetParameters(),tmp1,tmp2,tmp3);
#endif
#ifdef CVODE26
  return bds->Jacobian(N,t,x,fy,J,(void *)bds->GetParameters(),tmp1,tmp2,tmp3);
#endif
}

#ifdef CVODE25
int balDynamicalSystem::JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
						   DenseMat J, void *sys) {
#endif
#ifdef CVODE26
int balDynamicalSystem::JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
						   DlsMat J, void *sys) {
#endif
  balDynamicalSystem * bds = (balDynamicalSystem *) sys;
  N_Vector ref, perturb;
  double eps = 1e-8;
  int i, j;

  ref = N_VNew_Serial(N);
  perturb = N_VNew_Serial(N);

  bds->RHS(t,x,ref,(void *)bds->GetParameters());

  for(j=0; j<N; j++) {
    Ith(x,j) = Ith(x,j) + eps;
    bds->RHS(t,x,perturb,(void *)bds->GetParameters());
    for(i=0; i<N; i++)
      IJth(J,i,j) = (Ith(perturb,i)-Ith(ref,i)) / eps;
    Ith(x,j) = Ith(x,j) - eps;
  }

  N_VDestroy_Serial(perturb);
  N_VDestroy_Serial(ref);

  return CV_SUCCESS;
}
  

int balDynamicalSystem::Events (realtype t, N_Vector x, realtype * event, void * data) {
  return ! CV_SUCCESS;
}

int balDynamicalSystem::EventsWrapper (realtype t, N_Vector x, realtype * event, void * sys) {
  balDynamicalSystem * bds = (balDynamicalSystem *) sys;
  return bds->Events(t,x,event,(void *)bds->GetParameters());
}

void balDynamicalSystem::EventsConstraints (realtype t, N_Vector x, int * constraints, void * data) {
}

void balDynamicalSystem::SetParameters (balParameters * bp) throw (balException) {
  if(bp->GetNumber() != GetNumberOfParameters())
    throw balException("Wrong number of parameters in balDynamicalSystem::SetParameters");
  pars = bp;
}
 
balParameters * balDynamicalSystem::GetParameters () const {
  return pars;
}

void balDynamicalSystem::SetDimension(int n_) {
  if(n_ <= 0)
    return;
  n = n_;
  nExt = n*(n+1);
  if(_dealloc)
#ifdef CVODE25
    destroyMat(jac);
  jac = newDenseMat(n,n);
#endif
#ifdef CVODE26
    DestroyMat(jac);
  jac = NewDenseMat(n,n);
#endif
  _dealloc = true;
}

int balDynamicalSystem::GetDimension() const {
  return (ext ? nExt : n);
}

int balDynamicalSystem::GetOriginalDimension() const {
  return n;
}

void balDynamicalSystem::Extend(bool extend) {
  ext = extend;
}

bool balDynamicalSystem::IsExtended() const {
  return ext;
}

bool balDynamicalSystem::SpecialOptions(const void *opt) {
  printf("balDynamicalSystem::SpecialOptions\n");
  return false;
}

bool balDynamicalSystem::HasJacobian() const {
  return false;
}
  
bool balDynamicalSystem::HasEvents() const {
  return false;
}

bool balDynamicalSystem::HasEventsConstraints() const {
  return false;
}
  
void balDynamicalSystem::Reset() {
}

void balDynamicalSystem::ManageEvents(realtype t, N_Vector X, int * events, int * constraints) {
}
  
int balDynamicalSystem::GetNumberOfEvents() const {
  return nev;
}

int balDynamicalSystem::GetNumberOfParameters() const {
  return p;
}
  
void balDynamicalSystem::SetNumberOfParameters(int p_) {
  if(p_ >= 0)
    p = p_;
}

void balDynamicalSystem::SetNumberOfEvents(int nev_) {
  if(nev_ >= 0)
    nev = nev_;
}
