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

namespace bal {

DynamicalSystem::DynamicalSystem() {
  n = 0;
  p = 0;
  nev = 0;
  ext = false;
  nExt = 0;
  jac = NULL;
  dealloc_ = false;
}

DynamicalSystem::DynamicalSystem(const DynamicalSystem& system) : pars(system.pars) {
  n = system.n;
  nev = system.nev;
  p = pars.GetNumber();
  nExt = system.nExt;
  ext = system.ext;
#ifdef CVODE25
  jac = newDenseMat(n,n);
#endif
#ifdef CVODE26
  jac = NewDenseMat(n,n);
#endif
  dealloc_ = true;
}

DynamicalSystem::~DynamicalSystem() {
  if(dealloc_)
#ifdef CVODE25
    destroyMat(jac);
#endif
#ifdef CVODE26
    DestroyMat(jac);
#endif
}

int DynamicalSystem::RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void * sys) {
  DynamicalSystem * bds = (DynamicalSystem *) sys;

  // the first n components of the vector field
  int flag = bds->RHS(t,x,xdot,(void *)bds->GetParameters());

  if(! bds->IsExtended() || flag != CV_SUCCESS) {
    return flag;
  }
  
  // the Jacobian matrix
  if(bds->HasJacobian()) {
#ifdef CVODE25
    DynamicalSystem::JacobianWrapper(bds->n,bds->jac,t,x,NULL,(void *)bds,NULL,NULL,NULL);
#endif
#ifdef CVODE26
    DynamicalSystem::JacobianWrapper(bds->n,t,x,NULL,bds->jac,(void *)bds,NULL,NULL,NULL);
#endif
  }
  else {
    DynamicalSystem::JacobianFiniteDifferences(bds->n,t,x,bds->jac,(void *)bds);
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
int DynamicalSystem::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, 
				  N_Vector fy, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int DynamicalSystem::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, 
				  DlsMat J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
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
  DynamicalSystem * bds = (DynamicalSystem *) sys;
#ifdef CVODE25
  return bds->Jacobian(N,J,t,x,fy,(void *)bds->GetParameters(),tmp1,tmp2,tmp3);
#endif
#ifdef CVODE26
  return bds->Jacobian(N,t,x,fy,J,(void *)bds->GetParameters(),tmp1,tmp2,tmp3);
#endif
}

#ifdef CVODE25
int DynamicalSystem::JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
						   DenseMat J, void *sys) {
#endif
#ifdef CVODE26
int DynamicalSystem::JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
						   DlsMat J, void *sys) {
#endif
  DynamicalSystem * bds = (DynamicalSystem *) sys;
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
  

int DynamicalSystem::Events (realtype t, N_Vector x, realtype * event, void * data) {
  return ! CV_SUCCESS;
}

int DynamicalSystem::EventsWrapper (realtype t, N_Vector x, realtype * event, void * sys) {
  DynamicalSystem * bds = (DynamicalSystem *) sys;
  return bds->Events(t,x,event,(void *)bds->GetParameters());
}

void DynamicalSystem::EventsConstraints (realtype t, N_Vector x, int * constraints, void * data) {
}

void DynamicalSystem::SetParameters (const Parameters& bp) throw (Exception) {
  pars = bp;
}
 
Parameters* DynamicalSystem::GetParameters () const {
  return &pars;
}

void DynamicalSystem::SetDimension(int n_) {
  if(n_ <= 0)
    return;
  n = n_;
  nExt = n*(n+1);
  if(dealloc_)
#ifdef CVODE25
    destroyMat(jac);
  jac = newDenseMat(n,n);
#endif
#ifdef CVODE26
    DestroyMat(jac);
  jac = NewDenseMat(n,n);
#endif
  dealloc_ = true;
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

void DynamicalSystem::ManageEvents(realtype t, N_Vector X, int * events, int * constraints) {
}
  
int DynamicalSystem::GetNumberOfEvents() const {
  return nev;
}

int DynamicalSystem::GetNumberOfParameters() const {
  return p;
}
  
void DynamicalSystem::SetNumberOfParameters(int p_) {
  if(p_ >= 0)
    p = p_;
}

void DynamicalSystem::SetNumberOfEvents(int nev_) {
  if(nev_ >= 0)
    nev = nev_;
}

} // namespace bal

