/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balDynamicalSystem.cpp
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
  n = system.GetDimension();
  nev = system.GetNumberOfEvents();
  pars = system.GetParameters(); // non rialloco spazio parametri
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

int balDynamicalSystem::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) {
  return ! CV_SUCCESS;
}

int balDynamicalSystem::RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void * sys) {
  balDynamicalSystem * bds = (balDynamicalSystem *) sys;
  if(! bds->IsExtended()) {
    return bds->RHS(t,x,xdot,(void *)bds->GetParameters());
  }
  
  // the first n components of the vector field
  bds->RHS(t,x,xdot,(void *)bds->GetParameters());
  
  // the Jacobian matrix
#ifdef CVODE25
  balDynamicalSystem::JacobianWrapper(bds->n,bds->jac,t,x,NULL,(void *)bds,NULL,NULL,NULL);
#endif
#ifdef CVODE26
  balDynamicalSystem::JacobianWrapper(bds->n,t,x,NULL,bds->jac,(void *)bds,NULL,NULL,NULL);
#endif

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

int balDynamicalSystem::Events (realtype t, N_Vector x, realtype * event, void * data) {
  return ! CV_SUCCESS;
}

int balDynamicalSystem::EventsWrapper (realtype t, N_Vector x, realtype * event, void * sys) {
  balDynamicalSystem * bds = (balDynamicalSystem *) sys;
  return bds->Events(t,x,event,(void *)bds->GetParameters());
}

void balDynamicalSystem::EventsConstraints (realtype t, N_Vector x, int * constraints, void * data) {}

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

void balDynamicalSystem::Extend(bool extend) {
  ext = extend;
}

bool balDynamicalSystem::IsExtended() const {
  return ext;
}

bool balDynamicalSystem::SpecialOptions(void *opt) {
  return false;
}
