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
	pars = NULL;
}

balDynamicalSystem::balDynamicalSystem(const balDynamicalSystem& system) {
	n = system.GetDimension();
	nev = system.GetNumberOfEvents();
	pars = system.GetParameters(); // non rialloco spazio parametri
	p = (system.GetParameters())->GetNumber();
}

balDynamicalSystem::~balDynamicalSystem() {}

int balDynamicalSystem::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) {
	return ! CV_SUCCESS;
}

int balDynamicalSystem::RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void * sys) {
	balDynamicalSystem * bds = (balDynamicalSystem *) sys;
	return bds->RHS(t,x,xdot,(void *)bds->GetParameters());
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

