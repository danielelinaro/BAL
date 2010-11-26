/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balRossler.cpp
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
 * \file balRossler.cpp
 * \brief Implementation of the class balRossler
 */

#include "balRossler.h"

balDynamicalSystem* balRosslerFactory() {
  return balRossler::Create();
}

balRossler::balRossler() {
  SetDimension(3);
  SetNumberOfParameters(3);
  SetNumberOfEvents(GetDimension());
  xderiv = N_VNew_Serial(GetDimension());
}

balRossler::balRossler(const balRossler& ros) : balDynamicalSystem(ros) {
  xderiv = N_VNew_Serial(ros.GetDimension());
  for(int i = 0; i < ros.GetDimension(); i++)
    Ith(xderiv,i)=Ith(ros.xderiv,i);
}

balRossler::~balRossler() {
  N_VDestroy_Serial(xderiv);
}

balRossler * balRossler::Create () {
  return new balRossler;
}

balRossler * balRossler::Copy (balRossler *ros) {
  return new balRossler(*ros);
}

balDynamicalSystem * balRossler::Clone() const {
  return new balRossler(*this);
}

void balRossler::Destroy () {
  delete this;
}

const char * balRossler::GetClassName () const {
  return "balRossler";
}

int balRossler::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) {
  realtype x1, x2, x3;
  realtype a, b, c;
  balParameters * parameters;
  
  parameters = (balParameters *) data;
  a = parameters->At(0);
  b = parameters->At(1);
  c = parameters->At(2);
  
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
  
  Ith (xdot, 0) = - x2 - x3;
  Ith (xdot, 1) = x1 + a*x2;
  Ith (xdot, 2) = b + x3*(x1 - c);
  
  return CV_SUCCESS;
}

#ifdef CVODE25
int balRossler::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
			  void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int balRossler::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
			  void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  realtype a, b, c;
  realtype x1, x2, x3;
  balParameters * parameters;
  
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
  
  parameters = (balParameters *) jac_data;
  a = parameters->At(0);
  b = parameters->At(1);
  c = parameters->At(2);
  
  IJth (J, 0, 0) =  0.0;
  IJth (J, 0, 1) = -1.0;
  IJth (J, 0, 2) = -1.0;
  IJth (J, 1, 0) =  1.0;
  IJth (J, 1, 1) =  a;
  IJth (J, 1, 2) =  0.0;
  IJth (J, 2, 0) =  x3;
  IJth (J, 2, 1) =  0.0;
  IJth (J, 2, 2) =  x1-c;
  
  return CV_SUCCESS;
}
 
 int balRossler::Events (realtype t, N_Vector x, realtype * event, void * data) {
   RHS(t,x,xderiv,data);
   for(int i=0; i<GetNumberOfEvents(); i++)
     event[i] = Ith(xderiv,i);
   return CV_SUCCESS;
 }
 
void balRossler::EventsConstraints (realtype t, N_Vector x, int * constraints, void * data) {
  realtype a, b, c;
  realtype x1, x2, x3;
  realtype ris[3], xdot[3];
  balParameters * parameters;
  
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
 
  parameters = (balParameters *) data;
  a = parameters->At(0);
  b = parameters->At(1);
  c = parameters->At(2);

  RHS(t,x,xderiv,data);
  for(int i=0; i<GetDimension(); i++)
    xdot[i] = Ith(xderiv,i);
  
  ris[0] = - xdot[1] - xdot[2];
  ris[1] = xdot[0] + a*xdot[1];
  ris[2] = xdot[0]*x3 + xdot[2]*(x1-c);
	
  for(int i=0; i<GetNumberOfEvents(); i++)
    constraints[i] = (ris[i] < 0 ? 1 : 0);
}

 bool balRossler::HasJacobian() const {
   return (IsExtended() ? false : true);
 }
 
 bool balRossler::HasEvents() const {
   return true;
 }
 
 bool balRossler::HasEventsConstraints() const {
   return true;
 }
 
