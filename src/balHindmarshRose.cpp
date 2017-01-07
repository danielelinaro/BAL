/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balHindmarshRose.cpp
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
 * \file balHindmarshRose.cpp
 * \brief Implementation of the class HindmarshRose
 */

#include "balHindmarshRose.h"
#include <iostream>

bal::DynamicalSystem* HindmarshRoseFactory() {
  return new bal::HindmarshRose();
}

namespace bal {

const double HindmarshRose::xrest = -1.6;

HindmarshRose::HindmarshRose() : DynamicalSystem(3,4,3,false) {
#ifdef DEBUG
  std::cout << "HindmarshRose constructor.\n";
#endif
  xderiv = N_VNew_Serial(GetDimension());
}

HindmarshRose::HindmarshRose(const HindmarshRose& hr) : DynamicalSystem(hr) {
#ifdef DEBUG
  std::cout << "HindmarshRose copy constructor.\n";
#endif
  xderiv = N_VNew_Serial(hr.GetDimension());
  for(int i = 0; i < hr.GetDimension(); i++)
    Ith(xderiv,i) = Ith(hr.xderiv,i);
}

HindmarshRose::~HindmarshRose() {
#ifdef DEBUG
  std::cout << "HindmarshRose destructor.\n";
#endif
  N_VDestroy_Serial(xderiv);
}

int HindmarshRose::RHS (realtype t, N_Vector x, N_Vector xdot, void *sys) {
  realtype x1, x2, x3;
  realtype b, I, u, s;
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  const Parameters *parameters = ds->GetParameters();
  
  b = parameters->At(0);
  I = parameters->At(1);
  u = parameters->At(2);
  s = parameters->At(3);

  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);

  Ith (xdot, 0) = x2 - x1*x1*x1 + b*x1*x1 + I - x3;
  Ith (xdot, 1) = 1 - 5*x1*x1 - x2;
  Ith (xdot, 2) = u*(s*(x1 - xrest) - x3);

  return CV_SUCCESS;
}

#ifdef CVODE25
int HindmarshRose::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
				void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int HindmarshRose::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
				void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  realtype b, I, u, s;
  realtype x1, x2, x3;
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  const Parameters *parameters = ds->GetParameters();
  
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
 
  b = parameters->At(0);
  I = parameters->At(1);
  u = parameters->At(2);
  s = parameters->At(3);

  IJth (J, 0, 0) = -3.0*x1*x1 + 2.0*b*x1;
  IJth (J, 0, 1) =  1.0;
  IJth (J, 0, 2) = -1.0;
  IJth (J, 1, 0) = -10.0 * x1;
  IJth (J, 1, 1) = -1.0;
  IJth (J, 1, 2) =  0.0;
  IJth (J, 2, 0) =  u*s;
  IJth (J, 2, 1) =  0.0;
  IJth (J, 2, 2) = -u;

  return CV_SUCCESS;
}

int HindmarshRose::Events (realtype t, N_Vector x, realtype *event, void *sys) {
  RHS(t, x, xderiv, sys);
  for(int i=0; i<GetNumberOfEvents(); i++)
    event[i] = Ith(xderiv,i);
  return CV_SUCCESS;
}

void HindmarshRose::EventsConstraints (realtype t, N_Vector x, int *constraints, void *sys) {
  realtype b, u, s;
  realtype x1, x2, x3;
  realtype ris[3], xdot[3];
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  const Parameters *parameters = ds->GetParameters();
  
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
 
  b = parameters->At(0);
  u = parameters->At(2);
  s = parameters->At(3);

  RHS(t, x, xderiv, sys);
  for(int i=0; i<GetDimension(); i++)
    xdot[i] = Ith(xderiv,i);
  
  ris[0] = xdot[1] - xdot[2] - 3*x1*x1*xdot[0] + 2*b*x1*xdot[0];
  ris[1] = -10*x1*xdot[0] - xdot[1];
  ris[2] = u*s*xdot[0] - u*xdot[2];
	
  for(int i=0; i<GetNumberOfEvents(); i++)
    constraints[i] = (ris[i] < 0 ? 1 : 0);
}

bool HindmarshRose::HasJacobian() const {
  return true;
}
 
bool HindmarshRose::HasEvents() const {
  return true;
}
 
bool HindmarshRose::HasEventsConstraints() const {
  return true;
}
 
DynamicalSystem* HindmarshRose::Clone() const {
  return new HindmarshRose(*this);
}

}

