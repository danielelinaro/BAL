/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balLorenz.cpp
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
 * \file balLorenz.cpp
 * \brief Implementation of the class Lorenz
 */

#include "balLorenz.h"
#include <iostream>

bal::DynamicalSystem* LorenzFactory() {
  return new bal::Lorenz;
}

namespace bal {

Lorenz::Lorenz() : DynamicalSystem(3,3,0,false) {
}

Lorenz::Lorenz(const Lorenz& lor) : DynamicalSystem(lor) {
#ifdef DEBUG
  std::cout << "Lorenz copy constructor.\n";
#endif
}

Lorenz::~Lorenz() {}

bool Lorenz::HasJacobian() const {
  return (IsExtended() ? false : true);
}

int Lorenz::RHS (realtype t, N_Vector x, N_Vector xdot, void *sys) {
  // the state of the system
  realtype x1, x2, x3;
  // the parameters
  realtype sigma, rho, beta;
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  Parameters *parameters = ds->GetParameters();
  
  // parameters
  sigma = parameters->At(0);
  rho = parameters->At(1);
  beta = parameters->At(2);
  // state
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
  
  Ith (xdot, 0) = sigma*(x2-x1);
  Ith (xdot, 1) = -x1*x3 + rho*x1 - x2;
  Ith (xdot, 2) = x1*x2 - beta*x3;
  
  return CV_SUCCESS;
}

#ifdef CVODE25
int Lorenz::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
			 void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int Lorenz::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
			 void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  realtype x1, x2, x3;
  realtype sigma, rho, beta;
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  Parameters *parameters = ds->GetParameters();
  
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
 
  sigma = parameters->At(0);
  rho = parameters->At(1);
  beta = parameters->At(2);

  IJth (J, 0, 0) = -sigma;
  IJth (J, 0, 1) =  sigma;
  IJth (J, 0, 2) =  0.0;
  IJth (J, 1, 0) =  rho-x3;
  IJth (J, 1, 1) = -1.0;
  IJth (J, 1, 2) = -x1;
  IJth (J, 2, 0) =  x2;
  IJth (J, 2, 1) =  x1;
  IJth (J, 2, 2) = -beta;

  return CV_SUCCESS;
}

DynamicalSystem* Lorenz::Clone() const {
  return new Lorenz(*this);
}

} // namespace bal

