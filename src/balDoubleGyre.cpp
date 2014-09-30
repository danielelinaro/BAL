/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balDoubleGyre.cpp
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
 * \file balDoubleGyre.cpp
 * \brief Implementation of the class DoubleGyre
 */

#include <cmath>
#include "balDoubleGyre.h"

bal::DynamicalSystem* DoubleGyreFactory() {
  return new bal::DoubleGyre;
}

namespace bal {

DoubleGyre::DoubleGyre() {
  SetDimension(2);
  SetNumberOfParameters(3);
  pi = 4*atan(1.0);
}

DoubleGyre::DoubleGyre(const DoubleGyre& gyre) : DynamicalSystem(gyre) {
  pi = gyre.pi;
}

DoubleGyre::~DoubleGyre() {
}

realtype DoubleGyre::f(realtype t, realtype x, realtype omega, realtype eps) const {
  realtype s;
  s = sin(omega*t);
  return eps*s*x*x + (1-2*eps*s)*x;
}

realtype DoubleGyre::df(realtype t, realtype x, realtype omega, realtype eps) const {
  realtype s;
  s = sin(omega*t);
  return 2*eps*s*x + (1-2*eps*s);
}

int DoubleGyre::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) {
  realtype x1, x2;
  realtype A, omega, eps;
  realtype fx, dfx;
  Parameters * parameters;
  
  parameters = (Parameters *) data;
  A = parameters->At(0);
  omega = parameters->At(1);
  eps = parameters->At(2);

  x1 = Ith (x, 0);
  x2 = Ith (x, 1);

  fx = f(t,x1,omega,eps);
  dfx = df(t,x1,omega,eps);

  Ith (xdot, 0) = -pi*A*sin(pi*fx)*cos(pi*x2);
  Ith (xdot, 1) =  pi*A*cos(pi*fx)*sin(pi*x2)*dfx;

  return CV_SUCCESS;
}

#ifdef CVODE25
int DoubleGyre::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
				void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int DoubleGyre::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
				void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  return CV_SUCCESS;
}

bool DoubleGyre::HasJacobian() const {
  return false;
}

DynamicalSystem* DoubleGyre::Clone() const {
  return new DoubleGyre(*this);
}

} // namespace bal

