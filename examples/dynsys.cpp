/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    dynsys.cpp
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

#include <iostream>
#include <nvector/nvector_serial.h>
#include <cstdlib>
#include "balObject.h"
#include "balCommon.h"
#include "balDynamicalSystem.h"
#include "balHindmarshRose.h"
#include "balParameters.h"
using namespace bal;

#ifdef CVODE25
void PrintJacobian(int n, DenseMat J) {
#endif
#ifdef CVODE26
void PrintJacobian(int n, DlsMat J) {
#endif
  int i, j;
  for(i=0; i<n; i++) {
    printf("| ");
    for(j=0; j<n; j++)
      printf("%10.4f ", DENSE_ELEM(J,i,j));
    printf("|\n");
  }
}

// TEST balDynamicalSystem and balHindmarshRose
int main(int argc, char *argv[]) {
	 
  // HindmarshRose
  int i, n;
  HindmarshRose hr;
  n = hr.GetDimension();
  N_Vector x = N_VNew_Serial(n);
  N_Vector xdot = N_VNew_Serial(n);
  
  // parameters
  Parameters pars(hr.GetNumberOfParameters());
  hr.SetParameters(&pars);
  pars[0] = 3.0;
  pars[1] = 5.0;
  pars[2] = 0.01;
  pars[3] = 4.0;
 
  if(argc == n+1) {
    for(i=0; i<n; i++)
      NV_Ith_S(x,i) = atof(argv[i+1]);
  }
  else {
    for(i=0; i<n; i++)
      NV_Ith_S(x,i) = 0.0;
  }
  DynamicalSystem::RHSWrapper(0, x, xdot, (void *) &hr);

  std::cout << "xdot = (";
  for(i=0; i<n-1; i++) {
    std::cout << NV_Ith_S(xdot,i) << ",";
  }
  std::cout << NV_Ith_S(xdot,i) << ")" << std::endl;
  
#ifdef CVODE25
  DenseMat jac = newDenseMat(n,n);
  DynamicalSystem::JacobianWrapper(n, jac, 0, x, NULL, (void *) &hr, NULL, NULL, NULL);
#endif
#ifdef CVODE26
  DlsMat jac = NewDenseMat(n,n);
  DynamicalSystem::JacobianWrapper(n, 0, x, NULL, jac, (void *) &hr, NULL, NULL, NULL);
#endif
  printf(">> Exact Jacobian matrix <<\n");
  PrintJacobian(n,jac);
  DynamicalSystem::JacobianFiniteDifferences(n, 0, x, jac, (void *) &hr);
  printf(">> Approximated Jacobian matrix <<\n");
  PrintJacobian(n,jac);

#ifdef CVODE25
  destroyMat(jac);
#endif
#ifdef CVODE26
  DestroyMat(jac);
#endif

  N_VDestroy_Serial(x);
  N_VDestroy_Serial(xdot);

  HindmarshRose hrcopy(hr);
  std::cout << "hrcopy.n = " << hrcopy.GetDimension() << std::endl;
  std::cout << "hrcopy.xrest = " << hrcopy.xrest << std::endl;

  return 0;
}

