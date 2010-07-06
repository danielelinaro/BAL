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

#include <cstdio>
#include <iostream>
#include <nvector/nvector_serial.h>
#include "balObject.h"
#include "balCommon.h"
#include "balDynamicalSystem.h"
#include "balHindmarshRose.h"
#include "balParameters.h"
using namespace std;

#ifdef CVODE25
void PrintJacobian(int n, DenseMat J);
void PrintJacobian(int n, DenseMat J) {
#endif
#ifdef CVODE26
void PrintJacobian(int n, DlsMat J);
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
	
  // parameters
  balParameters * pars = balParameters::Create();
  pars->SetNumber(4);
  pars->At(0) = 3.0;
  pars->At(1) = 5.0;
  pars->At(2) = 0.01;
  pars->At(3) = 4.0;
  
  // balDynamicalSystem
  balDynamicalSystem * dynsys = balDynamicalSystem::Create();
  cout << dynsys->GetClassName() << endl;
  dynsys->Destroy();
  
  // balHindmarshRose
  int i, n;
  balDynamicalSystem * hr = balHindmarshRose::Create();
  n = hr->GetDimension();
  N_Vector x = N_VNew_Serial(n);
  N_Vector xdot = N_VNew_Serial(n);

  if(argc == n+1) {
    for(i=0; i<n; i++)
      NV_Ith_S(x,i) = atof(argv[i+1]);
  }
  else {
    for(i=0; i<n; i++)
      NV_Ith_S(x,i) = 0.0;
  }
  cout << hr->GetClassName() << endl;
  hr->SetParameters(pars);
  balDynamicalSystem::RHSWrapper(0,x,xdot,hr);

  cout << "xdot = (";
  for(i=0; i<n-1; i++) {
    cout << NV_Ith_S(xdot,i) << ",";
  }
  cout << NV_Ith_S(xdot,i) << ")" << endl;
  
#ifdef CVODE25
  DenseMat jac = newDenseMat(n,n);
  balDynamicalSystem::JacobianWrapper(n,jac,0,x,NULL,hr,NULL,NULL,NULL);
#endif
#ifdef CVODE26
  DlsMat jac = NewDenseMat(n,n);
  balDynamicalSystem::JacobianWrapper(n,0,x,NULL,jac,hr,NULL,NULL,NULL);
#endif
  printf("\n>> Exact Jacobian matrix <<\n");
  PrintJacobian(n,jac);
  balDynamicalSystem::JacobianFiniteDifferences(n,0,x,jac,hr);
  printf("\n>> Approximated Jacobian matrix <<\n");
  PrintJacobian(n,jac);

#ifdef CVODE25
  destroyMat(jac);
#endif
#ifdef CVODE26
  DestroyMat(jac);
#endif

  N_VDestroy_Serial(x);
  N_VDestroy_Serial(xdot);
  hr->Destroy();
  
  return 0;
}

