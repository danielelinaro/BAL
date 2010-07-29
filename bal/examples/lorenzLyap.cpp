/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    solver.cpp
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
#include "balObject.h"
#include "balLorenz.h"
#include "balParameters.h"
#include "balODESolver.h"
#include "balSolution.h"
using namespace std;

// TEST balODESolver SolveLyapunovExponents
int main(int argc, char *argv[]) {
	
  // parameters
  balParameters * pars = balParameters::Create();
  pars->SetNumber(3);
  pars->At(0) = 16.0;
  pars->At(1) = 45.92;
  pars->At(2) = 4;
  
  // Lorenz
  balLorenz *lor = balLorenz::Create();
  lor->SetParameters(pars);
  lor->Extend(true);
  
  realtype x0[] = {0,1,0};
  balODESolver * solver = balODESolver::Create();
  solver->SetDynamicalSystem(lor);
  solver->SetTransientDuration(1000);
  solver->SetFinalTime(2000);
  solver->SetTimeStep(0.01);
  solver->SetLyapunovTimeStep(2);
  solver->SetIntegrationMode(balLYAP);
  solver->SetX0(x0);
  solver->Solve();
  
  realtype * lp = solver->GetLyapunovExponents();
  printf("\nLyapunov Exponents: %e,%e,%e \n",lp[0],lp[1],lp[2]);
	
	
// balSolution *sol = solver->GetSolution();
//  int r,c,i,j;
//  double *buffer;
//  sol->GetSize(&r,&c);
//  buffer = sol->GetData();
//  for(i=0; i<r; i++) {
//    for(j=0; j<c-1; j++)
//      printf("%f ", buffer[i*c+j]);
//    printf("%d\n", (int) buffer[i*c+j]);
//  }
  
  //sol->Destroy();
  solver->Destroy();
  lor->Destroy();
  pars->Destroy();
  
  return 0;
}

