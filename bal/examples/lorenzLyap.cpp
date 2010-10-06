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

/**************** TEST balODESolver SolveLyapunovExponents **********************
 *										*
 *		lyapunov exponents calculated using base 2 logarithm		*
 *		in formula.							*
 *									        *
 *		test parameters [16.0 45.92 4.0] from Wolf's article		*
 *		'Determining Lyapunov Exponents From a Time Series'		*
 *		Lyapunov Spectrum: l1 = 2.16, l2 = 0.00, l3 = -32.4		*
 *		for t=10000							*
 *										*
 ********************************************************************************/


int main(int argc, char *argv[]) {
	
  // parameters
  balParameters * pars = balParameters::Create();
  pars->SetNumber(3);
  pars->At(0) = 16;
  pars->At(1) = 45.92;
  pars->At(2) = 4.0;
  
  // Lorenz
  balLorenz *lor = balLorenz::Create();
  lor->SetParameters(pars);
  
  // Setting ODESolver fields 
  realtype x0[] = {10,1,0};
  balODESolver * solver = balODESolver::Create();
  solver->SetDynamicalSystem(lor);
  solver->SetTransientDuration(1000);
  solver->SetFinalTime(1.1e4);
  solver->SetTimeStep(5);
  solver->SetLyapunovTimeStep(5);
  solver->SetIntegrationMode(balLYAP);
  solver->SetX0(x0);
	
  // Calculating Lyapunov Exponents
  solver->Solve();
  
  for(int i=0; i<lor->GetOriginalDimension(); i++){
    printf("%e ",solver->GetLyapunovExponents()[i]);
  }
  printf("\n");
  
  solver->Destroy();
  lor->Destroy();
  pars->Destroy();
  
  return 0;
}

