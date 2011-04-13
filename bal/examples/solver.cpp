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
#include "balHindmarshRose.h"
#include "balParameters.h"
#include "balODESolver.h"
using namespace bal;

// TEST ODESolver
int main(int argc, char *argv[]) {
	
  // parameters
  Parameters * pars = Parameters::Create();
  pars->SetNumber(4);
  pars->At(0) = 3.0;
  pars->At(1) = 5.0;
  pars->At(2) = 0.01;
  pars->At(3) = 4.0;
  
  // HindmarshRose
  HindmarshRose *hr = HindmarshRose::Create();
  hr->SetParameters(pars);
  
  ODESolver * solver = ODESolver::Create();
  solver->SetDynamicalSystem(hr);
  solver->SetTransientDuration(0.0);
  solver->SetFinalTime(1000.0);
  solver->HaltAtEquilibrium(true);
  solver->SetIntegrationMode(balTRAJ);
  printf("Computing the whole trajectory... ");
  solver->Solve();
  printf("done.\n");
  solver->SetIntegrationMode(balEVENTS);
  printf("Computing only the events... ");
  solver->Solve();
  printf("done.\n");
  solver->SetIntegrationMode(balBOTH);
  printf("Computing trajectory and events... ");
  solver->Solve();
  printf("done.\n");
  pars->At(1) = 1.0;
  solver->SetIntegrationMode(balTRAJ);
  printf("Computing the whole trajectory and stopping at an equilibrium point... ");
  solver->Solve();
  printf("done.\n");
  solver->Destroy();
  
  return 0;
}

