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
  Parameters pars(4);
  pars[0] = 2.5;
  pars[1] = 3.5;
  pars[2] = 0.01;
  pars[3] = 4.0;
  
  // HindmarshRose
  HindmarshRose hr;
  hr.SetParameters(&pars);
  
  ODESolver solver;
  solver.SetDynamicalSystem(&hr);
  solver.SetTransientDuration(1000.0);
  solver.SetFinalTime(2000.0);
  solver.SetMaxNumberOfIntersections(1000);
  solver.SetTimeStep(0.05);
  solver.HaltAtEquilibrium(true);
  solver.HaltAtCycle(false);
  solver.SetIntegrationMode(TRAJ);
  printf("Computing the whole trajectory... ");
  solver.Solve();
  printf("done.\n");
  solver.SetIntegrationMode(EVENTS);
  printf("Computing only the events... ");
  solver.Solve();
  printf("done.\n");
  solver.SetIntegrationMode(BOTH);
  printf("Computing trajectory and events... ");
  solver.Solve();
  solver.SaveOrbit("bursting.dat");
  printf("done.\n");

  return 0;
}

