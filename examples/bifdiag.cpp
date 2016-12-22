/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    bifdiag.cpp
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
#include <nvector/nvector_serial.h>
#include "balObject.h"
#include "balDynamicalSystem.h"
#include "balHindmarshRose.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"
#include "balBifurcationParameters.h"
using namespace bal;

// TEST BifurcationDiagram
int main(int argc, char *argv[]) {
  int steps[4] = {10,1,1,1};
  realtype x0[3] = {0.5,0.5,0.5};
  HindmarshRose hr;
  BifurcationDiagram bifd;
  BifurcationParameters *bp = hr.GetParameters();

  bp->SetIthParameterLowerBound(0,2.9);
  bp->SetIthParameterUpperBound(0,3.18);
  bp->SetIthParameter(1,3);
  bp->SetIthParameter(2,0.01);
  bp->SetIthParameter(3,4.0);
  bp->SetNumberOfSteps(steps);

  bifd.SetDynamicalSystem(&hr);
  bifd.RestartFromX0(true);
  bifd.GetODESolver()->SetIntegrationMode(BOTH);
  bifd.GetODESolver()->HaltAtEquilibrium(true);
  bifd.GetODESolver()->HaltAtCycle(true);
  bifd.GetODESolver()->SetInitialTime(10.);
  bifd.GetODESolver()->SetTransientDuration(10);
  bifd.GetODESolver()->SetFinalTime(1e4);
  bifd.GetODESolver()->SetMaxNumberOfIntersections(500);
  bifd.GetODESolver()->SetX0(x0);
  bifd.SetFilename("hr_comp.h5",false);
  bifd.SetNumberOfThreads(argc > 1 ? atoi(argv[1]) : 2);
  bifd.SetNumberOfThreads(1);

  bifd.ComputeDiagram();
  bifd.SaveSummaryData("hr_comp.classified");

  return 0;
}

