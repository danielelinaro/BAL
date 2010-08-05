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

// TEST balBifurcationDiagram
int main(int argc, char *argv[]) {

  int steps[4] = {1,1,1,1};
  realtype x0[3] = {0.5,0.5,0.5};
  balBifurcationParameters * bp = balBifurcationParameters::Create();
  bp->SetNumber(4);
  bp->SetIthParameter(0,3.065);
  bp->SetIthParameter(1,3);
  bp->SetIthParameter(2,0.01);
  bp->SetIthParameter(3,4.0);
  bp->SetNumberOfSteps(steps);
  balHindmarshRose * hr = balHindmarshRose::Create();
  hr->SetParameters(bp);
  balBifurcationDiagram * bifd = balBifurcationDiagram::Create();
  bifd->SetDynamicalSystem(hr);
  bifd->GetODESolver()->SetIntegrationMode(balBOTH);
  bifd->GetODESolver()->HaltAtEquilibrium(true);
  bifd->GetODESolver()->HaltAtCycle(false);
  bifd->GetODESolver()->SetTransientDuration(2e3);
  bifd->GetODESolver()->SetFinalTime(1e4);
  bifd->GetODESolver()->SetMaxNumberOfIntersections(200);
  bifd->GetODESolver()->SetX0(x0);
  bifd->SetFilename("hr.h5");
  bifd->SetNumberOfThreads(argc > 1 ? atoi(argv[1]) : 2);

  bifd->ComputeDiagram();
  bifd->SaveClassificationData("hr.classified");

  bifd->Destroy();
  hr->Destroy();
  bp->Destroy();
  
  return 0;
}

