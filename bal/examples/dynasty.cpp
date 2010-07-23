/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    dynasty.cpp
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
#include "balDynasty.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"
#include "balBifurcationParameters.h"

// TEST balBifurcationDiagram
int main(int argc, char *argv[]) {

  int steps[7] = {1001,1,1,1,1,1,1};
  realtype x0[3] = {0.5,0.5,0.5};
  balBifurcationParameters * bp = balBifurcationParameters::Create();
  bp->SetNumber(7);
  //bp->SetIthParameter(0,-3); // r
  bp->SetIthParameterLowerBound(0,-4); // r
  bp->SetIthParameterUpperBound(0,-1); // r
  bp->SetIthParameter(1,2.2); // e
  bp->SetIthParameter(2,0.17); // b
  bp->SetIthParameter(3,0.42); // d
  bp->SetIthParameter(4,0.09); // g
  bp->SetIthParameter(5,0.1); // h
  bp->SetIthParameter(6,0.4); // q
  bp->SetNumberOfSteps(steps);
  balDynasty * dynasty = balDynasty::Create();
  dynasty->SetParameters(bp);
  dynasty->SpecialOptions((const void *) "maxima");
  balBifurcationDiagram * bifd = balBifurcationDiagram::Create();
  bifd->SetDynamicalSystem(dynasty);
  bifd->GetODESolver()->SetIntegrationMode(balEVENTS);
  bifd->GetODESolver()->HaltAtEquilibrium(true);
  bifd->GetODESolver()->HaltAtCycle(false);
  bifd->GetODESolver()->SetTransientDuration(5e3);
  bifd->GetODESolver()->SetFinalTime(1e6);
  bifd->GetODESolver()->SetEquilibriumTolerance(1e-6);
  bifd->GetODESolver()->SetRelativeTolerance(1e-8);
  bifd->GetODESolver()->SetAbsoluteTolerance(1e-13);
  bifd->GetODESolver()->SetMaxNumberOfIntersections(500);
  bifd->GetODESolver()->SetX0(x0);
  bifd->SetFilename("dynasty.h5");
  bifd->SetNumberOfThreads(argc > 1 ? atoi(argv[1]) : 2);
  bifd->RestartFromX0(true);
  bifd->ComputeDiagram();
  bifd->SaveClassificationData("dynasty.classified");

  bifd->Destroy();
  dynasty->Destroy();
  bp->Destroy();
  
  return 0;
}

