/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    basin.cpp
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
#include <cstdlib>
#include <nvector/nvector_serial.h>
#include "balObject.h"
#include "balDynamicalSystem.h"
#include "balHindmarshRose.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"
#include "balParameters.h"

// TEST balBifurcationDiagram
int main(int argc, char *argv[]) {

  balParameters * par = balParameters::Create();
  par->SetNumber(4);
  par->At(0) = 2.96;
  par->At(1) = 3;
  par->At(2) = 0.01;
  par->At(3) = 4.0;

  balHindmarshRose * hr = balHindmarshRose::Create();
  hr->SetParameters(par);

  balBifurcationDiagram * bifd = balBifurcationDiagram::Create();
  bifd->SetDynamicalSystem(hr);
  bifd->SetFilename("hr-basin.h5");
  bifd->GetODESolver()->SetIntegrationMode(balTRAJ);
  bifd->GetODESolver()->HaltAtEquilibrium(true);
  bifd->GetODESolver()->HaltAtCycle(false);
  bifd->GetODESolver()->SetTransientDuration(0e3);
  bifd->GetODESolver()->SetFinalTime(1e3);
  //bifd->GetODESolver()->SetMaxNumberOfIntersections(200);

  bifd->SetNumberOfThreads(argc > 1 ? atoi(argv[1]) : 2);

  int nX0 = 100;
  double **X0 = new double*[nX0];
  for(int i=0; i<nX0; i++) {
    X0[i] = new double[3];
    for(int j=0; j<3; j++)
      X0[i][j] = -1 + 2*((double) random()/RAND_MAX);
  }
  bifd->SetMode(balIC);
  bifd->SetInitialConditions(nX0,X0);
  bifd->ComputeDiagram();
  bifd->Destroy();
  hr->Destroy();
  par->Destroy();
  
  return 0;
}

