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
using namespace bal;

// TEST BifurcationDiagram
int main(int argc, char *argv[]) {

  Parameters * par = Parameters::Create();
  par->SetNumber(4);
  par->At(0) = 2.88;
  par->At(1) = 2.6;
  par->At(2) = 0.01;
  par->At(3) = 4.0;

  HindmarshRose * hr = HindmarshRose::Create();
  hr->SetParameters(par);

  BifurcationDiagram * bifd = BifurcationDiagram::Create();
  bifd->SetDynamicalSystem(hr);
  bifd->SetFilename("hr-basin.h5");
  bifd->GetODESolver()->SetIntegrationMode(EVENTS);
  bifd->GetODESolver()->HaltAtEquilibrium(true);
  bifd->GetODESolver()->HaltAtCycle(true);
  bifd->GetODESolver()->SetTransientDuration(0e3);
  bifd->GetODESolver()->SetFinalTime(5e3);
  bifd->GetODESolver()->SetTimeStep(0.1);
  bifd->GetODESolver()->SetMaxNumberOfIntersections(300);

  bifd->SetNumberOfThreads(argc > 1 ? atoi(argv[1]) : 2);

  //double x0_chaos[] = {-0.882461371550183,-3.661932217696160,2.870154513826437};
  int nX0 = 1000;
  double **X0 = new double*[nX0];
  for(int i=0; i<nX0; i++) {
    X0[i] = new double[3];
    /*
    X0[i][0] = 2.5*((double) random()/RAND_MAX);
    X0[i][1] = -10 + 12*((double) random()/RAND_MAX);
    X0[i][2] = 1.5;
    */
    /*
    X0[i][0] = 2.5*((double) random()/RAND_MAX);
    X0[i][1] = -5;
    X0[i][2] = -1+3.0*((double) random()/RAND_MAX);
    */
    X0[i][0] = 1.;
    X0[i][1] = -10 + 12*((double) random()/RAND_MAX);
    X0[i][2] = -1+3.0*((double) random()/RAND_MAX);
  }
  bifd->SetMode(IC);
  bifd->SetInitialConditions(nX0,X0);
  bifd->ComputeDiagram();
  bifd->SaveSummaryData("hr-basin.classified");
  bifd->Destroy();
  hr->Destroy();
  par->Destroy();
  
  return 0;
}

