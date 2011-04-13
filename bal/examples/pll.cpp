/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    main.cpp
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

#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <nvector/nvector_serial.h>
#include "balCommon.h"
#include "balParameters.h"
#include "balLogger.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"
#include "balBifurcationParameters.h"
#include "balPLL.h"
#include "ConfigFile.h"
using namespace bal;

const realtype pi = 3.1415926535897931;

int main(int argc, char *argv[]) {
  if(argc == 1) {
    fprintf(stderr, "Usage: %s ConfigFile\n", argv[0]);
    exit(1);
  }
  ConfigFile config(argv[1]);
  char name[100];
  
  int steps[PLL::npar];
  BifurcationParameters * bp = BifurcationParameters::Create();
  bp->SetNumber(PLL::npar);
  
  for(int i=0; i<PLL::npar; i++) {
    sprintf(name,"%ssteps",PLL::parname[i]);
    steps[i] = config.read<int>(name);
    if(steps[i] == 1) {
      bp->SetIthParameter(i,config.read<double>(PLL::parname[i]));
    }
    else {
      sprintf(name,"%smin",PLL::parname[i]);
      bp->SetIthParameterLowerBound(i,config.read<double>(name));
      sprintf(name,"%smax",PLL::parname[i]);
      bp->SetIthParameterUpperBound(i,config.read<double>(name));
    }
  }
  
#ifdef WITHPHIERR
#ifdef EXTEND
  realtype x0[7] = {0.,0.,0.,0.,pi,0.,0.};
#else
  realtype x0[5] = {0.,0.,0.,0.,pi};
#endif
#else
#ifdef EXTEND
  realtype x0[6] = {0.,0.,0.,0.,0.,0.};
#else
  realtype x0[4] = {0.,0.,0.,0.};
#endif
#endif
  x0[1] = config.read<double>("vdd");	
  
  bp->SetNumberOfSteps(steps);
  PLL * pll = PLL::Create();
  pll->SetParameters(bp);
  BifurcationDiagram * bifd = BifurcationDiagram::Create();
  bifd->RestartFromX0(true);
  bifd->SetDynamicalSystem(pll);
  bifd->SetFilename((char *) config.read<string>("outputfile").c_str());
  if(config.read<bool>("trajectory"))
    bifd->GetODESolver()->SetIntegrationMode(balBOTH);
  else
    bifd->GetODESolver()->SetIntegrationMode(balEVENTS);
  bifd->GetODESolver()->SetTransientDuration(config.read<double>("ttran"));
  bifd->GetODESolver()->HaltAtEquilibrium(false);
  bifd->GetODESolver()->HaltAtCycle(false);
  bifd->GetODESolver()->SetFinalTime(config.read<double>("tout"));
  bifd->GetODESolver()->SetTimeStep(1e-11);
  bifd->GetODESolver()->SetMaxNumberOfIntersections((int) 1e7);
  bifd->GetODESolver()->SetX0(x0);
  bifd->GetODESolver()->SetRelativeTolerance(1e-10);
  bifd->SetNumberOfThreads(1);
  bifd->ComputeDiagram();
  bifd->Destroy();
  pll->Destroy();
  bp->Destroy();
  
  return EXIT_SUCCESS;
}

