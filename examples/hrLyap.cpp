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

// TEST balODESolver calculating Lyapunov Spectrum
int main(int argc, char *argv[]) {
 
  HindmarshRose hr;

  /** initial condition **/
  realtype x0[3] = {0.5,0.5,0.5};

  Parameters *pars = hr.GetParameters();
  /** parameters **/
  // chaos
  //pars->At(0) = 2.96;
  //pars->At(1) = 3.0;
  // limit cycle 
  //pars->At(0) = 4.0;
  //pars->At(1) = 5.0;
  // equilibrium
  pars->At(0) = 4.0;
  pars->At(1) = 1.0;
  pars->At(2) = 0.01;
  pars->At(3) = 4.0;
   
  ODESolver solver;
  solver.SetDynamicalSystem(&hr);
  solver.SetIntegrationMode(LYAP);
  solver.SetTimeStep(0.01);
  solver.SetTransientDuration(500);
  solver.SetLyapunovTimeStep(0.5);
  solver.SetFinalTime(2e3);
  solver.SetX0(x0);
  solver.Solve();
  
  for(int i=0; i<hr.GetOriginalDimension(); i++){
    printf("%e ",solver.GetLyapunovExponents()[i]);
  }
  printf("\n");
  
  return 0;
}

