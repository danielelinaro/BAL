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

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <nvector/nvector_serial.h>
#include "balCommon.h"
#include "balObject.h"
#include "balDynamicalSystem.h"
#include "balHindmarshRose.h"
#include "balParameters.h"
#include "balSolution.h"
#include "balLogger.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"
#include "balBifurcationParameters.h"
#include "balHeartNeuron.h"
#include "balPLL.h"
#include "balInterp1D.h"
using namespace std;

int main(int argc, char *argv[]) {

  balBaseInterp1D *interp[3];
  double x[10], y[10];
  int i;
  
  for(i=0; i<10; i++) {
    x[i] = i;
    y[i] = i*i;
  }
  
  interp[0] = balLinearInterp1D::Create(x,y,10);
  interp[1] = balPolyInterp1D::Create(x,y,10,3);
  interp[2] = balSplineInterp1D::Create(x,y,10);
  
  for(double xx=0.1; xx<=9.9; xx+=0.1) {
    printf("%e", xx);
    for(i=0; i<3; i++)
      printf(" %e", interp[i]->interp(xx));  
    printf("\n");
  }
  
  for(i=0; i<3; i++)
    interp[i]->Destroy();
  
  return 0;
}

