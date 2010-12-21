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

  balBaseInterp1D *interp[4];
  double *x, **y;
  int i;
  double * yy;
  double ** dyy;
 
  int n = 10;
  int nf = 1;

  x = new double [n];
  y = new double * [nf];
  for (i=0; i<nf; i++)
    y[i] = new double [n];

  for(i=0; i<n; i++) {
    x[i] = i;
    y[0][i] = i*i*i;
  }
  yy = new double [nf];
  dyy = new double * [nf];
  for (i=0; i<nf; i++)
    dyy[i] = new double[1];

  
  balLinearInterp1D *lin = balLinearInterp1D::Create();
  lin -> SetInterpolationPoints(x,y,n,nf);
  //lin -> Init();
  interp[0] = lin;

  balPolyInterp1D *pol = balPolyInterp1D::Create();
  pol -> SetInterpolationPoints(x,y,n,nf);
  pol -> SetInterpolationOrder(3);
  pol -> Init();
  interp[1] = pol;

  balSplineInterp1D *spl = balSplineInterp1D::Create();
  spl -> SetInterpolationPoints(x,y,n,nf);
  spl -> Init();
  interp[2] = spl;

  balSmoothingSplineInterp1D *smooth = balSmoothingSplineInterp1D::Create();
  smooth -> SetInterpolationPoints(x,y,n,nf);
  smooth -> SetSmoothingParameters(10);
  smooth -> Init();
  interp[3] = smooth;
 
  FILE *fout = fopen("interp.dat","w");
  for(double xx=0.1; xx<=9.9; xx+=0.1) {
    fprintf(fout,"%e", xx);
    for(i=0; i<4; i++) {
      interp[i]->Evaluate(&xx,yy);
      fprintf(fout," %e", yy[0]);
      interp[i]->EvaluateJacobian(&xx,dyy);
      fprintf(fout," %e", dyy[0][0]);
    }
    fprintf(fout,"\n");
  }
  fclose(fout);
  
  for(i=0; i<4; i++)
    interp[i]->Destroy();

  delete [] x;
  for (i=0; i<nf; i++)
    delete [] y[i];
  delete [] y;
  delete [] yy;
  for(i=0; i<nf; i++)
    delete [] dyy[i];
  delete [] dyy;
  
  return 0;
}

