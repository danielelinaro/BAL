/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp3D.cpp
 *
 *   Copyright (C) 2009,2010 Daniele Linaro
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

/**
 *  \file balInterp3D.cpp
 *  \brief Implementation of the classes BaseInterp3D LinearInterp3D 
 */

#include <cstring>
#include "balInterp3D.h"

namespace bal {

BaseInterp3D::BaseInterp3D() {
  nnx1 = 0;
  nnx2 = 0;
  nnx3 = 0;
  yy = NULL;
  xx1 = NULL;
  xx2 = NULL;
  xx3 = NULL;
  nnd = 3;
}

BaseInterp3D::BaseInterp3D(const BaseInterp3D &interp) : Interpolator(interp) {
  nnx1 = interp.nnx1;
  nnx2 = interp.nnx2;
  nnx3 = interp.nnx3;
  yy = interp.yy;
  xx1 = interp.xx1;
  xx2 = interp.xx2;
  xx3 = interp.xx3;
}

BaseInterp3D::~BaseInterp3D() { }

void BaseInterp3D::SetInterpolationPoints(double * xi1, double *xi2, double *xi3, double **yi, int nx1, int nx2, int nx3, int nf) {
  nnx1 = nx1;
  nnx2 = nx2;
  nnx3 = nx3;
  nnf = nf;
  yy = yi;
  xx1 = xi1;
  xx2 = xi2;
  xx3 = xi3;
}

// LinearInterp3D

LinearInterp3D::LinearInterp3D(): BaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  yloc = NULL;
}

LinearInterp3D::LinearInterp3D(const LinearInterp3D &interp):BaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = new LinearInterp1D(*interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new LinearInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = new LinearInterp2D(*interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = new LinearInterp1D(*interp.interpz);
  else
    interpz = NULL;
  if (interp.yloc != NULL) {
    int nnx1x2 = nnx1*nnx2;
    yloc = new double ** [nnx3];
    for (int i=0; i<nnx3; i++) {
      yloc[i] = new double * [nnf]; 
      for (int j=0; j<nnf; j++) {
          yloc[i][j] = new double[nnx1x2];
	  memcpy(yloc[i][j],interp.yloc[i][j],nnx1x2*sizeof(double));
      }
    }
  }
  else
    yloc = NULL;
}

LinearInterp3D::~LinearInterp3D() {
  if (x3terp != NULL)
    delete x3terp;
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      delete interpsxy[i];
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    delete interpz;
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

int LinearInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    std::cerr << "LinearInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if (x3terp != NULL)
    delete x3terp;
  x3terp = new LinearInterp1D;
  x3terp->SetInterpolationPoints(xx3,&xx3,nnx3,1);

  int i,j;
  int nnx1x2 = nnx1*nnx2;

  if (yloc != NULL){
    for (i=0; i<nnx3; i++) {
      for (j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
  
  
  yloc = new double ** [nnx3];
  for (i=0; i<nnx3; i++) {
    yloc[i] = new double * [nnf];
    for (j=0; j<nnf; j++)
      yloc[i][j] = new double[nnx1x2];
  }
  interpsxy = new LinearInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = new LinearInterp2D;
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    delete interpz;
  interpz = new LinearInterp1D;

  return 0;
}

int LinearInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    std::cerr << "LinearInterp3D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  
  int i,j;
  int idx3;
  double **yrow, *tmp;
  
  tmp = new double [nnf];
  yrow = new double * [nnf];
  for (i=0; i<nnf; i++)
    yrow[i] = new double[2];
  
  idx3 = x3terp->nextHunt() ? x3terp->Hunt(x[2]) : x3terp->Locate(x[2]);
  
  for(i=0; i<2; i++) {
    interpsxy[idx3+i]->Evaluate(x,tmp);
    for (j=0; j<nnf; j++)
      yrow[j][i] = tmp[j];
  }
  interpz->SetInterpolationPoints(&xx3[idx3],yrow,2,nnf);
  interpz->Evaluate(&x[2],y);
  delete [] tmp;
  for (i=0; i<nnf; i++)
    delete [] yrow[i];
  delete [] yrow;
  
  return 0;
}

int LinearInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "LinearInterp3D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  for (int i=0; i<nnf; i++) {
    y[i][0] = (yr[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][1] = (yf[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][2] = (yt[i]-yc[i])/FINITE_DIFFERENCES_STEP;
  }
  return 0;
}

int LinearInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "LinearInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    std::cerr << "LinearInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  y[0] = (yr[0]-yc[0])/FINITE_DIFFERENCES_STEP + (yf[1]-yc[1])/FINITE_DIFFERENCES_STEP + (yt[2]-yc[2])/FINITE_DIFFERENCES_STEP;
  
  return 0;

}


// PolyInterp3D

PolyInterp3D::PolyInterp3D(): BaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  mm1 = 2;
  mm2 = 2;
  mm3 = 2;
  yloc = NULL;
}

PolyInterp3D::PolyInterp3D(const PolyInterp3D &interp):BaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = new PolyInterp1D(*interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new PolyInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = new PolyInterp2D(*interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = new PolyInterp1D(*interp.interpz);
  else
    interpz = NULL;
  mm1 = interp.mm1;
  mm2 = interp.mm2;
  mm3 = interp.mm3;
  if (interp.yloc != NULL) {
    int nnx1x2 = nnx1*nnx2;
    yloc = new double ** [nnx3];
    for (int i=0; i<nnx3; i++) {
      yloc[i] = new double * [nnf]; 
      for (int j=0; j<nnf; j++) {
          yloc[i][j] = new double[nnx1x2];
	  memcpy(yloc[i][j],interp.yloc[i][j],nnx1x2*sizeof(double));
      }
    }
  }
  else
    yloc = NULL;
}

PolyInterp3D::~PolyInterp3D() {
  if (x3terp != NULL)
    delete x3terp;
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      delete interpsxy[i];
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    delete interpz;
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

int PolyInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    std::cerr << "PolyInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((mm1 < 0) || (mm1 > nnx1) || (mm2 < 0) || (mm2 > nnx2) || (mm3 < 0) || (mm3 > nnx3))  {
    std::cerr << "PolyInterp3D::Init() - Invalid interpolation order\n"; 
  }
  if (x3terp != NULL)
    delete x3terp;
  x3terp = new PolyInterp1D;
  x3terp->SetInterpolationPoints(xx3,&xx3,nnx3,1);
  x3terp->SetInterpolationOrder(mm3);

  int i,j;
  int nnx1x2 = nnx1*nnx2;

  if (yloc != NULL){
    for (i=0; i<nnx3; i++) {
      for (j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
  
  
  yloc = new double ** [nnx3];
  for (i=0; i<nnx3; i++) {
    yloc[i] = new double * [nnf];
    for (j=0; j<nnf; j++)
      yloc[i][j] = new double[nnx1x2];
  }
  interpsxy = new PolyInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = new PolyInterp2D;
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->SetInterpolationOrder(mm1,mm2);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    delete interpz;
  interpz = new PolyInterp1D;

  return 0;
}

void PolyInterp3D::SetInterpolationOrder(int m1, int m2, int m3) {
  mm1 = m1;
  mm2 = m2;
  mm3 = m3;
}

int PolyInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    std::cerr << "PolyInterp3D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  
  int i,j;
  int idx3;
  double **yrow, *tmp;
  
  tmp = new double [nnf];
  yrow = new double * [nnf];
  for (i=0; i<nnf; i++)
    yrow[i] = new double[mm3];
  
  idx3 = x3terp->nextHunt() ? x3terp->Hunt(x[2]) : x3terp->Locate(x[2]);
  
  for(i=0; i<mm3; i++) {
    interpsxy[idx3+i]->Evaluate(x,tmp);
    for (j=0; j<nnf; j++)
      yrow[j][i] = tmp[j];
  }
  interpz->SetInterpolationPoints(&xx3[idx3],yrow,mm3,nnf);
  interpz->SetInterpolationOrder(mm3);
  interpz->Evaluate(&x[2],y);
  delete [] tmp;
  for (i=0; i<nnf; i++)
    delete [] yrow[i];
  delete [] yrow;
  
  return 0;
}

int PolyInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "PolyInterp3D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  for (int i=0; i<nnf; i++) {
    y[i][0] = (yr[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][1] = (yf[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][2] = (yt[i]-yc[i])/FINITE_DIFFERENCES_STEP;
  }
  return 0;
}

int PolyInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "PolyInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    std::cerr << "PolyInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  y[0] = (yr[0]-yc[0])/FINITE_DIFFERENCES_STEP + (yf[1]-yc[1])/FINITE_DIFFERENCES_STEP + (yt[2]-yc[2])/FINITE_DIFFERENCES_STEP;
  
  return 0;

}

// SplineInterp3D

SplineInterp3D::SplineInterp3D(): BaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  window = 2;
  yloc = NULL;
}

SplineInterp3D::SplineInterp3D(const SplineInterp3D &interp):BaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = new PolyInterp1D(*interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new SplineInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = new SplineInterp2D(*interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = new SplineInterp1D(*interp.interpz);
  else
    interpz = NULL;
  window = interp.window;
  if (interp.yloc != NULL) {
    int nnx1x2 = nnx1*nnx2;
    yloc = new double ** [nnx3];
    for (int i=0; i<nnx3; i++) {
      yloc[i] = new double * [nnf]; 
      for (int j=0; j<nnf; j++) {
          yloc[i][j] = new double[nnx1x2];
	  memcpy(yloc[i][j],interp.yloc[i][j],nnx1x2*sizeof(double));
      }
    }
  }
  else
    yloc = NULL;
}

SplineInterp3D::~SplineInterp3D() {
  if (x3terp != NULL)
    delete x3terp;
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      delete interpsxy[i];
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    delete interpz;
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

int SplineInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    std::cerr << "SplineInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((window < 2) || (window > nnx3)) {
    std::cerr << "SplineInterp3D::Init() - Windows size must be positive and smaller than number of points along z axis\n";
    return -1;
  }
  if (x3terp != NULL)
    delete x3terp;
  x3terp = new PolyInterp1D;
  x3terp->SetInterpolationPoints(xx3,&xx3,nnx3,1);
  x3terp->SetInterpolationOrder(window);

  int i,j;
  int nnx1x2 = nnx1*nnx2;

  if (yloc != NULL){
    for (i=0; i<nnx3; i++) {
      for (j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
  
  
  yloc = new double ** [nnx3];
  for (i=0; i<nnx3; i++) {
    yloc[i] = new double * [nnf];
    for (j=0; j<nnf; j++)
      yloc[i][j] = new double[nnx1x2];
  }
  interpsxy = new SplineInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = new SplineInterp2D;
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    delete interpz;
  interpz = new SplineInterp1D;

  return 0;
}

void SplineInterp3D::SetWindow(int w) {
  window = w;
}

int SplineInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    std::cerr << "SplineInterp3D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  
  int i,j;
  int idx3;
  double **yrow, *tmp;
  
  tmp = new double [nnf];
  yrow = new double * [nnf];
  for (i=0; i<nnf; i++)
    yrow[i] = new double[window];
  
  idx3 = x3terp->nextHunt() ? x3terp->Hunt(x[2]) : x3terp->Locate(x[2]);
  
  for(i=0; i<window; i++) {
    interpsxy[idx3+i]->Evaluate(x,tmp);
    for (j=0; j<nnf; j++)
      yrow[j][i] = tmp[j];
  }
  interpz->SetInterpolationPoints(&xx3[idx3],yrow,window,nnf);
  interpz->Init();
  interpz->Evaluate(&x[2],y);
  delete [] tmp;
  for (i=0; i<nnf; i++)
    delete [] yrow[i];
  delete [] yrow;
  
  return 0;
}

int SplineInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "SplineInterp3D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  for (int i=0; i<nnf; i++) {
    y[i][0] = (yr[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][1] = (yf[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][2] = (yt[i]-yc[i])/FINITE_DIFFERENCES_STEP;
  }
  return 0;
}

int SplineInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "SplineInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    std::cerr << "SplineInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  y[0] = (yr[0]-yc[0])/FINITE_DIFFERENCES_STEP + (yf[1]-yc[1])/FINITE_DIFFERENCES_STEP + (yt[2]-yc[2])/FINITE_DIFFERENCES_STEP;
  
  return 0;

}

// SmoothingSplineInterp3D

SmoothingSplineInterp3D::SmoothingSplineInterp3D(): BaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  window = 2;
  yloc = NULL;
  SS = 0;
}

SmoothingSplineInterp3D::SmoothingSplineInterp3D(const SmoothingSplineInterp3D &interp):BaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = new PolyInterp1D(*interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new SmoothingSplineInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = new SmoothingSplineInterp2D(*interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = new SmoothingSplineInterp1D(*interp.interpz);
  else
    interpz = NULL;
  window = interp.window;
  if (interp.yloc != NULL) {
    int nnx1x2 = nnx1*nnx2;
    yloc = new double ** [nnx3];
    for (int i=0; i<nnx3; i++) {
      yloc[i] = new double * [nnf]; 
      for (int j=0; j<nnf; j++) {
          yloc[i][j] = new double[nnx1x2];
	  memcpy(yloc[i][j],interp.yloc[i][j],nnx1x2*sizeof(double));
      }
    }
  }
  else
    yloc = NULL;
  SS = interp.SS;
}

SmoothingSplineInterp3D::~SmoothingSplineInterp3D() {
  if (x3terp != NULL)
    delete x3terp;
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      delete interpsxy[i];
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    delete interpz;
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

int SmoothingSplineInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    std::cerr << "SmoothingSplineInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((window < 2) || (window > nnx3)) {
    std::cerr << "SmoothingSplineInterp3D::Init() - Windows size must be positive and < than number of points along third dimension\n";
    return -1;
  }
  if (x3terp != NULL)
    delete x3terp;
  x3terp = new PolyInterp1D;
  x3terp->SetInterpolationPoints(xx3,&xx3,nnx3,1);
  x3terp->SetInterpolationOrder(window);

  int i,j;
  int nnx1x2 = nnx1*nnx2;

  if (yloc != NULL){
    for (i=0; i<nnx3; i++) {
      for (j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
  
  
  yloc = new double ** [nnx3];
  for (i=0; i<nnx3; i++) {
    yloc[i] = new double * [nnf];
    for (j=0; j<nnf; j++)
      yloc[i][j] = new double[nnx1x2];
  }
  interpsxy = new SmoothingSplineInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = new SmoothingSplineInterp2D;
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->SetSmoothingParameter(SS);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    delete interpz;
  interpz = new SmoothingSplineInterp1D;

  return 0;
}

void SmoothingSplineInterp3D::SetWindow(int w) {
  window = w;
}

void SmoothingSplineInterp3D::SetSmoothingParameter(double S) {
  SS = S;
}

int SmoothingSplineInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    std::cerr << "SmoothingSplineInterp3D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  
  int i,j;
  int idx3;
  double **yrow, *tmp;
  
  tmp = new double [nnf];
  yrow = new double * [nnf];
  for (i=0; i<nnf; i++)
    yrow[i] = new double[window];
  
  idx3 = x3terp->nextHunt() ? x3terp->Hunt(x[2]) : x3terp->Locate(x[2]);
  
  for(i=0; i<window; i++) {
    interpsxy[idx3+i]->Evaluate(x,tmp);
    for (j=0; j<nnf; j++)
      yrow[j][i] = tmp[j];
  }
  interpz->SetInterpolationPoints(&xx3[idx3],yrow,window,nnf);
  interpz->SetSmoothingParameters(SS);
  interpz->Init();
  interpz->Evaluate(&x[2],y);
  delete [] tmp;
  for (i=0; i<nnf; i++)
    delete [] yrow[i];
  delete [] yrow;
  
  return 0;
}

int SmoothingSplineInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "SmoothingSplineInterp3D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  for (int i=0; i<nnf; i++) {
    y[i][0] = (yr[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][1] = (yf[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][2] = (yt[i]-yc[i])/FINITE_DIFFERENCES_STEP;
  }
  return 0;
}

int SmoothingSplineInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "SmoothingSplineInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    std::cerr << "SmoothingSplineInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  double xr[3],xf[3],xt[3],yc[nnf],yr[nnf],yf[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xr[2] = x[2];
  xf[0] = x[0];
  xf[1] = x[1]+FINITE_DIFFERENCES_STEP;
  xf[2] = x[2];
  xt[0] = x[0];
  xt[1] = x[1];
  xt[2] = x[2]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xf,yf);
  Evaluate(xt,yt);

  y[0] = (yr[0]-yc[0])/FINITE_DIFFERENCES_STEP + (yf[1]-yc[1])/FINITE_DIFFERENCES_STEP + (yt[2]-yc[2])/FINITE_DIFFERENCES_STEP;
  
  return 0;

}

} // namespace bal

