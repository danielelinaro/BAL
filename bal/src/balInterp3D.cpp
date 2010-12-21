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
 *  \brief Implementation of the classes balBaseInterp3D balLinearInterp3D 
 */

#include "balInterp3D.h"

balBaseInterp3D::balBaseInterp3D() {
  nnx1 = 0;
  nnx2 = 0;
  nnx3 = 0;
  yy = NULL;
  xx1 = NULL;
  xx2 = NULL;
  xx3 = NULL;
  nnd = 3;
}

balBaseInterp3D::balBaseInterp3D(const balBaseInterp3D &interp) : balInterpolator(interp) {
  nnx1 = interp.nnx1;
  nnx2 = interp.nnx2;
  nnx3 = interp.nnx3;
  yy = interp.yy;
  xx1 = interp.xx1;
  xx2 = interp.xx2;
  xx3 = interp.xx3;
}

const char * balBaseInterp3D::GetClassName() const {
  return "balBaseInterp3D";
}

void balBaseInterp3D::Destroy() {
  delete this;
}


balBaseInterp3D::~balBaseInterp3D() { }

void balBaseInterp3D::SetInterpolationPoints(double * xi1, double *xi2, double *xi3, double **yi, int nx1, int nx2, int nx3, int nf) {
  nnx1 = nx1;
  nnx2 = nx2;
  nnx3 = nx3;
  nnf = nf;
  yy = yi;
  xx1 = xi1;
  xx2 = xi2;
  xx3 = xi3;
}

// balLinearInterp3D

balLinearInterp3D::balLinearInterp3D(): balBaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  yloc = NULL;
}

balLinearInterp3D::balLinearInterp3D(const balLinearInterp3D &interp):balBaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = balLinearInterp1D::Copy(interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new balLinearInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = balLinearInterp2D::Copy(interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = balLinearInterp1D::Copy(interp.interpz);
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

const char * balLinearInterp3D::GetClassName() const {
  return "balLinearInterp3D";
}

void balLinearInterp3D::Destroy() {
  delete this;
}

balLinearInterp3D::~balLinearInterp3D() {
  if (x3terp != NULL)
    x3terp->Destroy();
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      interpsxy[i]->Destroy();
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    interpz->Destroy();
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

balLinearInterp3D * balLinearInterp3D::Create() {
  return new balLinearInterp3D();
}

balLinearInterp3D * balLinearInterp3D::Copy(balLinearInterp3D *interp) {
  return new balLinearInterp3D(*interp);
}

balLinearInterp3D * balLinearInterp3D::Clone() const {
  return new balLinearInterp3D(*this);
}

int balLinearInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    cerr<<"balLinearInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if (x3terp != NULL)
    x3terp->Destroy();
  x3terp = balLinearInterp1D::Create();
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
  interpsxy = new balLinearInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = balLinearInterp2D::Create();
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    interpz->Destroy();
  interpz = balLinearInterp1D::Create();

  return 0;
}

int balLinearInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    cerr<<"balLinearInterp3D::Evaluate() - Interpolation points not set\n";
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

int balLinearInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balLinearInterp3D::EvaluateJacobian() - Interpolation points not set\n";
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

int balLinearInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balLinearInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    cerr<<"balLinearInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
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


// balPolyInterp3D

balPolyInterp3D::balPolyInterp3D(): balBaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  mm1 = 2;
  mm2 = 2;
  mm3 = 2;
  yloc = NULL;
}

balPolyInterp3D::balPolyInterp3D(const balPolyInterp3D &interp):balBaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = balPolyInterp1D::Copy(interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new balPolyInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = balPolyInterp2D::Copy(interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = balPolyInterp1D::Copy(interp.interpz);
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

const char * balPolyInterp3D::GetClassName() const {
  return "balPolyInterp3D";
}

void balPolyInterp3D::Destroy() {
  delete this;
}

balPolyInterp3D::~balPolyInterp3D() {
  if (x3terp != NULL)
    x3terp->Destroy();
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      interpsxy[i]->Destroy();
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    interpz->Destroy();
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

balPolyInterp3D * balPolyInterp3D::Create() {
  return new balPolyInterp3D();
}

balPolyInterp3D * balPolyInterp3D::Copy(balPolyInterp3D *interp) {
  return new balPolyInterp3D(*interp);
}

balPolyInterp3D * balPolyInterp3D::Clone() const {
  return new balPolyInterp3D(*this);
}

int balPolyInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    cerr<<"balPolyInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((mm1 < 0) || (mm1 > nnx1) || (mm2 < 0) || (mm2 > nnx2) || (mm3 < 0) || (mm3 > nnx3))  {
    cerr<<"balPolyInterp3D::Init() - Invalid interpolation order\n"; 
  }
  if (x3terp != NULL)
    x3terp->Destroy();
  x3terp = balPolyInterp1D::Create();
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
  interpsxy = new balPolyInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = balPolyInterp2D::Create();
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->SetInterpolationOrder(mm1,mm2);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    interpz->Destroy();
  interpz = balPolyInterp1D::Create();

  return 0;
}

void balPolyInterp3D::SetInterpolationOrder(int m1, int m2, int m3) {
  mm1 = m1;
  mm2 = m2;
  mm3 = m3;
}

int balPolyInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    cerr<<"balPolyInterp3D::Evaluate() - Interpolation points not set\n";
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

int balPolyInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balPolyInterp3D::EvaluateJacobian() - Interpolation points not set\n";
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

int balPolyInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balPolyInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    cerr<<"balPolyInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
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

// balSplineInterp3D

balSplineInterp3D::balSplineInterp3D(): balBaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  window = 2;
  yloc = NULL;
}

balSplineInterp3D::balSplineInterp3D(const balSplineInterp3D &interp):balBaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = balPolyInterp1D::Copy(interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new balSplineInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = balSplineInterp2D::Copy(interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = balSplineInterp1D::Copy(interp.interpz);
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

const char * balSplineInterp3D::GetClassName() const {
  return "balSplineInterp3D";
}

void balSplineInterp3D::Destroy() {
  delete this;
}

balSplineInterp3D::~balSplineInterp3D() {
  if (x3terp != NULL)
    x3terp->Destroy();
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      interpsxy[i]->Destroy();
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    interpz->Destroy();
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

balSplineInterp3D * balSplineInterp3D::Create() {
  return new balSplineInterp3D();
}

balSplineInterp3D * balSplineInterp3D::Copy(balSplineInterp3D *interp) {
  return new balSplineInterp3D(*interp);
}

balSplineInterp3D * balSplineInterp3D::Clone() const {
  return new balSplineInterp3D(*this);
}

int balSplineInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    cerr<<"balSplineInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((window < 2) || (window > nnx3)) {
    cerr<<"balSplineInterp3D::Init() - Windows size must be positive and smaller than number of points along z axis\n";
    return -1;
  }
  if (x3terp != NULL)
    x3terp->Destroy();
  x3terp = balPolyInterp1D::Create();
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
  interpsxy = new balSplineInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = balSplineInterp2D::Create();
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    interpz->Destroy();
  interpz = balSplineInterp1D::Create();

  return 0;
}

void balSplineInterp3D::SetWindow(int w) {
  window = w;
}

int balSplineInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    cerr<<"balSplineInterp3D::Evaluate() - Interpolation points not set\n";
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

int balSplineInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balSplineInterp3D::EvaluateJacobian() - Interpolation points not set\n";
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

int balSplineInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balSplineInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    cerr<<"balSplineInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
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

// balSmoothingSplineInterp3D

balSmoothingSplineInterp3D::balSmoothingSplineInterp3D(): balBaseInterp3D() {
  x3terp = NULL;
  interpsxy = NULL;
  interpz = NULL;
  window = 2;
  yloc = NULL;
  SS = 0;
}

balSmoothingSplineInterp3D::balSmoothingSplineInterp3D(const balSmoothingSplineInterp3D &interp):balBaseInterp3D(interp) {
  if (interp.x3terp!=NULL)
    x3terp = balPolyInterp1D::Copy(interp.x3terp);
  else
    x3terp = NULL;
  if (interp.interpsxy != NULL) {
    interpsxy = new balSmoothingSplineInterp2D * [nnx3];
    for (int i=0; i<nnx3; i++) 
      interpsxy[i] = balSmoothingSplineInterp2D::Copy(interp.interpsxy[i]);
  }
  if (interp.interpz!=NULL)
    interpz = balSmoothingSplineInterp1D::Copy(interp.interpz);
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

const char * balSmoothingSplineInterp3D::GetClassName() const {
  return "balSmoothingSplineInterp3D";
}

void balSmoothingSplineInterp3D::Destroy() {
  delete this;
}

balSmoothingSplineInterp3D::~balSmoothingSplineInterp3D() {
  if (x3terp != NULL)
    x3terp->Destroy();
  if (interpsxy != NULL) {
    for (int i=0; i<nnx3; i++)
      interpsxy[i]->Destroy();
    delete [] interpsxy;
  }
  if (interpz != NULL) 
    interpz->Destroy();
  if (yloc != NULL){
    for (int i=0; i<nnx3; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }

}

balSmoothingSplineInterp3D * balSmoothingSplineInterp3D::Create() {
  return new balSmoothingSplineInterp3D();
}

balSmoothingSplineInterp3D * balSmoothingSplineInterp3D::Copy(balSmoothingSplineInterp3D *interp) {
  return new balSmoothingSplineInterp3D(*interp);
}

balSmoothingSplineInterp3D * balSmoothingSplineInterp3D::Clone() const {
  return new balSmoothingSplineInterp3D(*this);
}

int balSmoothingSplineInterp3D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(xx3==NULL)||(yy==NULL)) {
    cerr<<"balSmoothingSplineInterp3D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((window < 2) || (window > nnx3)) {
    cerr<<"balSmoothingSplineInterp3D::Init() - Windows size must be positive and < than number of points along third dimension\n";
    return -1;
  }
  if (x3terp != NULL)
    x3terp->Destroy();
  x3terp = balPolyInterp1D::Create();
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
  interpsxy = new balSmoothingSplineInterp2D * [nnx3];

  for (i = 0; i<nnx3; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1x2],nnx1x2*sizeof(double));
     }
     interpsxy[i] = balSmoothingSplineInterp2D::Create();
     interpsxy[i]->SetInterpolationPoints(xx1,xx2,yloc[i],nnx1,nnx2,nnf);
     interpsxy[i]->SetSmoothingParameter(SS);
     interpsxy[i]->Init();
  } 
  if (interpz != NULL)
    interpz->Destroy();
  interpz = balSmoothingSplineInterp1D::Create();

  return 0;
}

void balSmoothingSplineInterp3D::SetWindow(int w) {
  window = w;
}

void balSmoothingSplineInterp3D::SetSmoothingParameter(double S) {
  SS = S;
}

int balSmoothingSplineInterp3D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    cerr<<"balSmoothingSplineInterp3D::Evaluate() - Interpolation points not set\n";
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

int balSmoothingSplineInterp3D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balSmoothingSplineInterp3D::EvaluateJacobian() - Interpolation points not set\n";
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

int balSmoothingSplineInterp3D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    cerr<<"balSmoothingSplineInterp3D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 3) {
    cerr<<"balSmoothingSplineInterp3D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
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
