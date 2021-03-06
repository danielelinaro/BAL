/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp2D.cpp
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
 *  \file balInterp2D.cpp
 *  \brief Implementation of the classes BaseInterp2D LinearInterp2D 
 */

#include "balInterp2D.h"

namespace bal {

BaseInterp2D::BaseInterp2D() {
  nnx1 = 0;
  nnx2 = 0;
  yy = NULL;
  xx1 = NULL;
  xx2 = NULL;
  nnd = 2;
}

BaseInterp2D::BaseInterp2D(const BaseInterp2D &interp) : Interpolator(interp) {
  nnx1 = interp.nnx1;
  nnx2 = interp.nnx2;
  yy = interp.yy;
  xx1 = interp.xx1;
  xx2 = interp.xx2;
}

const char * BaseInterp2D::GetClassName() const {
  return "BaseInterp2D";
}

void BaseInterp2D::Destroy() {
  delete this;
}


BaseInterp2D::~BaseInterp2D() { }

void BaseInterp2D::SetInterpolationPoints(double * xi1, double *xi2, double **yi, int nx1, int nx2, int nf) {
  nnx1 = nx1;
  nnx2 = nx2;
  nnf = nf;
  yy = yi;
  xx1 = xi1;
  xx2 = xi2;
}


// LinearInterp2D

LinearInterp2D::LinearInterp2D(): BaseInterp2D() {
  x1terp = NULL;
  x2terp = NULL;
}

LinearInterp2D::LinearInterp2D(const LinearInterp2D &interp):BaseInterp2D(interp) {
  if (interp.x1terp!=NULL)
    x1terp = LinearInterp1D::Copy(interp.x1terp);
  else
    x1terp = NULL;
  if (interp.x2terp!=NULL)
    x2terp = LinearInterp1D::Copy(interp.x2terp);
  else
    x2terp = NULL;
}

const char * LinearInterp2D::GetClassName() const {
  return "LinearInterp2D";
}

void LinearInterp2D::Destroy() {
  delete this;
}

LinearInterp2D::~LinearInterp2D() {
  if (x1terp != NULL)
    x1terp->Destroy();
  if (x2terp != NULL)
    x2terp->Destroy();
}

LinearInterp2D * LinearInterp2D::Create() {
  return new LinearInterp2D();
}

LinearInterp2D * LinearInterp2D::Copy(LinearInterp2D *interp) {
  return new LinearInterp2D(*interp);
}

LinearInterp2D * LinearInterp2D::Clone() const {
  return new LinearInterp2D(*this);
}

int LinearInterp2D::Init() {
  if ((xx1==NULL)||(xx2==NULL)||(yy==NULL)) {
    std::cerr << "LinearInterp2D::Init() - Interpolation points not set\n";
    return -1;
  }
  if (x1terp != NULL)
    x1terp->Destroy();
  if (x2terp != NULL)
    x2terp->Destroy();
  x1terp = LinearInterp1D::Create();
  x2terp = LinearInterp1D::Create();
  x1terp->SetInterpolationPoints(xx1,&xx1,nnx1,1);
  x2terp->SetInterpolationPoints(xx2,&xx2,nnx2,1);
  return 0;
}

int LinearInterp2D::Evaluate(double *x, double *y) {

  if ((x1terp == NULL)||(x2terp == NULL)) {
    std::cerr << "LinearInterp2D::Evaluate() - Interpolator not initialized. Call method Init()\n";
    return -1;
  }
  
  int i;
  int idx1, idx2;
  double t, u, d1, d2;

  idx1 = x1terp->nextHunt() ? x1terp->Hunt(x[0]) : x1terp->Locate(x[0]);
  idx2 = x2terp->nextHunt() ? x2terp->Hunt(x[1]) : x2terp->Locate(x[1]);

  if ((d1 = xx1[idx1+1]-xx1[idx1]) == 0) {
    std::cerr << "Bad input to LinearInterp2D::Evaluate()\n";
    return -1;
  }
  if ((d2 = xx2[idx2+1]-xx2[idx2]) == 0) {
    std::cerr << "Bad input to LinearInterp2D::Evaluate()\n";
    return -1;
  }

  t = (x[0]-xx1[idx1])/(d1);
  u = (x[1]-xx2[idx2])/(d2); 
  for (i=0; i<nnf; i++) {
    y[i] = (1.-t)*(1.-u)*yy[i][idx1+nnx1*idx2] + t*(1.-u)*yy[i][(idx1+1)+nnx1*idx2] + (1.-t)*u*yy[i][idx1+nnx1*(idx2+1)] + t*u*yy[i][(idx1+1)+nnx1*(idx2+1)];
  }
  return 0;
}

int LinearInterp2D::EvaluateJacobian(double *x, double **y) {
  if ((x1terp == NULL)||(x2terp == NULL)) {
    std::cerr << "LinearInterp2D::EvaluateJacobian() - Interpolator not initialized. Call method Init()\n";
    return -1;
  }
  
  int idx1, idx2;
  idx1 = x1terp->nextHunt() ? x1terp->Hunt(x[0]) : x1terp->Locate(x[0]);
  idx2 = x2terp->nextHunt() ? x2terp->Hunt(x[1]) : x2terp->Locate(x[1]);
  
  for (int i=0; i<nnf; i++) {
    y[i][0] = (yy[i][(idx1+1)+nnx1*idx2]-yy[i][idx1+nnx1*idx2])/(xx1[idx1+1]-xx1[idx1]);
    y[i][1] = (yy[i][idx1+nnx1*(idx2+1)]-yy[i][idx1+nnx1*idx2])/(xx2[idx2+1]-xx2[idx2]);
  }
  return 0;
}
  
int LinearInterp2D::EvaluateDivergence(double *x, double *y) {
  if ((x1terp == NULL)||(x2terp == NULL)) {
    std::cerr << "LinearInterp2D::EvaluateDivergence() - Interpolator not initialized. Call method Init()\n";
    return -1;
  }
  if (nnf != 2) {
    std::cerr << "LinearInterp2D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  int idx1, idx2;
  idx1 = x1terp->nextHunt() ? x1terp->Hunt(x[0]) : x1terp->Locate(x[0]);
  idx2 = x2terp->nextHunt() ? x2terp->Hunt(x[1]) : x2terp->Locate(x[1]);
  
  y[0] = (yy[0][(idx1+1)+nnx1*idx2]-yy[0][idx1+nnx1*idx2])/(xx1[idx1+1]-xx1[idx1]) + (yy[1][idx1+nnx1*(idx2+1)]-yy[1][idx1+nnx1*idx2])/(xx2[idx2+1]-xx2[idx2]);
  return 0;
}


// PolyInterp2D

PolyInterp2D::PolyInterp2D(): BaseInterp2D() {
  mm1 = 2;
  mm2 = 2;
  interpsx = NULL;
  interpy = NULL;
  x2terp = NULL;
  yloc = NULL;
}

PolyInterp2D::PolyInterp2D(const PolyInterp2D &interp):BaseInterp2D(interp) {
  mm1 = interp.mm1;
  mm2 = interp.mm2;
  if (interp.interpsx != NULL) {
    interpsx = new PolyInterp1D * [nnx2];
    for (int i=0; i<nnx2; i++) 
      interpsx[i] = PolyInterp1D::Copy(interp.interpsx[i]);
  }
  else
    interpsx = NULL;
  if (interp.interpy != NULL) {
    interpy = PolyInterp1D::Copy(interp.interpy);
  }
  else
    interpy = NULL;
  if (interp.x2terp!=NULL)
    x2terp = PolyInterp1D::Copy(interp.x2terp);
  else 
    x2terp = NULL;
  if (interp.yloc != NULL) {
    yloc = new double ** [nnx2];
    for (int i=0; i<nnx2; i++) {
      yloc[i] = new double * [nnf];
      for (int j=0; j<nnf; j++) {
        yloc[i][j] = new double [nnx1];
	memcpy(yloc[i][j],interp.yloc[i][j],nnx1*sizeof(double));
      }
    }
  }
  else
    yloc = NULL;
}

const char * PolyInterp2D::GetClassName() const {
  return "PolyInterp2D";
}

void PolyInterp2D::Destroy() {
  delete this;
}

PolyInterp2D::~PolyInterp2D() {
  if (interpsx != NULL) {
    for (int i=0; i<nnx2; i++)
      interpsx[i]->Destroy();
    delete [] interpsx;
  }
  if (interpy != NULL)
    interpy->Destroy();
  if (x2terp != NULL)
    x2terp->Destroy();
  if (yloc != NULL){
    for (int i=0; i<nnx2; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
}

PolyInterp2D * PolyInterp2D::Copy(PolyInterp2D *interp) {
  return new PolyInterp2D(*interp);
}

PolyInterp2D * PolyInterp2D::Clone() const {
  return new PolyInterp2D(*this);
}

PolyInterp2D * PolyInterp2D::Create() {
  return new PolyInterp2D();
}

void PolyInterp2D::SetInterpolationOrder(int m1, int m2) {
  mm1 = m1;
  mm2 = m2;
}

int PolyInterp2D::Init() {
  if ((xx1 == NULL)||(xx2 == NULL) || (yy == NULL)) {
    std::cerr << "PolyInterp2D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((mm1 < 0) || (mm1 > nnx1) || (mm2 < 0) || (mm2 > nnx2))  {
    std::cerr << "PolyInterp2D::Init() - Invalid interpolation order\n"; 
  }
  int i,j;

  if (yloc != NULL){
    for (int i=0; i<nnx2; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
  
  yloc = new double ** [nnx2];
  for (i=0; i<nnx2; i++) {
    yloc[i] = new double * [nnf];
    for (j=0; j<nnf; j++)
      yloc[i][j] = new double[nnx1];
  }
  interpsx = new PolyInterp1D * [nnx2];
  for (i = 0; i<nnx2; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1],nnx1*sizeof(double));
     }
     interpsx[i] = PolyInterp1D::Create();
     interpsx[i]->SetInterpolationPoints(xx1,yloc[i],nnx1,nnf);
     interpsx[i]->SetInterpolationOrder(mm1);
  } 
  if (interpy != NULL)
    interpy->Destroy();
  interpy = PolyInterp1D::Create();
  interpy->SetInterpolationOrder(mm2);
  
  if (x2terp != NULL)
    x2terp->Destroy();
  
  x2terp = PolyInterp1D::Create();
  x2terp->SetInterpolationPoints(xx2,&xx2,nnx2,1);
  x2terp->SetInterpolationOrder(mm2);
  return 0;
}

int PolyInterp2D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    std::cerr << "PolyInterp2D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  
  int i,j;
  int idx2;
  double **yrow, *tmp;
  
  tmp = new double [nnf];
  yrow = new double * [nnf];
  for (i=0; i<nnf; i++)
    yrow[i] = new double[mm1];
  
  idx2 = x2terp->nextHunt() ? x2terp->Hunt(x[1]) : x2terp->Locate(x[1]);
  
  for(i=0; i<mm2; i++) {
    interpsx[idx2+i]->Evaluate(&x[0],tmp);
    for (j=0; j<nnf; j++)
      yrow[j][i] = tmp[j];
  }
  interpy->SetInterpolationPoints(&xx2[idx2],yrow,mm2,nnf);
  interpy->Evaluate(&x[1],y);
  delete [] tmp;
  for (i=0; i<nnf; i++)
    delete [] yrow[i];
  delete [] yrow;
  return 0;
}

int PolyInterp2D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  
  if (xx1 == NULL) {
    std::cerr << "PolyInterp2D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }

  double xr[2],xt[2],yc[nnf],yr[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xt[0] = x[0];
  xt[1] = x[1]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xt,yt);

  for (int i=0; i<nnf; i++) {
    y[i][0] = (yr[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][1] = (yt[i]-yc[i])/FINITE_DIFFERENCES_STEP;
  }

  return 0;
}

int PolyInterp2D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "PolyInterp2D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 2) {
    std::cerr << "PolyInterp2D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  double xr[2],xt[2],yc[nnf],yr[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xt[0] = x[0];
  xt[1] = x[1]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xt,yt);

  y[0] = (yr[0]-yc[0])/FINITE_DIFFERENCES_STEP + (yt[1]-yc[1])/FINITE_DIFFERENCES_STEP;

  return 0;
}

// SplineInterp2D

SplineInterp2D::SplineInterp2D(): BaseInterp2D() {
  x1terp = NULL;
  x2terp = NULL;
  y1d = NULL;
  y2d = NULL;
  y12d = NULL;
  c = NULL;
  
  wt[0][0] = 1;
  wt[0][1] = 0;
  wt[0][2] = 0;
  wt[0][3] = 0;
  wt[0][4] = 0;
  wt[0][5] = 0;
  wt[0][6] = 0;
  wt[0][7] = 0;
  wt[0][8] = 0;
  wt[0][9] = 0;
  wt[0][10] = 0;
  wt[0][11] = 0;
  wt[0][12] = 0;
  wt[0][13] = 0;
  wt[0][14] = 0;
  wt[0][15] = 0;
  wt[1][0] = 0;
  wt[1][1] = 0;
  wt[1][2] = 0;
  wt[1][3] = 0;
  wt[1][4] = 0;
  wt[1][5] = 0;
  wt[1][6] = 0;
  wt[1][7] = 0;
  wt[1][8] = 1;
  wt[1][9] = 0;
  wt[1][10] = 0;
  wt[1][11] = 0;
  wt[1][12] = 0;
  wt[1][13] = 0;
  wt[1][14] = 0;
  wt[1][15] = 0;
  wt[2][0] = -3;
  wt[2][1] = 0;
  wt[2][2] = 0;
  wt[2][3] = 3;
  wt[2][4] = 0;
  wt[2][5] = 0;
  wt[2][6] = 0;
  wt[2][7] = 0;
  wt[2][8] = -2;
  wt[2][9] = 0;
  wt[2][10] = 0;
  wt[2][11] = -1;
  wt[2][12] = 0;
  wt[2][13] = 0;
  wt[2][14] = 0;
  wt[2][15] = 0;
  wt[3][0] = 2;
  wt[3][1] = 0;
  wt[3][2] = 0;
  wt[3][3] = -2;
  wt[3][4] = 0;
  wt[3][5] = 0;
  wt[3][6] = 0;
  wt[3][7] = 0;
  wt[3][8] = 1;
  wt[3][9] = 0;
  wt[3][10] = 0;
  wt[3][11] = 1;
  wt[3][12] = 0;
  wt[3][13] = 0;
  wt[3][14] = 0;
  wt[3][15] = 0;
  wt[4][0] = 0;
  wt[4][1] = 0;
  wt[4][2] = 0;
  wt[4][3] = 0;
  wt[4][4] = 1;
  wt[4][5] = 0;
  wt[4][6] = 0;
  wt[4][7] = 0;
  wt[4][8] = 0;
  wt[4][9] = 0;
  wt[4][10] = 0;
  wt[4][11] = 0;
  wt[4][12] = 0;
  wt[4][13] = 0;
  wt[4][14] = 0;
  wt[4][15] = 0;
  wt[5][0] = 0;
  wt[5][1] = 0;
  wt[5][2] = 0;
  wt[5][3] = 0;
  wt[5][4] = 0;
  wt[5][5] = 0;
  wt[5][6] = 0;
  wt[5][7] = 0;
  wt[5][8] = 0;
  wt[5][9] = 0;
  wt[5][10] = 0;
  wt[5][11] = 0;
  wt[5][12] = 1;
  wt[5][13] = 0;
  wt[5][14] = 0;
  wt[5][15] = 0;
  wt[6][0] = 0;
  wt[6][1] = 0;
  wt[6][2] = 0;
  wt[6][3] = 0;
  wt[6][4] = -3;
  wt[6][5] = 0;
  wt[6][6] = 0;
  wt[6][7] = 3;
  wt[6][8] = 0;
  wt[6][9] = 0;
  wt[6][10] = 0;
  wt[6][11] = 0;
  wt[6][12] = -2;
  wt[6][13] = 0;
  wt[6][14] = 0;
  wt[6][15] = -1;
  wt[7][0] = 0;
  wt[7][1] = 0;
  wt[7][2] = 0;
  wt[7][3] = 0;
  wt[7][4] = 2;
  wt[7][5] = 0;
  wt[7][6] = 0;
  wt[7][7] = -2;
  wt[7][8] = 0;
  wt[7][9] = 0;
  wt[7][10] = 0;
  wt[7][11] = 0;
  wt[7][12] = 1;
  wt[7][13] = 0;
  wt[7][14] = 0;
  wt[7][15] = 1;
  wt[8][0] = -3;
  wt[8][1] = 3;
  wt[8][2] = 0;
  wt[8][3] = 0;
  wt[8][4] = -2;
  wt[8][5] = -1;
  wt[8][6] = 0;
  wt[8][7] = 0;
  wt[8][8] = 0;
  wt[8][9] = 0;
  wt[8][10] = 0;
  wt[8][11] = 0;
  wt[8][12] = 0;
  wt[8][13] = 0;
  wt[8][14] = 0;
  wt[8][15] = 0;
  wt[9][0] = 0;
  wt[9][1] = 0;
  wt[9][2] = 0;
  wt[9][3] = 0;
  wt[9][4] = 0;
  wt[9][5] = 0;
  wt[9][6] = 0;
  wt[9][7] = 0;
  wt[9][8] = -3;
  wt[9][9] = 3;
  wt[9][10] = 0;
  wt[9][11] = 0;
  wt[9][12] = -2;
  wt[9][13] = -1;
  wt[9][14] = 0;
  wt[9][15] = 0;
  wt[10][0] = 9;
  wt[10][1] = -9;
  wt[10][2] = 9;
  wt[10][3] = -9;
  wt[10][4] = 6;
  wt[10][5] = 3;
  wt[10][6] = -3;
  wt[10][7] = -6;
  wt[10][8] = 6;
  wt[10][9] = -6;
  wt[10][10] = -3;
  wt[10][11] = 3;
  wt[10][12] = 4;
  wt[10][13] = 2;
  wt[10][14] = 1;
  wt[10][15] = 2;
  wt[11][0] = -6;
  wt[11][1] = 6;
  wt[11][2] = -6;
  wt[11][3] = 6;
  wt[11][4] = -4;
  wt[11][5] = -2;
  wt[11][6] = 2;
  wt[11][7] = 4;
  wt[11][8] = -3;
  wt[11][9] = 3;
  wt[11][10] = 3;
  wt[11][11] = -3;
  wt[11][12] = -2;
  wt[11][13] = -1;
  wt[11][14] = -1;
  wt[11][15] = -2;
  wt[12][0] = 2;
  wt[12][1] = -2;
  wt[12][2] = 0;
  wt[12][3] = 0;
  wt[12][4] = 1;
  wt[12][5] = 1;
  wt[12][6] = 0;
  wt[12][7] = 0;
  wt[12][8] = 0;
  wt[12][9] = 0;
  wt[12][10] = 0;
  wt[12][11] = 0;
  wt[12][12] = 0;
  wt[12][13] = 0;
  wt[12][14] = 0;
  wt[12][15] = 0;
  wt[13][0] = 0;
  wt[13][1] = 0;
  wt[13][2] = 0;
  wt[13][3] = 0;
  wt[13][4] = 0;
  wt[13][5] = 0;
  wt[13][6] = 0;
  wt[13][7] = 0;
  wt[13][8] = 2;
  wt[13][9] = -2;
  wt[13][10] = 0;
  wt[13][11] = 0;
  wt[13][12] = 1;
  wt[13][13] = 1;
  wt[13][14] = 0;
  wt[13][15] = 0;
  wt[14][0] = -6;
  wt[14][1] = 6;
  wt[14][2] = -6;
  wt[14][3] = 6;
  wt[14][4] = -3;
  wt[14][5] = -3;
  wt[14][6] = 3;
  wt[14][7] = 3;
  wt[14][8] = -4;
  wt[14][9] = 4;
  wt[14][10] = 2;
  wt[14][11] = -2;
  wt[14][12] = -2;
  wt[14][13] = -2;
  wt[14][14] = -1;
  wt[14][15] = -1;
  wt[15][0] = 4;
  wt[15][1] = -4;
  wt[15][2] = 4;
  wt[15][3] = -4;
  wt[15][4] = 2;
  wt[15][5] = 2;
  wt[15][6] = -2;
  wt[15][7] = -2;
  wt[15][8] = 2;
  wt[15][9] = -2;
  wt[15][10] = -2;
  wt[15][11] = 2;
  wt[15][12] = 1;
  wt[15][13] = 1;
  wt[15][14] = 1;
  wt[15][15] = 1;
}

SplineInterp2D::SplineInterp2D(const SplineInterp2D &interp):BaseInterp2D(interp) { 
  int i,j;
  for (i=0; i<16; i++)
    for (j=0; j<16; j++)
      wt[i][j] = interp.wt[i][j];
  if (interp.x1terp != NULL)
    x1terp = PolyInterp1D::Copy(interp.x1terp);
  else
    x1terp = NULL;
  if (interp.x2terp != NULL)
    x2terp = PolyInterp1D::Copy(interp.x2terp);
  else
    x2terp = NULL;
  if (interp.y1d != NULL) {
    y1d = new double ** [nnf];
    for (i=0; i<nnf; i++) {
      y1d[i] = new double * [nnx1];
      for(j=0; j<nnx1; j++) {
	y1d[i][j] = new double [nnx2];
	memcpy(y1d[i][j],interp.y1d[i][j],nnx2*sizeof(double));
      }
    }
  }
  else
    y1d = NULL;
  if (interp.y2d != NULL) {
    y2d = new double ** [nnf];
    for (i=0; i<nnf; i++) {
      y2d[i] = new double * [nnx1];
      for(j=0; j<nnx1; j++) {
	y2d[i][j] = new double [nnx2];
	memcpy(y2d[i][j],interp.y2d[i][j],nnx2*sizeof(double));
      }
    }
  }
  else
    y2d = NULL;
  if (interp.y12d != NULL) {
    y12d = new double ** [nnf];
    for (i=0; i<nnf; i++) {
      y12d[i] = new double * [nnx1];
      for(j=0; j<nnx1; j++) {
	y12d[i][j] = new double [nnx2];
	memcpy(y12d[i][j],interp.y12d[i][j],nnx2*sizeof(double));
      }
    }
  }
  else
    y12d = NULL;
  
  if (interp.c != NULL) {
    c = new double ** [nnf];
    for (i=0; i<nnf; i++) {
      c[i] = new double * [4];
      for(j=0; j<4; j++) {
	c[i][j] = new double [4];
	memcpy(c[i][j],interp.c[i][j],4*sizeof(double));
      }
    }
  }
  else
    c = NULL;

}


const char * SplineInterp2D::GetClassName() const {
  return "SplineInterp2D";
}

void SplineInterp2D::Destroy() {
  delete this;
}

SplineInterp2D::~SplineInterp2D() {
  if (x1terp != NULL)
    x1terp->Destroy();
  if (x2terp != NULL)
    x2terp->Destroy();
  if (y1d != NULL) {
    for (int i=0; i<nnf; i++) {
      for (int j=0; j<nnx1; j++) {
        delete [] y1d[i][j];
      }
      delete [] y1d[i];
    }
    delete [] y1d;
  }
  if (y2d != NULL) {
    for (int i=0; i<nnf; i++) {
      for (int j=0; j<nnx1; j++) {
        delete [] y2d[i][j];
      }
      delete [] y2d[i];
    }
    delete [] y2d;
  }
  if (y12d != NULL) {
    for (int i=0; i<nnf; i++) {
      for (int j=0; j<nnx1; j++) {
        delete [] y12d[i][j];
      }
      delete [] y12d[i];
    }
    delete [] y12d;
  }
  
  if (c != NULL) {
    for (int i=0; i<nnf; i++) {
      for (int j=0; j<4; j++) 
        delete [] c[i][j];
      delete [] c[i];
    }
    delete [] c;
  }

}

SplineInterp2D * SplineInterp2D::Copy(SplineInterp2D *interp) {
  return new SplineInterp2D(*interp);
}

SplineInterp2D * SplineInterp2D::Clone() const {
  return new SplineInterp2D(*this);
}

SplineInterp2D * SplineInterp2D::Create() {
  return new SplineInterp2D();
}

int SplineInterp2D::Init() {

  int i, j, k;
  
  if (y1d != NULL) {
    for (int i=0; i<nnf; i++) {
      for (int j=0; j<nnx1; j++) {
        delete [] y1d[i][j];
      }
      delete [] y1d[i];
    }
    delete [] y1d;
  }
  if (y2d != NULL) {
    for (int i=0; i<nnf; i++) {
      for (int j=0; j<nnx1; j++) {
        delete [] y2d[i][j];
      }
      delete [] y2d[i];
    }
    delete [] y2d;
  }
  if (y12d != NULL) {
    for (int i=0; i<nnf; i++) {
      for (int j=0; j<nnx1; j++) {
        delete [] y12d[i][j];
      }
      delete [] y12d[i];
    }
    delete [] y12d;
  }
  
  y1d = new double ** [nnf];
  y2d = new double ** [nnf];
  y12d = new double ** [nnf];

  for (i=0; i<nnf; i++) {
    y1d[i] = new double * [nnx1];
    y2d[i] = new double * [nnx1];
    y12d[i] = new double * [nnx1];
    y1d[i][0] = new double [nnx2];
    y2d[i][0] = new double [nnx2];
    y12d[i][0] = new double [nnx2];
    y1d[i][nnx1-1] = new double [nnx2];
    y2d[i][nnx1-1] = new double [nnx2];
    y12d[i][nnx1-1] = new double [nnx2];
    for(j=1; j<nnx1-1; j++) {
      y1d[i][j] = new double [nnx2];
      y2d[i][j] = new double [nnx2];
      y12d[i][j] = new double [nnx2];
      for(k=1; k<nnx2-1; k++) {
      	y1d[i][j][k] = (yy[i][(j+1)+nnx1*k]-yy[i][(j-1)+nnx1*k])/(xx1[j+1]-xx1[j-1]); 
	y2d[i][j][k] = (yy[i][j+nnx1*(k+1)]-yy[i][j+nnx1*(k-1)])/(xx2[k+1]-xx2[k-1]); 
	y12d[i][j][k] = (yy[i][(j+1)+nnx1*(k+1)]-yy[i][(j+1)+nnx1*(k-1)]-yy[i][(j-1)+nnx1*(k+1)]+yy[i][(j-1)+nnx1*(k-1)]) / ((xx1[j+1]-xx1[j-1])*(xx2[k+1]-xx2[k-1])); 
      }
    }
  }
  
  if (x1terp != NULL)
    x1terp->Destroy();
  if (x2terp != NULL)
    x2terp->Destroy();
  
  x1terp = PolyInterp1D::Create();
  x1terp->SetInterpolationPoints(xx1,&xx1,nnx1,1);
  x1terp->SetInterpolationOrder(2);
  x2terp = PolyInterp1D::Create();
  x2terp->SetInterpolationPoints(xx2,&xx2,nnx2,1);
  x2terp->SetInterpolationOrder(2);
  
  if (c != NULL) {
    for (i=0; i<nnf; i++) {
      for (j=0; j<4; j++) 
        delete [] c[i][j];
      delete [] c[i];
    }
    delete [] c;
  }
  
  c = new double ** [nnf];
  for (i=0; i<nnf; i++) {
    c[i] = new double * [4];
    for(j=0; j<4; j++) 
      c[i][j] = new double [4];
  }
  
  return 0;
}

int SplineInterp2D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    std::cerr << "SplineInterp2D::Evaluate() - Interpolation points not set\n";
    return -1;
  }

  if (y1d == NULL) {
    std::cerr << "SplineInterp2D::Evaluate() - Spline not initialized. Call method Init()\n";
    return -1;
  }
  
  int i,j,k;
  int idx1, idx2, l;
  double x1l, x1u, x2l, x2u;
  double **ya, **y1a, **y2a, **y12a;
  double t, u, d1, d2; 
  double *xxa, d1d2; 
  double **cl, **xa; 

  // Find the grid square. 
  idx1 = x1terp->nextHunt() ? x1terp->Hunt(x[0]) : x1terp->Locate(x[0]);
  idx2 = x2terp->nextHunt() ? x2terp->Hunt(x[1]) : x2terp->Locate(x[1]);

  x1l = xx1[idx1];
  x1u = xx1[idx1+1];
  x2l = xx2[idx2];
  x2u = xx2[idx2+1];
  
  if ((x1u == x1l) || (x2u == x2l)) {
    std::cerr << "Bad input to SplineInterp2D::Evaluate()\n";
    return -1;
  }
  ya = new double * [nnf];
  y1a = new double * [nnf];
  y2a = new double * [nnf];
  y12a = new double * [nnf];
  
  for (i=0; i<nnf; i++) {
    ya[i] = new double [4];
    y1a[i] = new double [4];
    y2a[i] = new double [4];
    y12a[i] = new double [4];
    ya[i][0] = yy[i][idx1+nnx1*idx2];
    ya[i][1] = yy[i][(idx1+1)+nnx1*idx2];
    ya[i][2] = yy[i][(idx1+1)+nnx1*(idx2+1)];
    ya[i][3] = yy[i][idx1+nnx1*(idx2+1)];
    y1a[i][0] = y1d[i][idx1][idx2];
    y1a[i][1] = y1d[i][idx1+1][idx2];
    y1a[i][2] = y1d[i][idx1+1][idx2+1];
    y1a[i][3] = y1d[i][idx1][idx2+1];
    y2a[i][0] = y2d[i][idx1][idx2];
    y2a[i][1] = y2d[i][idx1+1][idx2];
    y2a[i][2] = y2d[i][idx1+1][idx2+1];
    y2a[i][3] = y2d[i][idx1][idx2+1];
    y12a[i][0] = y12d[i][idx1][idx2];
    y12a[i][1] = y12d[i][idx1+1][idx2];
    y12a[i][2] = y12d[i][idx1+1][idx2+1];
    y12a[i][3] = y12d[i][idx1][idx2+1];
  
  }
  
  d1 = x1u - x1l;
  d2 = x2u - x2l;

  xxa = new double [nnf];

  cl = new double * [nnf];
  xa = new double * [nnf];
  for (i=0; i<nnf; i++) {
    cl[i] = new double [16];
    xa[i] = new double [16];
  }

  d1d2 = d1*d2;

  for(i=0; i<4; i++) { // Pack a temporary vector x.
    for (j=0; j<nnf; j++) {
      xa[j][i] = ya[j][i]; 
      xa[j][i+4] = y1a[j][i]*d1; 
      xa[j][i+8] = y2a[j][i]*d2; 
      xa[j][i+12] = y12a[j][i]*d1d2;
    }
  } 
  for(i=0; i<16; i++) { // Matrix-multiply by the stored table. 
    for (j=0; j<nnf; j++) {
      xxa[j] = 0.0; 
      for(k=0; k<16; k++) 
	xxa[j] += wt[i][k]*xa[j][k]; 
      cl[j][i] = xxa[j]; 
    }
  } 
  l = 0; 
  for (i=0; i<4; i++) { // Unpack the result into the output table. 
    for (j=0;j<4;j++) {
      for (k=0; k<nnf; k++)
	c[k][i][j]=cl[k][l];
      l++;
    }
  }
  for (i=0; i<nnf; i++) {
    delete [] cl[i];
    delete [] xa[i];
  }
  delete [] cl;
  delete [] xa;

  for (i=0; i<nnf; i++) {
    delete [] ya[i];
    delete [] y1a[i];
    delete [] y2a[i];
    delete [] y12a[i];
  }
  delete [] ya;
  delete [] y1a;
  delete [] y2a;
  delete [] y12a;


  t = (x[0]-x1l) / d1; // Equation (3.6.4) in NR. 
  u = (x[1]-x2l) / d2; 
  for (j=0; j<nnf; j++)
    y[j] = 0.0; 
  for (i=3;i>=0;i--) { // Equation (3.6.6) in NR. 
    for (j=0; j<nnf; j++)
      y[j]  = t*(y[j])+((c[j][i][3]*u+c[j][i][2])*u+c[j][i][1])*u+c[j][i][0]; 
  } 

  delete [] xxa;

  return 0;
}

int SplineInterp2D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  
  if (xx1 == NULL) {
    std::cerr << "SplineInterp2D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }

  double xr[2],xt[2],yc[nnf],yr[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xt[0] = x[0];
  xt[1] = x[1]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xt,yt);

  for (int i=0; i<nnf; i++) {
    y[i][0] = (yr[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][1] = (yt[i]-yc[i])/FINITE_DIFFERENCES_STEP;
  }
  return 0;
}

int SplineInterp2D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "SplineInterp2D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 2) {
    std::cerr << "SplineInterp2D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  double xr[2],xt[2],yc[nnf],yr[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xt[0] = x[0];
  xt[1] = x[1]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xt,yt);

  y[0] = (yr[0]-yc[0])/FINITE_DIFFERENCES_STEP + (yt[1]-yc[1])/FINITE_DIFFERENCES_STEP;

  return 0;
}

// SmoothingSplineInterp2D

SmoothingSplineInterp2D::SmoothingSplineInterp2D(): BaseInterp2D() {
  interpsx = NULL;
  interpy = NULL;
  x2terp = NULL;
  yloc = NULL;
  SS = 0;
  window = 2;
}

SmoothingSplineInterp2D::SmoothingSplineInterp2D(const SmoothingSplineInterp2D &interp):BaseInterp2D(interp) {
  if (interp.interpsx != NULL) {
    interpsx = new SmoothingSplineInterp1D * [nnx2];
    for (int i=0; i<nnx2; i++) 
      interpsx[i] = SmoothingSplineInterp1D::Copy(interp.interpsx[i]);
  }
  else
    interpsx = NULL;
  if (interp.interpy != NULL) {
    interpy = SmoothingSplineInterp1D::Copy(interp.interpy);
  }
  else
    interpy = NULL;
  if (interp.x2terp!=NULL)
    x2terp = PolyInterp1D::Copy(interp.x2terp);
  else
    x2terp = NULL;
  if (interp.yloc != NULL) {
    yloc = new double ** [nnx2];
    for (int i=0; i<nnx2; i++) {
      yloc[i] = new double * [nnf];
      for (int j=0; j<nnf; j++) {
        yloc[i][j] = new double [nnx1];
	memcpy(yloc[i][j],interp.yloc[i][j],nnx1*sizeof(double));
      }
    }
  }
  else
    yloc = NULL;
  SS = interp.SS;
  window = interp.window;
}

const char * SmoothingSplineInterp2D::GetClassName() const {
  return "SmoothingSplineInterp2D";
}

void SmoothingSplineInterp2D::Destroy() {
  delete this;
}

SmoothingSplineInterp2D::~SmoothingSplineInterp2D() {
  if (interpsx != NULL) {
    for (int i=0; i<nnx2; i++)
      interpsx[i]->Destroy();
    delete [] interpsx;
  }
  if (interpy != NULL)
    interpy->Destroy();
  if (x2terp != NULL)
    x2terp->Destroy();
  if (yloc != NULL){
    for (int i=0; i<nnx2; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
}

SmoothingSplineInterp2D * SmoothingSplineInterp2D::Copy(SmoothingSplineInterp2D *interp) {
  return new SmoothingSplineInterp2D(*interp);
}

SmoothingSplineInterp2D * SmoothingSplineInterp2D::Clone() const {
  return new SmoothingSplineInterp2D(*this);
}

SmoothingSplineInterp2D * SmoothingSplineInterp2D::Create() {
  return new SmoothingSplineInterp2D();
}

void SmoothingSplineInterp2D::SetSmoothingParameter(double S) {
  SS = S;
}

void SmoothingSplineInterp2D::SetWindow(int w) {
  window = w;
}

int SmoothingSplineInterp2D::Init() {
  if ((xx1 == NULL)||(xx2 == NULL) || (yy == NULL)) {
    std::cerr << "SmoothingSplineInterp2D::Init() - Interpolation points not set\n";
    return -1;
  }
  if ((window < 0) || (window > nnx2)) {
    std::cerr << "SmoothingSplineInterp2D::Init() - Windows size must be smaller than number of interpolation points along x axis\n";
    return -1;
  }
  int i,j;

  if (yloc != NULL){
    for (int i=0; i<nnx2; i++) {
      for (int j=0; j<nnf; j++)
        delete [] yloc[i][j];
      delete [] yloc[i];
    }
    delete [] yloc;
  }
  
  yloc = new double ** [nnx2];
  for (i=0; i<nnx2; i++) {
    yloc[i] = new double * [nnf];
    for (j=0; j<nnf; j++)
      yloc[i][j] = new double[nnx1];
  }

  if (interpsx != NULL) {
    for (int i=0; i<nnx2; i++)
      interpsx[i]->Destroy();
    delete [] interpsx;
  }

  interpsx = new SmoothingSplineInterp1D * [nnx2];
  for (i = 0; i<nnx2; i++) {
     for (j=0; j<nnf; j++) {
       memcpy(yloc[i][j],&yy[j][i*nnx1],nnx1*sizeof(double));
     }
     interpsx[i] = SmoothingSplineInterp1D::Create();
     interpsx[i]->SetInterpolationPoints(xx1,yloc[i],nnx1,nnf);
     interpsx[i]->SetSmoothingParameters(SS);
     interpsx[i]->Init();
  } 
  
  if (interpy != NULL)
    interpy->Destroy();
  interpy = SmoothingSplineInterp1D::Create();
  
  if (x2terp != NULL)
    x2terp->Destroy();
  x2terp = PolyInterp1D::Create();
  x2terp->SetInterpolationPoints(xx2,&xx2,nnx2,1);
  x2terp->SetInterpolationOrder(window);
  return 0;
}

int SmoothingSplineInterp2D::Evaluate(double *x, double *y) {

  if (xx1 == NULL) {
    std::cerr << "SmoothingSplineInterp2D::Evaluate() - Interpolation points not set\n";
    return -1;
  }
  
  int i,j;
  int idx2;
  double **yrow, *tmp;
  
  tmp = new double [nnf];
  yrow = new double * [nnf];
  for (i=0; i<nnf; i++)
    yrow[i] = new double[window];
  
  idx2 = x2terp->nextHunt() ? x2terp->Hunt(x[1]) : x2terp->Locate(x[1]);
  
  for(i=0; i<window; i++) {
    interpsx[idx2+i]->Evaluate(&x[0],tmp);
    for (j=0; j<nnf; j++)
      yrow[j][i] = tmp[j];
  }
  interpy->SetInterpolationPoints(&xx2[idx2],yrow,window,nnf);
  interpy->Init();
  interpy->Evaluate(&x[1],y);
  delete [] tmp;
  for (i=0; i<nnf; i++)
    delete [] yrow[i];
  delete [] yrow;
  return 0;
}

int SmoothingSplineInterp2D::EvaluateJacobian(double *x, double **y) {
  // Implemented with finite differences
  
  if (xx1 == NULL) {
    std::cerr << "SmoothingSplineInterp2D::EvaluateJacobian() - Interpolation points not set\n";
    return -1;
  }

  double xr[2],xt[2],yc[nnf],yr[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xt[0] = x[0];
  xt[1] = x[1]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xt,yt);

  for (int i=0; i<nnf; i++) {
    y[i][0] = (yr[i]-yc[i])/FINITE_DIFFERENCES_STEP;
    y[i][1] = (yt[i]-yc[i])/FINITE_DIFFERENCES_STEP;
  }
  return 0;
}

int SmoothingSplineInterp2D::EvaluateDivergence(double *x, double *y) {
  // Implemented with finite differences
  if (xx1 == NULL) {
    std::cerr << "SmoothingSplineInterp2D::EvaluateDivergence() - Interpolation points not set\n";
    return -1;
  }
  if (nnf != 2) {
    std::cerr << "SmoothingSplineInterp2D::EvaluateDivergence() - Invalid vector field, dimension of codomain must be 2.\n";
    return -1;
  }
  
  double xr[2],xt[2],yc[nnf],yr[nnf],yt[nnf];
  xr[0] = x[0]+FINITE_DIFFERENCES_STEP;
  xr[1] = x[1];
  xt[0] = x[0];
  xt[1] = x[1]+FINITE_DIFFERENCES_STEP;
  Evaluate(x,yc);
  Evaluate(xr,yr);
  Evaluate(xt,yt);

  y[0] = (yr[0]-yc[0])/FINITE_DIFFERENCES_STEP + (yt[1]-yc[1])/FINITE_DIFFERENCES_STEP;

  return 0;
}

} // namespace bal

