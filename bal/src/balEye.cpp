/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balEye.cpp
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

/** 
 * \file balEye.cpp
 * \brief Implementation of the class balEye
 */

#include "balEye.h"

balDynamicalSystem* balEyeFactory() {
  return balEye::Create();
}

balEye::balEye() {
  SetDimension(2);
  SetNumberOfParameters(0);
  SetNumberOfEvents(GetDimension());
  xderiv = N_VNew_Serial(GetDimension());
  _dealloc = false;
  nx = ny = -1;
  x = y = NULL;
  u = v = NULL;
  u_interp = v_interp = NULL;
}

balEye::balEye(const balEye& eye) : balDynamicalSystem( eye ) {
  xderiv = N_VNew_Serial(eye.GetDimension());
  for(int i = 0; i < eye.GetDimension(); i++)
    Ith(xderiv,i)=Ith(eye.xderiv,i);
  _dealloc = eye._dealloc;
  if(_dealloc) {
    int i, j;
    nx = eye.nx;
    ny = eye.ny;
    x = new double[nx];
    y = new double[ny];
    u = new double*[nx];
    v = new double*[nx];
    for(i=0; i<nx; i++) {
      u[i] = new double[ny];
      v[i] = new double[ny];
      x[i] = eye.x[i];
      for(j=0; j<ny; j++) {
	u[i][j] = eye.u[i][j];
	v[i][j] = eye.v[i][j];
      }
    }
    for(j=0; j<ny; j++)
      y[j] = eye.y[j];
  }
}

balEye::~balEye() {
  N_VDestroy_Serial(xderiv);
  DeleteVectorField();
}

balEye * balEye::Create () {
  return new balEye;
}

void balEye::Destroy () {
  delete this;
}

balDynamicalSystem * balEye::Copy() {
  return new balEye(*this);
}
  
const char * balEye::GetClassName () const { 
  return "balEye";
}

bool balEye::HasEvents() const {
  return true;
}

bool balEye::SpecialOptions(const void *opt) {
  printf("balEye::SpecialOptions %s\n", (const char *) opt);
  return ReadVectorField((const char*) opt);
}

void balEye::DeleteVectorField() {
  if(_dealloc) {
    delete x;
    delete y;
    for(int i=0; i<nx; i++) {
      delete u[i];
      delete v[i];
    }
    delete u;
    delete v;
    u_interp->Destroy();
    v_interp->Destroy();
    _dealloc = false;
  }
}

void balEye::AllocateVectorField() {
  DeleteVectorField();
  x = new double[nx];
  y = new double[ny];
  u = new double*[nx];
  v = new double*[nx];
  for(int i=0; i<nx; i++) {
    u[i] = new double[ny];
    v[i] = new double[ny];
  }
  _dealloc = true;
}

bool balEye::ReadVectorField(const char *filename) {
  DeleteVectorField();

  char line[256], *token;
  int i, j;
  FILE *fid;
  
  printf("Reading from %s\n", filename);
  fid = fopen(filename,"r");
  if(fid == NULL)
    return false;

  if(fgets(line,256,fid) == NULL) {
    fclose(fid);
    return false;
  }

  if(fgets(line,256,fid) == NULL) {
    fclose(fid);
    return false;
  }
  token = strtok(line," ");
  token = strtok(NULL,", ");
  nx = atoi(token+2);
  token = strtok(NULL,", ");
  ny = atoi(token+2);
  printf("nx = %d ny = %d\n", nx, ny);

  AllocateVectorField();

  for(j=0; j<ny; j++) {
    for(i=0; i<nx; i++) {
      fscanf(fid, "%lf%lf%lf%lf", &x[i], &y[j], &u[i][j], &v[i][j]);
      fgets(line,256,fid);
    }
  }

  u_interp = balBilinearInterp2D::Create(x,y,u,nx,ny);
  v_interp = balBilinearInterp2D::Create(x,y,v,nx,ny);

  fclose(fid);
  return true;
}

int balEye::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) {
  Ith(xdot,0) = u_interp->interp(Ith(x,0),Ith(x,1));
  Ith(xdot,1) = v_interp->interp(Ith(x,0),Ith(x,1));
  return CV_SUCCESS;
}

int balEye::Events (realtype t, N_Vector x, realtype * event, void * data) {
  RHS(t,x,xderiv,data);
  for(int i=0; i<GetNumberOfEvents(); i++)
    event[i] = Ith(xderiv,i);
  return CV_SUCCESS;
}
