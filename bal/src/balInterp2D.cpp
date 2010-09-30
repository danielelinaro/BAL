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
 *  \brief Implementation of the class balBilinearInterp2D 
 */

#include "balInterp2D.h"

balBilinearInterp2D::balBilinearInterp2D(double *x1v, double *x2v, double **yy, int mm, int nn)
  : m(mm), n(nn), y(yy) { 
  x1terp = balLinearInterp1D::Create(x1v,x1v,m);
  x2terp = balLinearInterp1D::Create(x2v,x2v,n);
}

balBilinearInterp2D::~balBilinearInterp2D() {
  x1terp->Destroy();
  x2terp->Destroy();
}

const char * balBilinearInterp2D::GetClassName() const {
  return "balBilinearInterp2D";
}

void balBilinearInterp2D::Destroy() {
  delete this;
}

balBilinearInterp2D * balBilinearInterp2D::Create(double *x1v, double *x2v, double **yy, int mm, int nn) {
  return new balBilinearInterp2D(x1v,x2v,yy,mm,nn);
}

double balBilinearInterp2D::interp(double x1p, double x2p) {
  int i, j;
  double yy, t, u;

  i = x1terp->cor ? x1terp->hunt(x1p) : x1terp->locate(x1p);
  j = x2terp->cor ? x2terp->hunt(x2p) : x2terp->locate(x2p);

  t = (x1p-x1terp->xx[i])/(x1terp->xx[i+1]-x1terp->xx[i]); 
  u = (x2p-x2terp->xx[j])/(x2terp->xx[j+1]-x2terp->xx[j]); 
  yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j]
    + (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
  return yy;
}
