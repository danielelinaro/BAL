/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp2D.h
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
 *  \file balInterp2D.h
 *  \brief Classes for two dimensional interpolation.
 *  
 */

#ifndef _BALINTERP2D_
#define _BALINTERP2D_

#include <cmath>
#include <cstdio>
#include "balObject.h"
#include "balInterp1D.h"
using namespace std;

class balBilinearInterp2D : public balObject {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const { return "balBilinearInterp2D" ; }
  /** Destroys a balBilinInterp2D. */
  virtual void Destroy() { delete this; }
  /** Creates a balBilinInterp2D. */
  static balBilinearInterp2D * Create(double *x1v, double *x2v, double **yy, int mm, int nn) {
    return new balBilinearInterp2D(x1v,x2v,yy,mm,nn);
  }

  double interp(double x1p, double x2p);
  
 protected:
   balBilinearInterp2D(double *x1v, double *x2v, double **yy, int mm, int nn)
     : m(mm), n(nn), y(yy) { 
     x1terp = balLinearInterp1D::Create(x1v,x1v,m);
     x2terp = balLinearInterp1D::Create(x2v,x2v,n);
   }

   ~balBilinearInterp2D() {
     x1terp->Destroy();
     x2terp->Destroy();
   }
  
 private:
  int m, n;
  double **y;
  balLinearInterp1D *x1terp, *x2terp;
};


#endif


