/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balEye.h
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

// .NAME balEye
// 
// .SECTION Description
//
// .SECTION See also
// balDynamicalSystem balInterp1D balInterp2D

#ifndef _BALEYE_
#define _BALEYE_

#include <cmath>
#include <cstdio>
#include <cstring>

#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#ifdef CVODE26
#include <sundials/sundials_direct.h>
#endif
#include <cvode/cvode.h>

#include "balCommon.h"
#include "balObject.h"
#include "balParameters.h"
#include "balDynamicalSystem.h"
#include "balInterp2D.h"

class balEye : public balDynamicalSystem {
 public:
  virtual const char * GetClassName () const { return "balEye"; }
  static balEye * Create () { return new balEye; }
  virtual balDynamicalSystem * Copy() { return new balEye(*this); }
  virtual void Destroy () { this->~balEye(); }
  
  int RHS (realtype t, N_Vector x, N_Vector xdot, void * data);
  int Events (realtype t, N_Vector x, realtype * event, void * data);
  
  bool HasEvents() const { return true; }
  bool SpecialOptions(void *opt);
  bool ReadVectorField(const char *filename);

 protected:
  balEye();
  balEye(const balEye& hr);
  virtual ~balEye();
  
 private:
  void DeleteVectorField();
  void AllocateVectorField();

  N_Vector xderiv;

  int nx, ny;
  double *x, *y, **u, **v;
  bool alloc_mem;

  balBilinearInterp2D *u_interp, *v_interp;
};

#ifdef __cplusplus
extern "C" {
#endif

balDynamicalSystem* balEyeFactory();
	
#ifdef __cplusplus
}
#endif

#endif


