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

/** 
 * \file balEye.h
 * \brief Definition of the class balEye
 */

#ifndef _BALEYE_
#define _BALEYE_

#include <cmath>
#include <cstdio>
#include <cstdlib>
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

/**
 * \class balEye
 * \brief Implementation of a dynamical system that represents the flow
 * of particles inside the human eye
 * 
 * \sa balDynamicalSystem
 */
class balEye : public balDynamicalSystem {
 public:
  virtual const char * GetClassName () const;
  static balEye * Create ();
  virtual balDynamicalSystem * Copy();
  virtual void Destroy ();
  
  int RHS (realtype t, N_Vector x, N_Vector xdot, void * data);
  int Events (realtype t, N_Vector x, realtype * event, void * data);
  
  bool HasEvents() const;
  bool SpecialOptions(void *opt);
  bool ReadVectorField(const char *filename);

 protected:
  balEye();
  balEye(const balEye& eye);
  virtual ~balEye();
  
 private:
  void DeleteVectorField();
  void AllocateVectorField();

  N_Vector xderiv;

  int nx, ny;
  double *x, *y, **u, **v;
  bool _dealloc;

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


