/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterpSystem.h
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
 * \file balInterpSystem.h
 * \brief Definition of the class balInterpSystem
 */

#ifndef _BALINTERPSYSTEM_
#define _BALINTERPSYSTEM_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

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
#include "balInterpolator.h"

/**
 * \class balInterpSystem
 * \brief Implementation of a dynamical system whos vector vield is known only n a regular grid
 * 
 * \sa balDynamicalSystem
 */
class balInterpSystem : public balDynamicalSystem {
 public:
  virtual const char * GetClassName () const;
  static balInterpSystem * Create ();
  virtual balDynamicalSystem * Copy();
  virtual void Destroy ();
  
  int RHS (realtype t, N_Vector x, N_Vector xdot, void * data);
  int Events (realtype t, N_Vector x, realtype * event, void * data);
  
  bool HasEvents() const;
  int SetInterpolator(balInterpolator *interp);

 protected:
  balInterpSystem();
  balInterpSystem(const balInterpSystem& interpsystem);
  virtual ~balInterpSystem();
  
 private:

  N_Vector xderiv;
  balInterpolator *interpolator;
};

#ifdef __cplusplus
extern "C" {
#endif

balDynamicalSystem* balInterpSystemFactory();
	
#ifdef __cplusplus
}
#endif

#endif


