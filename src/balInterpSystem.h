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
 * \brief Definition of the class InterpSystem
 */

#ifndef _BALINTERPSYSTEM_
#define _BALINTERPSYSTEM_

#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#ifdef CVODE26
#include <sundials/sundials_direct.h>
#endif
#include <cvode/cvode.h>

#include "balDynamicalSystem.h"
#include "balInterpolator.h"

namespace bal {

/**
 * \class InterpSystem
 * \brief Implementation of a dynamical system whos vector field is known only on a regular grid.
 * 
 * \sa DynamicalSystem
 */
class InterpSystem : public DynamicalSystem {
 public:
  InterpSystem();
  InterpSystem(const InterpSystem& interpsystem);
  virtual ~InterpSystem();

  int RHS (realtype t, N_Vector x, N_Vector xdot, void *sys);
  int Events (realtype t, N_Vector x, realtype * event, void * data);
  
  bool HasEvents() const;
  int SetInterpolator(Interpolator *interp);
  
  // Used to set the direction of the integration (forward or backward in time)
  bool SpecialOptions(const void *opt); 
  
  virtual DynamicalSystem* Clone() const;

 private:

  N_Vector xderiv;
  Interpolator *interpolator;
  bool backward, arclength;
  bool _dealloc;
};

} // namespace bal

#ifdef __cplusplus
extern "C" {
#endif

bal::DynamicalSystem* InterpSystemFactory();
	
#ifdef __cplusplus
}
#endif

#endif

