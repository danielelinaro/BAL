/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balRossler.h
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
 * \file balRossler.h
 * \brief Definition of the example class Rossler 
 */

#ifndef _BALROSSLER_
#define _BALROSSLER_

#include "balDynamicalSystem.h"
#include <cvode/cvode.h>

namespace bal {

/**
 * \class Rossler
 * \brief This class defines the methods needed to characterize and use Rossler dynamical system with BAL.
 * \sa DynamicalSystem
 */
class Rossler : public DynamicalSystem {
 public:
  Rossler();
  Rossler(const Rossler& hr);
  virtual ~Rossler();
  
  int RHS (realtype t, N_Vector x, N_Vector xdot, void *sys);
#ifdef CVODE25
  int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef CVODE26
  int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
		void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
/** Events functions are system equations themselves.*/
  int Events (realtype t, N_Vector x, realtype * event, void * data);
/** Allows to locate a maximum of the \f$i_{th}\f$ component. */
  void EventsConstraints (realtype t, N_Vector x, int * constraints, void * data);

  bool HasJacobian() const;
  bool HasEvents() const;
  bool HasEventsConstraints() const;
  
  virtual DynamicalSystem* Clone() const;

 private:
  N_Vector xderiv;
};

} // namespace bal

#ifdef __cplusplus
extern "C" {
#endif
	
bal::DynamicalSystem* RosslerFactory();
  
#ifdef __cplusplus
}
#endif

#endif

