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
 * \brief Definition of the class Rossler
 */

#ifndef _BALROSSLER_
#define _BALROSSLER_

#include "balObject.h"
#include "balParameters.h"
#include "balDynamicalSystem.h"
#include <cvode/cvode.h>

namespace bal {

/**
 * \class Rossler
 * \brief Implementation of the Rossler system
 *
 * \sa DynamicalSystem
 */
class Rossler : public DynamicalSystem {
 public:
  virtual const char * GetClassName () const;
  static Rossler * Create ();
  static Rossler * Copy (Rossler *ros);
  virtual DynamicalSystem * Clone() const;
  virtual void Destroy ();
  
  int RHS (realtype t, N_Vector x, N_Vector xdot, void * data);
#ifdef CVODE25
  int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef CVODE26
  int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
		void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  int Events (realtype t, N_Vector x, realtype * event, void * data);
  void EventsConstraints (realtype t, N_Vector x, int * constraints, void * data);

  bool HasJacobian() const;
  bool HasEvents() const;
  bool HasEventsConstraints() const;
  
 protected:
  Rossler();
  Rossler(const Rossler& hr);
  virtual ~Rossler();
  
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

