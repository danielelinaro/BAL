/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balHindmarshRose.h
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
 * \file balHindmarshRose.h
 * \brief Definition of the class balHindmarshRose
 */

#ifndef _BALHINDMARSHROSE_
#define _BALHINDMARSHROSE_

#include "balObject.h"
#include "balParameters.h"
#include "balDynamicalSystem.h"
#include <cvode/cvode.h>

#define XREST (-1.6)

/**
 * \class balHindmarshRose
 * \brief Implementation of a dynamical system that describes a the
 * Hindmarsh-Rose neuron model
 * 
 * \sa balDynamicalSystem
 */
class balHindmarshRose : public balDynamicalSystem {
 public:
  virtual const char * GetClassName () const;
  static balHindmarshRose * Create ();
  virtual balDynamicalSystem * Copy();
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
  
  //const double xrest;

 protected:
  balHindmarshRose();
  balHindmarshRose(const balHindmarshRose& hr);
  virtual ~balHindmarshRose();
  
 private:
  N_Vector xderiv;
};

#ifdef __cplusplus
extern "C" {
#endif

balDynamicalSystem* balHindmarshRoseFactory();
	
#ifdef __cplusplus
}
#endif

#endif
