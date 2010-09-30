/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balDynamicalSystem.h
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
 * \file balDynamicalSystem.h
 * \brief Definition of the class balDynamicalSystem
 */

#ifndef _BALDYNAMICALSYSTEM_
#define _BALDYNAMICALSYSTEM_

#include <cmath>
#include <cstdio>

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

/** 
 * \class balDynamicalSystem
 * \brief Base class to represent dynamical systems.
 * 
 * balDynamicalSystem is the base class for all dynamical systems that are
 * to be integrated using the Bifurcation Analysis Library. Every class
 * that is inherited from balDynamicalSystem must have at least the RHS
 * method (implementing the vector field of a dynamical system) and
 * possibly a method that implements the Jacobian matrix of the system,
 * even though this is not strictly necessary for the integration of the
 * model.
 *
 * \sa balParameters
 */
class balDynamicalSystem : public balObject {
 public:
  virtual const char * GetClassName () const;
  static balDynamicalSystem * Create ();
  virtual balDynamicalSystem * Copy ();
  virtual void Destroy ();
  
  virtual int RHS (realtype t, N_Vector x, N_Vector xdot, void * data);
  static int RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void * sys);
  
#ifdef CVODE25
  static int JacobianWrapper (long int N, DenseMat J, realtype t, N_Vector x, 
			      N_Vector fy, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
					DenseMat J, void *sys);
  virtual int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, 
			N_Vector fy, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  
#ifdef CVODE26
  static int JacobianWrapper (int N, realtype t, N_Vector x, N_Vector fy, 
			      DlsMat J, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
					DlsMat J, void *sys);
  virtual int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, 
			DlsMat J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  
  static int EventsWrapper (realtype t, N_Vector x, realtype * event, void * sys);
  virtual	int Events (realtype t, N_Vector x, realtype * event, void * data);
  
  virtual void EventsConstraints (realtype t, N_Vector x, int * constraints, void * data);
  
  virtual bool HasJacobian() const;
  virtual bool HasEvents() const;
  virtual bool HasEventsConstraints() const;
  
  virtual void Reset();
  virtual bool SpecialOptions(const void *opt);
  virtual void ManageEvents(realtype t, N_Vector X, int * events, int * constraints = NULL);
  
  int GetNumberOfEvents() const;
  int GetDimension() const;
  int GetOriginalDimension() const;
  int GetNumberOfParameters() const;
  void SetParameters(balParameters *) throw(balException);
  balParameters * GetParameters() const;
  
  void Extend(bool extend);
  bool IsExtended() const;
  
 protected:
  balDynamicalSystem();
  balDynamicalSystem(const balDynamicalSystem& system);
  virtual ~balDynamicalSystem();
  void SetDimension(int n_);
  void SetNumberOfParameters(int p_);
  void SetNumberOfEvents(int nev_);
  
 private:
  int n;
  int p;
  int nev;
  
  int nExt;
  bool ext;
  
  bool _dealloc;

#ifdef CVODE25
  DenseMat jac;
#endif
#ifdef CVODE26
  DlsMat jac;
#endif
  balParameters * pars;
};

#endif
