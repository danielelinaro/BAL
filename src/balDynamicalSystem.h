/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    DynamicalSystem.h
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
 * \brief Definition of the class DynamicalSystem
 */

#ifndef _BALDYNAMICALSYSTEM_
#define _BALDYNAMICALSYSTEM_

#include <string>
#include <boost/shared_ptr.hpp>

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

namespace bal {

/** 
 * \class DynamicalSystem
 * \brief Base class to represent dynamical systems.
 * 
 * DynamicalSystem is the base class for all dynamical systems that are
 * to be integrated using the Bifurcation Analysis Library. Every class
 * that inherits from DynamicalSystem must have at least the RHS
 * method (implementing the vector field of a dynamical system) and
 * possibly a method that implements the Jacobian matrix of the system,
 * even though this is not strictly necessary for the integration of the
 * model.
 *
 * \example dynsys.cpp
 *
 * \sa Parameters
 */
class DynamicalSystem : public Object {
public:
  DynamicalSystem();
  DynamicalSystem(const DynamicalSystem& system);
  virtual ~DynamicalSystem();

  /** Definition of system's right hand side equations.
      The prototype of this method is defined as required by CVODE integration library.
      For more informations refer to the
      <a href="https://computation.llnl.gov/casc/sundials/documentation/cv_guide/node5.html#SECTION00561000000000000000">
      CVODE documentation</a>.

      \param t Time instant in which the system has to be evaluated (only for time-variant systems)
      \param x Current state vector
      \param xdot Vector in which is stored the resulting evolution of the system
      \param sys Pointer to the dynamical system
  */
  virtual int RHS (realtype t, N_Vector x, N_Vector xdot, void *sys) = 0;
  static int RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void *sys);

  
#ifdef CVODE25
  static int JacobianWrapper (long int N, DenseMat J, realtype t, N_Vector x, 
			      N_Vector fy, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
					DenseMat J, void *sys);
  virtual int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, 
			N_Vector fy, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  
#ifdef CVODE26
  /** Jacobian static wrapper for CVODE C static library. */
  static int JacobianWrapper (int N, realtype t, N_Vector x, N_Vector fy, 
			      DlsMat J, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int JacobianFiniteDifferences (long int N, realtype t, N_Vector x,
					DlsMat J, void *sys);
					
  /** Definition of system's Jacobian matrix calculated analytically.
      The prototype of this method is defined as required by CVODE integration library.
      The other arguments are intentionally left blank.
      For more informations refer to <a href="https://computation.llnl.gov/casc/sundials/documentation/cv_guide/node5.html#SECTION00565000000000000000">CVODE documentation</a>.
      \param t Time instant in which the jacobian matrix has to be evaluated (only for time-variant systems).
      \param x Current state vector.
      \param J Result of jacobian matrix evaluation.
      \param jac_data Vector of system's parameters values.
  */
  virtual int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, 
			DlsMat J, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  
  /** Events static wrapper for CVODE C static library. */
  static int EventsWrapper (realtype t, N_Vector x, realtype * event, void * sys);
  
  /** With this method it is possible to define the event functions (i.e.\ implicit equations
   *  of Poincare' section) which will be evaluated at every integration step by CVODE
   *  rootfinding algorithm.
   *	
   *  The prototype of this method is defined as required by CVODE integration library.
   *  For more informations refer to the 
   *  <a href="https://computation.llnl.gov/casc/sundials/documentation/cv_guide/node5.html#ss:rootFn">
   *  CVODE documentation</a>.
   */
  virtual int Events (realtype t, N_Vector x, realtype * event, void * data);
  
  /** Defines additional constraints that will be evaluated when a root of an event functions is found.
   *	\remark The implementation of this method is foundamental for the correct detection of limit cycle performed
   *					in \ref bal::EVENTS and \ref bal::BOTH \ref integration_mode.
   *					So it is possible to univocally detect limit cycles by their intersection with 
   *					the defined Poincare' surface using \ref Events and \ref EventsConstraints methods.
   *  \param t Current time.
   *  \param x Current state vector.
   *	\param constraints Vector of dimension \ref nev in which \f$constraints[i]=1\f$ if the \f$i_{th}\f$ constraint is active otherwise is \f$0\f$.
   *  \param data Vector of the system's parameters.
   */
  virtual void EventsConstraints (realtype t, N_Vector x, int * constraints, void * data);
  
  /** Allows to define more complex dynamics of a system, managing occuring event. For istance
   *	a change of parameters or even a switch of the entire vector field (switching systems).
   *
   * \sa PLL
   */
  virtual void ManageEvents(realtype t, N_Vector X, int * events, int * constraints = NULL);
  
  virtual bool HasJacobian() const;
  virtual bool HasEvents() const;
  virtual bool HasEventsConstraints() const;
  
  virtual void Reset();
  virtual bool SpecialOptions(const void *opt);
  
  int GetNumberOfEvents() const;
  int GetDimension() const;
  int GetOriginalDimension() const;
  int GetNumberOfParameters() const;
  // makes a copy of p
  void SetParameters(const Parameters& params);
  // just storestruct p
  void SetParameters(boost::shared_ptr<Parameters>& params);
  boost::shared_ptr<Parameters> GetParameters() const;

  /** Extends dynamical system dimensionality to calculate Lyapunov exponents.
   *  The algorithm used is described in Alan Wolf et al.
   *  "Determining Lyapunov exponents from a time series", Physica D: Nonlinear Phenomena, 16(3):285 – 317, 1985.
   */
  void Extend(bool extend);
  /** Tells if the system is extended. */
  bool IsExtended() const;

  virtual DynamicalSystem* Clone() const = 0;

protected:
  void SetDimension(int n_);
  void SetNumberOfParameters(int p_);
  void SetNumberOfEvents(int nev_);
  
private:
  int n;
  int p;
  int nev;
  
  int nExt;
  bool ext;
  
  bool dealloc_;
#ifdef CVODE25
  DenseMat jac;
#endif
#ifdef CVODE26
  DlsMat jac;
#endif

  boost::shared_ptr<Parameters> pars;
};

} // namespace bal

#endif

