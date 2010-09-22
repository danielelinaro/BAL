/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balPLL.h
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
 * \file balPLL.h
 * \brief Definition of the class balPLL
 */

#ifndef _BALPLL_
#define _BALPLL_

#include "balObject.h"
#include "balParameters.h"
#include "balDynamicalSystem.h"
#include <cvode/cvode.h>

//#define WITHPHIERR
#define FRACTIONAL
//#define MISMATCH
//#define EXTEND

#ifdef FRACTIONAL
const int ndiv = 2;
#endif

/**
 * \class balPLL
 * \brief Implementation of a dynamical system that describes a
 * Phase-Locked Loop
 *
 * This PLL is modelled as a switch system.
 *
 * \sa balDynamicalSystem
 */
class balPLL : public balDynamicalSystem {
 public:
  virtual const char * GetClassName () const;
  static balPLL * Create ();
  virtual void Destroy ();
  
  int RHS (realtype t, N_Vector X, N_Vector Xdot, void * data);
#ifdef CVODE25
  int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef CVODE26
  int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
		void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  int Events (realtype t, N_Vector X, realtype * event, void * data);
  void EventsConstraints (realtype t, N_Vector X, int * constraints, void * data);
  
  bool HasJacobian() const;
  bool HasEvents() const;
  bool HasEventsConstraints() const;
  
  void Reset();
  void ManageEvents(realtype t, N_Vector x, int * events, int * constraints = NULL);
  
  static const int npar;
  static const char * parname[14];
  
 protected:
  balPLL();
  virtual ~balPLL();
  
 private:
  const realtype pi;
  
  int n;
#ifndef WITHPHIERR
  bool zu, zd;
#else
  bool S;
#endif
  
#ifndef WITHPHIERR
  // variable to be reset at every new integration
  bool wait_reset;
#endif
  
  /*************************************************/
  /****** CIRCUIT PARAMETERS WITHOUT THE VCO *******/
  /*************************************************/
  realtype C0;			// the state variable w is associated to C0
  realtype C1, R1;	// make up the filter, and the state variable r is associated to C1
  realtype Vdd;			// supply voltage
  realtype Aud;			// current in the charge-pump
#ifdef MISMATCH
  realtype Aud_mismatch;
#endif
  realtype fREF, T;	// frequency (and corresponding period) of the square wave input to one of the flip-flops
  realtype treset;	// time at which both Zu and Zd are equal to 1
  realtype dt;			// interval after which the circuit is reset (i.e., treset+dt)
  
  
  /*****************************/
  /******* VCO PARAMETERS ******/
  /*****************************/
  
  /***** VCO modelled as a polar oscillator *****/
  // variable parameters
  realtype omega0, rhoap, rho0;
  realtype k0, Krho, KVCOa, KVCOb, KVCOc, Kap, alpha;
  realtype tuning_coeff;
  
  // coefficient of division of the frequency divider
#ifndef FRACTIONAL
  int N;
#else
  int N[ndiv], idx;
#endif

#ifndef WITHPHIERR
  // counter
  int cnt;
  // output of the frequency divider
  int divout;
#endif

};

#ifdef __cplusplus
extern "C" {
#endif

balDynamicalSystem* balPLLFactory();
	
#ifdef __cplusplus
}
#endif


#endif

