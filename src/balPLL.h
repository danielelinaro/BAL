/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balPLL.h
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
 * \file balPLL.h
 * \brief Definition of the class PLL
 */

#ifndef _BALPLL_
#define _BALPLL_

#include "balDynamicalSystem.h"
#include <cvode/cvode.h>

//#define WITHPHIERR
//#define FRACTIONAL
#define MISMATCH
//#define EXTEND

#ifdef FRACTIONAL
const int ndiv = 2;
#endif

namespace bal {

/**
 * \class PLL
 * \brief Implementation of a dynamical system that describes a
 * Phase-Locked Loop.\ This PLL is modelled as a switch system.
 *
 * \example pll.cpp
 * \sa DynamicalSystem
 */
class PLL : public DynamicalSystem {
 public:
  PLL();
  virtual ~PLL();
  
  int RHS (realtype t, N_Vector X, N_Vector Xdot, void *sys);
#ifdef CVODE25
  int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef CVODE26
  int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
		void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  int Events (realtype t, N_Vector X, realtype * event, void * data);
  void EventsConstraints (realtype t, N_Vector X, int * constraints, void * data);
  
  bool HasJacobian() const;
  bool HasEvents() const;
  bool HasEventsConstraints() const;
  
  void Reset();
  void ManageEvents(realtype t, N_Vector x, int * events, int * constraints = NULL);

  virtual DynamicalSystem* Clone() const;
  
  static const int npar;
  static const char * parname[14];
  
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
  realtype dt;		// duration of the falling edge of the output of the flip-flops
  realtype tau_d;       // delay of the output of the AND gate
  
  
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

} // namespace bal

#ifdef __cplusplus
extern "C" {
#endif

bal::DynamicalSystem* PLLFactory();
	
#ifdef __cplusplus
}
#endif

#endif

