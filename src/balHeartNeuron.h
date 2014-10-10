/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balHeartNeuron.h
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
 * \file balHeartNeuron.h
 * \brief Definition of the class HeartNeuron
 */

#ifndef _BALHEARTNEURON_
#define _BALHEARTNEURON_

#include "balDynamicalSystem.h"
#include <cvode/cvode.h>

namespace bal {

/**
 * \class HeartNeuron
 * \brief Implementation of a dynamical system that describes a heart
 * neuron model
 * 
 * \sa DynamicalSystem HindmarshRose
 */
class HeartNeuron : public DynamicalSystem {
 public:
  HeartNeuron();
  HeartNeuron(const HeartNeuron& hn);
  virtual ~HeartNeuron();
  
  int RHS (realtype t, N_Vector x, N_Vector xdot, void *sys);
#ifdef CVODE25
  int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef CVODE26
  int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
		void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  int Events (realtype t, N_Vector x, realtype * event, void * data);
  void EventsConstraints (realtype t, N_Vector x, int * constraints, void * data);
  
  bool HasJacobian() const;
  bool HasEvents() const;
  bool HasEventsConstraints() const;
  
  virtual DynamicalSystem* Clone() const;

  inline static realtype BoltzmannF(realtype a, realtype b, realtype V) {
    return 1./(1. + exp(a*(b + V)));
  }
  
  inline static realtype BoltzmannDFDV(realtype a, realtype b, realtype V) {
    realtype e = exp(a*(b+V));
    return - (a*e) / ((1.+e)*(1.+e));			
  }
  
 private:
  N_Vector xderiv;
  const realtype C,gK2,EK,ENa,gNa,E1,g1,tauNa;	// Constant parameters
  realtype A[3], B[3];
};


} // namespace bal

#ifdef __cplusplus
extern "C" {
#endif

bal::DynamicalSystem* HeartNeuronFactory();
	
#ifdef __cplusplus
}
#endif

#endif
