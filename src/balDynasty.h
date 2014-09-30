/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balDynasty.h
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
 * \file balDynasty.h
 * \brief Definition of the class Dynasty
 */

#ifndef _BALDYNASTY_
#define _BALDYNASTY_

#include "balObject.h"
#include "balParameters.h"
#include "balDynamicalSystem.h"
#include <cstring>
#include <cvode/cvode.h>

#define EPS (1.0E-4)

namespace bal {

enum {MINIMA, MAXIMA, ANY};

/**
 * \class Dynasty
 * \brief Implementation of a dynamical system that describes a
 * model of Chinese dynasty.
 * 
 * \example dynasty.cpp
 *
 * \sa DynamicalSystem
 */
class Dynasty : public DynamicalSystem {
 public:
  Dynasty();
  Dynasty(const Dynasty& hr);
  virtual ~Dynasty();
  
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
  bool SpecialOptions(const void *opt);

  virtual DynamicalSystem* Clone() const;

 private:
  N_Vector xderiv;
  int constraint_type;
};

} // namespace bal

#ifdef __cplusplus
extern "C" {
#endif

bal::DynamicalSystem* DynastyFactory();
	
#ifdef __cplusplus
}
#endif

#endif

