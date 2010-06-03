/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balDynamicalSystem.h
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

// .NAME balDynamicalSystem - base class to represent dynamical systems.
// 
// .SECTION Description
// balDynamicalSystem is the base class for all dynamical systems that are
// to be integrated using the Bifurcation Analysis Library. Every class
// that is inherited from balDynamicalSystem must have at least the RHS
// method (implementing the vector field of a dynamical system) and
// possibly a method that implements the Jacobian matrix of the system,
// even though this is not strictly necessary for the integration of the
// model.
//
// .SECTION See also
// balParameters

#ifndef _BALDYNAMICALSYSTEM_
#define _BALDYNAMICALSYSTEM_

#include <cmath>
#include <cstdio>

#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <cvode/cvode.h>

#include "balCommon.h"
#include "balObject.h"
#include "balParameters.h"

class balDynamicalSystem : public balObject {
	public:
		virtual const char * GetClassName () const { return "balDynamicalSystem"; }
		static balDynamicalSystem * Create () { return new balDynamicalSystem; }
		virtual balDynamicalSystem * Copy () { return new	balDynamicalSystem(*this); }
		virtual void Destroy () { this->~balDynamicalSystem(); }

		virtual int RHS (realtype t, N_Vector x, N_Vector xdot, void * data);
		static int RHSWrapper (realtype t, N_Vector x, N_Vector xdot, void * sys);

#ifdef CVODE25
		static int JacobianWrapper (long int N, DenseMat J, realtype t, N_Vector x, 
													N_Vector fy, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
		virtual int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, 
													N_Vector fy, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

#ifdef CVODE26
		static int JacobianWrapper (int N, realtype t, N_Vector x, N_Vector fy, 
													DlsMat J, void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
		virtual int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, 
													DlsMat J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

		static int EventsWrapper (realtype t, N_Vector x, realtype * event, void * sys);
		virtual	int Events (realtype t, N_Vector x, realtype * event, void * data);

		virtual void EventsConstraints (realtype t, N_Vector x, int * constraints, void * data);

		virtual bool HasJacobian() const { return false; }
		virtual bool HasEvents() const { return false; }
		virtual bool HasEventsConstraints() const { return false; }

		virtual void Reset() {}
		virtual void ManageEvents(realtype t, N_Vector X, int * events, int * constraints = NULL) {}

		int GetNumberOfEvents() const { return nev; }
		int GetDimension() const { return n; }
		int GetNumberOfParameters() const { return p; }
		void SetParameters(balParameters *) throw(balException);
		balParameters * GetParameters() const;

	protected:
		balDynamicalSystem();
		balDynamicalSystem(const balDynamicalSystem& system);
		virtual ~balDynamicalSystem();
		void SetDimension(int n_) { if(n_ > 0) n = n_; }
		void SetNumberOfParameters(int p_) { if(p_ >= 0) p = p_; }
		void SetNumberOfEvents(int nev_) { if(nev_ >= 0) nev = nev_; }

	private:
		int n;
		int p;
		int nev;
		balParameters * pars;
};

#endif


