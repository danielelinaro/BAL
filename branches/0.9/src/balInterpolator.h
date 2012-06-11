/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterpolator.h
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
 *  \file balInterpolator.h
 *  \brief Definition of abstract class Interpolator. 
 */

#ifndef _BALINTERPOLATOR_
#define _BALINTERPOLATOR_

#include <cmath>
#include <cstdio>
#include "balObject.h"

namespace bal {

/**
 * \class Interpolator 
 * \brief Base class for generic interpolator (also smooth approximator) of functions \f$f: R^{n} \rightarrow R^{m}\f$.
 *
 * \sa Interp1D, Interp2D, Interp3D
 */

class Interpolator : public Object {
  public:
  
    /** Returns the name of the class. */
    virtual const char * getClassName() const;
    
    /** Destroys a Interpolator. */
    virtual void Destroy();
    
    virtual Interpolator * Clone() const = 0;
    
  /** Evaluates the function in a generic point.
	 * \param x Point of \f$R^{n}\f$.
	 * \param y Value of \f$f(x)\f$. 
	 * \retval -1 An error has occurred.
	 * \retval 0 Default value.
	 */ 
    virtual int Evaluate(double *x,double *y) = 0;
    
	/** Evaluates the Jacobian matrix of the function in a generic point.
	 * \param x Point of \f$R^{n}\f$.
	 * \param y \f$m \times n\f$ matrix of evaluated \f$J(x)\f$ value.
	 * \retval -1 An error has occurred.
	 * \retval 0 Default value.
	 */ 
	virtual int EvaluateJacobian(double *x, double **y) = 0;
  
  /** Evaluates the divergence of the vector field in a generic point (if \f$n = m\f$).
	 * \param x Point of \f$R^{n}\f$.
	 * \param y Value of \f$div(x)\f$.
	 * \retval -1 An error has occurred.
	 * \retval 0 Default value.
	 */
    virtual int EvaluateDivergence(double *x, double *y) = 0;  
    
	/** Initializes interpolator, the user has to perform this operation before evaluating the function. 
	 * \retval -1 An error has occurred.
	 * \retval 0 Default value.
	 */
    virtual int Init();

    /** Gets the dimension of the function domain (\f$n\f$). */
    int GetDomainDimensions();
    
    /** Gets the dimension of the function codomain (\f$m\f$). */
    int GetCodomainDimensions();

  protected:
    /** Constructor of Interpolator. */
    Interpolator();
    
    /** Destructor of Interpolator. */
    virtual ~Interpolator();

    /** Copy constructor of Interpolator. */
    Interpolator(const Interpolator & interp);

    int nnd, nnf;

};

} // namespace bal

#endif

