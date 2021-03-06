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
 *  \file balInterp1D.h
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
 * \brief Base class for generic interpolator (also smooth approximator) of functions f: R^{nd} -> R^{nf}.
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
    
    /** Evaluates the function in point x. The result is stored in array y. 
     *  If an error has occurred the return value is -1, otherwise it is 0. */ 
    virtual int Evaluate(double *x,double *y) = 0;
    
    /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x nd). 
     *  If an error has occurred the return value is -1, otherwise it is 0. */ 
    virtual int EvaluateJacobian(double *x, double **y) = 0;
  
    /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
     *  If an error has occurred the return value is -1, otherwise it is 0. */
    virtual int EvaluateDivergence(double *x, double *y) = 0;  
    
    /** Initializes the Interpolator. You have to perform this operation before evaluating the function. 
     *  If an error has occurred the return value is -1, otherwise it is 0. */
    virtual int Init();

    /** Gets the dimension of the function domain (nd) */
    int GetDomainDimensions();
    
    /** Gets the dimension of the function codomain (nf) */
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

