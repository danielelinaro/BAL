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
 *  \brief Definition of abstract class balInterpolator. 
 */

#ifndef _BALINTERPOLATOR_
#define _BALINTERPOLATOR_

#include <cmath>
#include <cstdio>
#include "balObject.h"

/**
 * \class balInterpolator 
 * \brief Base class for generic interpolator (also smooth approximator) of functions f: R^{nd} -> R^{nf}.
 *
 * \sa balInterp1D, balInterp2D, balInterp3D
 */

class balInterpolator : public balObject {
  public:
  
    /** Returns the name of the class. */
    virtual const char * getClassName() const;
    
    /** Destroys a balInterpolator. */
    virtual void Destroy();
    
    virtual balInterpolator * Clone() const = 0;
    
    /** Evaluates the function in point x. The result is stored in array y. 
     *  If an error has occurred the return value is -1, otherwise it is 0. */ 
    virtual int Evaluate(double *x,double *y) = 0;
    
    /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x nd). 
     *  If an error has occurred the return value is -1, otherwise it is 0. */ 
    virtual int EvaluateDerivative(double *x, double **y) = 0;
    
    /** Initializes the balInterpolator. You have to perform this operation before evaluating the function. 
     *  If an error has occurred the return value is -1, otherwise it is 0. */
    virtual int Init();

    /** Gets the dimension of the function domain (nd) */
    int GetDomainDimensions();
    
    /** Gets the dimension of the function codomain (nf) */
    int GetCodomainDimensions();

  protected:
    /** Constructor of balInterpolator. */
    balInterpolator();
    
    /** Destructor of balInterpolator. */
    virtual ~balInterpolator();

    /** Copy constructor of balInterpolator. */
    balInterpolator(const balInterpolator & interp);

    int nnd, nnf;

};
#endif
