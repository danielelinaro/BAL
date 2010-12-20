/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp1D.h
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
 *  \brief Definition of classes balBaseInterp1D, balLinearInterp1D, balPolyInterp1D and balSplineInterp1D. 
 */

#ifndef _BALINTERP1D_
#define _BALINTERP1D_

#include <cmath>
#include <iostream>
#include "balInterpolator.h"
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define NATURAL_SPLINE 1.0e99
#define FINITE_DIFFERENCES_STEP 1.0e-6
/**
 * \class balBaseInterp1D 
 * \brief Base class for one dimensional interpolation of vector functions f: R -> R^{nf}.
 *
 * \sa balInterpolator
 */
class balBaseInterp1D : public balInterpolator {
		
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;

  /** Copies a balBaseInterp1D. */
  static balBaseInterp1D * Copy(const balBaseInterp1D & interp);
  

  /** Destroys a balBaseInterp1D. */
  virtual void Destroy();
  
  /** Sets the interpolation points.
   *  xi is a vector containing the abscissas of the points
   *  yi is a matrix containing the different ordinates corresponding to each ascissa. yi must have dimensions nf x length.
   *  length is the number of interpolation points.
   *  nf is the number of ordinates for each abscissa (dimension of the function codomain) */
  virtual void SetInterpolationPoints(double * xi, double **yi, int length, int nf);
  
  /**
   * Given a value x, return a value j such that x is (insofar as possible) centered in the subrange 
   * xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either 
   * increasing or decreasing. The returned value is not less than 0, nor greater than n-1. 
   *
   * Taken from "Numerical Recipes: The art of scientific computing". Third
   * Edition, page 115.
   */
  int Locate(const double x);
  
  /**
   * Given a value x, return a value j such that x is (insofar as possible) centered in the subrange 
   * xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either 
   * increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
   *
   * Taken from "Numerical Recipes: The art of scientific computing". Third
   * Edition, page 116.
   */
  int Hunt(const double x);

  /** Returns true if it is better you use Hunt() instead of Locate() */
  bool nextHunt();
  
 protected:

 /** Constructor of balBaseInterp1D. */
  balBaseInterp1D();

 /** Copy constructor of balBaseInterp1D. */
  balBaseInterp1D(const balBaseInterp1D & interp);
 
 /** Destructor of balBaseInterp1D. */
  ~balBaseInterp1D();
  
  int n, mm, jsav, cor, dj; 
  double *xx, **yy; 
};


/**
 * \class balLinearInterp1D 
 * \brief Class for one dimensional linear interpolation of vector functions f: R -> R^{nf}.
 *
 * \sa balBaseInterp1D
 */
class balLinearInterp1D : public balBaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balLinearInterp1D. */
  virtual void Destroy();
  
  /** Creates a balLinearInterp1D */
  static balLinearInterp1D * Create();
  
  /** Copies a balLinearInterp1D */
  static balLinearInterp1D * Copy(balLinearInterp1D *interp);
  virtual balLinearInterp1D * Clone() const;
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 1). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);

 protected:
 
  /** Constructor of balLinearInterp1D. */
  balLinearInterp1D();
  
  /** Copy constructor of balLinearInterp1D. */
  balLinearInterp1D(const balLinearInterp1D & interp);

  /** Destructor of balLinearInterp1D. */
  ~balLinearInterp1D();
  
};

/**
 * \class balPolyInterp1D 
 * \brief Class for one dimensional polynomial interpolation of vector functions f: R -> R^{nf}.
.
 * \sa balBaseInterp1D
 */
class balPolyInterp1D : public balBaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balPolyInterp1D. */
  virtual void Destroy();
  
  /** Creates a balPolyInterp1D */
  static balPolyInterp1D * Create();
  
  /** Copies a balPolyInterp1D */
  static balPolyInterp1D * Copy(balPolyInterp1D *interp);
  virtual balPolyInterp1D * Clone() const;

  /** Sets the interpolation order, i.e. the order of the polinomial used to interpolate. Default: 2 (linear interpolation) */
  void SetInterpolationOrder(int m);

  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 1).
   *  If an error has occurred the return value is -1, otherwise it is 0. 
   *  Implemented with finite differences. */
  virtual int EvaluateDerivative(double *x, double **y);
  
 protected:
  
  /** Constructor of balPolyInterp1D. */
  balPolyInterp1D();

  /** Copy constructor of balPolyInterp1D. */
  balPolyInterp1D(const balPolyInterp1D & interp);
  
  /** Destructor of balPolyInterp1D. */
  ~balPolyInterp1D();
  
 private:
  double dy;
};

/**
 * \class balSplineInterp1D 
 * \brief Class for one dimensional interpolation of vector functions f: R -> R^{nf} by using splines.
 *
 * \sa balBaseInterp1D
 */
class balSplineInterp1D : public balBaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balSplineInterp1D. */
  virtual void Destroy();
  
  /** Creates a balSplineInterp1D */
  static balSplineInterp1D * Create();
  
  /** Copies a balSplineInterp1D */
  static balSplineInterp1D * Copy(balSplineInterp1D *interp);
  virtual balSplineInterp1D * Clone() const;
  
  /** Sets the value of the first derivative at the boundaries of the domain. There is a single value for each component of the codomain.
   *  Default value: natural splines (second derivative equal to 0 in the boundaries) */
  void SetBoundaryConditions(double yp1, double ypn);

  /** Initializes the balSpline1D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 1). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);

  
 protected:
  
  /** Constructor of balSplineInterp1D. */
  balSplineInterp1D();
  
  /** Copy constructor of balSplineInterp1D. */
  balSplineInterp1D(const balSplineInterp1D & interp);
  
  /** Destructor of balSplineInterp1D. */
  ~balSplineInterp1D();
 
  /** Method used by Init(). */
  void Sety2();

 private:
  double **y2,yyp1,yypn;
};

/**
 * \class balSmoothingSplineInterp1D 
 * \brief Class for one dimensional approximation of vector functions f: R -> R^{nf} by using smoothing Reinsch splines.
 *        Notice that an interpolation is not performed.
 * \sa balBaseInterp1D
 */
class balSmoothingSplineInterp1D : public balBaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balSmoothingSplineInterp1D. */
  virtual void Destroy();

  /** Creates a balSmoothingSplineInterp1D. */
  static balSmoothingSplineInterp1D * Create();
  virtual balSmoothingSplineInterp1D * Clone() const;
  
  /** Copies a balSmoothingSplineInterp1D. */
  static balSmoothingSplineInterp1D * Copy(balSmoothingSplineInterp1D *interp);
  
  /** Sets the smoothing factors. The smoothing spline is computed by finding a third order polinomial which minimizes the
   *  integral of the square of the second derivative under these constraints:
   *  __
   *  \
   *  /_ ( (f(x_i) - y_i) / dy_i )^2 <= S
   *
   * If S = 0, or all dy values = 0, an interpolation is performed.
   * By increasing S or dy values the function becomes smoother and smoother.
   */
  void SetSmoothingParameters(double *dy, double S);
   
  /** Sets the smoothing factors. The smoothing spline is computed by finding a third order polinomial which minimizes the
   *  integral of the square of the second derivative under these constraints:
   *  __
   *  \
   *  /_ (f(x_i) - y_i)^2 <= smooth^2
   *
   * If smooth = 0 an interpolation is performed.
   * By increasing smooth the function becomes smoother and smoother.
   */
  void SetSmoothingParameters(double smooth);
  
  /** Initializes the balSmoothingSplineInterp1D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 1). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);
  
 protected:
  
  /** Constructor of balSmoothingSplineInterp1D. */
  balSmoothingSplineInterp1D();
  
  /** Copy constructor of balSmoothingSplineInterp1D. */
  balSmoothingSplineInterp1D(const balSmoothingSplineInterp1D & interp);
  
  /** Destructor of balSmoothingSplineInterp1D. */
  ~balSmoothingSplineInterp1D();
  
  /** Method used by Init(). */
  int ComputeCoefficients(int i);
  
 private:
  double *ddy,SS,**a,**b,**c,**d;
  int _DEALLOC_DDY;
};
#endif


