/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp2D.h
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
 *  \file balInterp2D.h
 *  \brief Definition of classes balBaseInterp2D balLinearInterp2D 
 */

#ifndef _BALINTERP2D_
#define _BALINTERP2D_

#include <cmath>
#include <cstdio>
#include "balObject.h"
#include "balInterp1D.h"
using namespace std;

/**
 * \class balBaseInterp2D 
 * \brief Base class for two dimensional interpolation.
 * \sa balInterpolator
 */
class balBaseInterp2D : public balInterpolator {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balBaseInterp2D. */
  virtual void Destroy();
 
    /** Sets the interpolation points.
   *  xi1 is a vector containing the first components of the interpolation points
   *  xi2 is a vector containing the second components of the interpolation points
   *  yi is a matrix containing the different values of the funcion in correspondence to each input point. 
   *  yi must have dimensions nf x (nx1 x nx2) such as yi[k][i+nx1*j] = fk(xi1[i],xi2[j]).
   *  xi1 is the number of interpolation points along the first dimension (number of element of array xi1).
   *  xi2 is the number of interpolation points along the second dimension (number of element of array xi2).
   *  nf is the number of function values corresponding to each input point (dimension of the function codomain) */
  virtual void SetInterpolationPoints(double * xi1, double * xi2, double **yi, int nx1, int nx2, int nf);
  
 protected:
 
 /** Constructor of balBaseInterp2D. */
  balBaseInterp2D();
 
 /** Copy constructor of balBaseInterp2D. */
  balBaseInterp2D(const balBaseInterp2D & interp);
 
 /** Destructor of balBaseInterp2D. */
  ~balBaseInterp2D();
  
  int nnx1, nnx2;
  double *xx1, *xx2, **yy;
};

/**
 * \class balLinearInterp2D 
 * \brief Class for two dimensional linear interpolation.
 * \sa balBaseInterp2D
 */
class balLinearInterp2D : public balBaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balLinearInterp2D. */
  virtual void Destroy();
  
  /** Creates a balLinearInterp2D. */
  static balLinearInterp2D * Create();
  
  balLinearInterp2D * Clone() const;
 
  /** Copies a balLinearInterp2D */
  static balLinearInterp2D * Copy(balLinearInterp2D *interp);
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
  
  /** Initializes the balLinearInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balLinearInterp2D. */
  balLinearInterp2D();
 
 /** Copy constructor of balLinearInterp2D. */
  balLinearInterp2D(const balLinearInterp2D & interp);
 
 /** Destructor of balLinearInterp2D. */
  ~balLinearInterp2D();
  
  balLinearInterp1D *x1terp, *x2terp;
  
};

/**
 * \class balPolyInterp2D 
 * \brief Class for two dimensional polinomial interpolation.
 * \sa balBaseInterp2D
 */
class balPolyInterp2D : public balBaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balPolyInterp2D. */
  virtual void Destroy();
  
  /** Creates a balPolyInterp2D. */
  static balPolyInterp2D * Create();
 
  /** Copies a balPolyInterp2D */
  static balPolyInterp2D * Copy(balPolyInterp2D *interp);
  
  balPolyInterp2D * Clone() const;
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
 
  /** Sets the interpolation order for each dimension, i.e. the order of the polinomials used to interpolate. Default: m1 = m2 = 2 (linear interpolation) */
  void SetInterpolationOrder(int m1, int m2);
  
  /** Initializes the balPolyInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balPolyInterp2D. */
  balPolyInterp2D();
 
 /** Copy constructor of balPolyInterp2D. */
  balPolyInterp2D(const balPolyInterp2D & interp);
 
 /** Destructor of balPolyInterp2D. */
  ~balPolyInterp2D();

  int mm1, mm2;
  balPolyInterp1D ** interpsx, *interpy;
  balPolyInterp1D *x2terp;
  double ***yloc;
  
};

/**
 * \class balSplineInterp2D 
 * \brief Class for two dimensional spline interpolation.
 * \sa balBaseInterp2D
 */
class balSplineInterp2D : public balBaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balSplineInterp2D. */
  virtual void Destroy();
  
  /** Creates a balSplineInterp2D. */
  static balSplineInterp2D * Create();
 
  /** Copies a balSplineInterp2D */
  static balSplineInterp2D * Copy(balSplineInterp2D *interp);
  
  balSplineInterp2D * Clone() const;
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
 
  /** Initializes the balSplineInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balSplineInterp2D. */
  balSplineInterp2D();
 
 /** Copy constructor of balSplineInterp2D. */
  balSplineInterp2D(const balSplineInterp2D & interp);
 
 /** Destructor of balSplineInterp2D. */
  ~balSplineInterp2D();

  balPolyInterp1D *x1terp, *x2terp;
  int wt[16][16];
  double ***y1d, ***y2d, ***y12d, ***c;
};

/**
 * \class balSmoothingSplineInterp2D 
 * \brief Class for two dimensional approximation with smoothing splines.
 * \sa balBaseInterp2D
 */
class balSmoothingSplineInterp2D : public balBaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balSmoothingSplineInterp2D. */
  virtual void Destroy();
  
  /** Creates a balSmoothingSplineInterp2D. */
  static balSmoothingSplineInterp2D * Create();
 
  /** Copies a balSmoothingSplineInterp2D */
  static balSmoothingSplineInterp2D * Copy(balSmoothingSplineInterp2D *interp);
  
  balSmoothingSplineInterp2D * Clone() const;
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
  
  /** Sets the smoothing factor. If S = 0 an interpolation is performed.
  * By increasing S the function becomes smoother and smoother.
  * Default: S = 0 (interpolation).
  */
  void SetSmoothingParameter(double S);

  /** Sets the dimension of the interpolation window (how many points you want to use to perform the interpolation) */
  void SetWindow(int w);
  
  /** Initializes the balSmoothingSplineInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balSmoothingSplineInterp2D. */
  balSmoothingSplineInterp2D();
 
 /** Copy constructor of balSmoothingSplineInterp2D. */
  balSmoothingSplineInterp2D(const balSmoothingSplineInterp2D & interp);
 
 /** Destructor of balSmoothingSplineInterp2D. */
  ~balSmoothingSplineInterp2D();

  int window;
  double SS;
  balSmoothingSplineInterp1D ** interpsx, *interpy;
  balPolyInterp1D *x2terp;
  double ***yloc;
  
};
#endif


