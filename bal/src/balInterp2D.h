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
 *  \file Interp2D.h
 *  \brief Definition of classes BaseInterp2D LinearInterp2D 
 */

#ifndef _BALINTERP2D_
#define _BALINTERP2D_

#include <cmath>
#include <cstdio>
#include "balObject.h"
#include "balInterp1D.h"

namespace bal {

/**
 * \class BaseInterp2D 
 * \brief Base class for two dimensional interpolation.
 * \sa Interpolator
 */
class BaseInterp2D : public Interpolator {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a BaseInterp2D. */
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
 
 /** Constructor of BaseInterp2D. */
  BaseInterp2D();
 
 /** Copy constructor of BaseInterp2D. */
  BaseInterp2D(const BaseInterp2D & interp);
 
 /** Destructor of BaseInterp2D. */
  ~BaseInterp2D();
  
  int nnx1, nnx2;
  double *xx1, *xx2, **yy;
};

/**
 * \class LinearInterp2D 
 * \brief Class for two dimensional linear interpolation.
 * \sa BaseInterp2D
 */
class LinearInterp2D : public BaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a LinearInterp2D. */
  virtual void Destroy();
  
  /** Creates a LinearInterp2D. */
  static LinearInterp2D * Create();
  
  LinearInterp2D * Clone() const;
 
  /** Copies a LinearInterp2D */
  static LinearInterp2D * Copy(LinearInterp2D *interp);
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateJacobian(double *x, double **y);  
  
  /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDivergence(double *x, double *y);  
  
  /** Initializes the LinearInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of LinearInterp2D. */
  LinearInterp2D();
 
 /** Copy constructor of LinearInterp2D. */
  LinearInterp2D(const LinearInterp2D & interp);
 
 /** Destructor of LinearInterp2D. */
  ~LinearInterp2D();
  
  LinearInterp1D *x1terp, *x2terp;
  
};

/**
 * \class PolyInterp2D 
 * \brief Class for two dimensional polinomial interpolation.
 * \sa BaseInterp2D
 */
class PolyInterp2D : public BaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a PolyInterp2D. */
  virtual void Destroy();
  
  /** Creates a PolyInterp2D. */
  static PolyInterp2D * Create();
 
  /** Copies a PolyInterp2D */
  static PolyInterp2D * Copy(PolyInterp2D *interp);
  
  PolyInterp2D * Clone() const;
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateJacobian(double *x, double **y);  
  
  /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDivergence(double *x, double *y);  
 
  /** Sets the interpolation order for each dimension, i.e. the order of the polinomials used to interpolate. Default: m1 = m2 = 2 (linear interpolation) */
  void SetInterpolationOrder(int m1, int m2);
  
  /** Initializes the PolyInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of PolyInterp2D. */
  PolyInterp2D();
 
 /** Copy constructor of PolyInterp2D. */
  PolyInterp2D(const PolyInterp2D & interp);
 
 /** Destructor of PolyInterp2D. */
  ~PolyInterp2D();

  int mm1, mm2;
  PolyInterp1D ** interpsx, *interpy;
  PolyInterp1D *x2terp;
  double ***yloc;
  
};

/**
 * \class SplineInterp2D 
 * \brief Class for two dimensional spline interpolation.
 * \sa BaseInterp2D
 */
class SplineInterp2D : public BaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a SplineInterp2D. */
  virtual void Destroy();
  
  /** Creates a SplineInterp2D. */
  static SplineInterp2D * Create();
 
  /** Copies a SplineInterp2D */
  static SplineInterp2D * Copy(SplineInterp2D *interp);
  
  SplineInterp2D * Clone() const;
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateJacobian(double *x, double **y);  
  
  /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDivergence(double *x, double *y);  
 
  /** Initializes the SplineInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of SplineInterp2D. */
  SplineInterp2D();
 
 /** Copy constructor of SplineInterp2D. */
  SplineInterp2D(const SplineInterp2D & interp);
 
 /** Destructor of SplineInterp2D. */
  ~SplineInterp2D();

  PolyInterp1D *x1terp, *x2terp;
  int wt[16][16];
  double ***y1d, ***y2d, ***y12d, ***c;
};

/**
 * \class SmoothingSplineInterp2D 
 * \brief Class for two dimensional approximation with smoothing splines.
 * \sa BaseInterp2D
 */
class SmoothingSplineInterp2D : public BaseInterp2D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a SmoothingSplineInterp2D. */
  virtual void Destroy();
  
  /** Creates a SmoothingSplineInterp2D. */
  static SmoothingSplineInterp2D * Create();
 
  /** Copies a SmoothingSplineInterp2D */
  static SmoothingSplineInterp2D * Copy(SmoothingSplineInterp2D *interp);
  
  SmoothingSplineInterp2D * Clone() const;
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 2). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateJacobian(double *x, double **y);  
  
  /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDivergence(double *x, double *y);  
  
  /** Sets the smoothing factor. If S = 0 an interpolation is performed.
  * By increasing S the function becomes smoother and smoother.
  * Default: S = 0 (interpolation).
  */
  void SetSmoothingParameter(double S);

  /** Sets the dimension of the interpolation window (how many points you want to use to perform the interpolation) */
  void SetWindow(int w);
  
  /** Initializes the SmoothingSplineInterp2D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of SmoothingSplineInterp2D. */
  SmoothingSplineInterp2D();
 
 /** Copy constructor of SmoothingSplineInterp2D. */
  SmoothingSplineInterp2D(const SmoothingSplineInterp2D & interp);
 
 /** Destructor of SmoothingSplineInterp2D. */
  ~SmoothingSplineInterp2D();

  int window;
  double SS;
  SmoothingSplineInterp1D ** interpsx, *interpy;
  PolyInterp1D *x2terp;
  double ***yloc;
  
};

} // namespace bal

#endif


