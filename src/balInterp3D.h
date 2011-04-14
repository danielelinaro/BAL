/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp3D.h
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
 *  \file balInterp3D.h
 *  \brief Definition of classes BaseInterp3D LinearInterp3D 
 */

#ifndef _BALINTERP3D_
#define _BALINTERP3D_

#include <cmath>
#include <cstdio>
#include "balObject.h"
#include "balInterp2D.h"

namespace bal {

/**
 * \class BaseInterp3D 
 * \brief Base class for three dimensional interpolation.
 * \sa Interpolator
 */
class BaseInterp3D : public Interpolator {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a BaseInterp3D. */
  virtual void Destroy();
  
  
   /** Sets the interpolation points.
   *  xi1 is a vector containing the first components of the interpolation points
   *  xi2 is a vector containing the second components of the interpolation points
   *  xi3 is a vector containing the third components of the interpolation points
   *  yi is a matrix containing the different values of the funcion in correspondence to each input point. 
   *  yi must have dimensions nf x (nx1 x nx2 x nx3) such as yi[k][i+nx1*j+nx1*nx2*l] = fk(xi1[i],xi2[j],xi3[l]).
   *  xi1 is the number of interpolation points along the first dimension (number of element of array xi1).
   *  xi2 is the number of interpolation points along the second dimension (number of element of array xi2).
   *  xi3 is the number of interpolation points along the third dimension (number of element of array xi3).
   *  nf is the number of function values corresponding to each input point (dimension of the function codomain) */
  virtual void SetInterpolationPoints(double * xi1, double * xi2, double * xi3, double **yi, int nx1, int nx2, int nx3, int nf);
  
 protected:
 
 /** Constructor of BaseInterp3D. */
  BaseInterp3D();
 
 /** Copy constructor of BaseInterp3D. */
  BaseInterp3D(const BaseInterp3D & interp);
 
 /** Destructor of BaseInterp3D. */
  ~BaseInterp3D();
  
  int nnx1, nnx2, nnx3;
  double *xx1, *xx2, *xx3, **yy;
};

/**
 * \class LinearInterp3D 
 * \brief Class for three dimensional linear interpolation.
 * \sa BaseInterp3D
 */
class LinearInterp3D : public BaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a LinearInterp3D. */
  virtual void Destroy();
  
  /** Creates a LinearInterp3D. */
  static LinearInterp3D * Create();
 
  /** Copies a LinearInterp3D */
  static LinearInterp3D * Copy(LinearInterp3D *interp);
  
  virtual LinearInterp3D * Clone() const; 
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateJacobian(double *x, double **y);  
  
  /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDivergence(double *x, double *y);  
  
  /** Initializes the LinearInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of LinearInterp3D. */
  LinearInterp3D();
 
 /** Copy constructor of LinearInterp3D. */
  LinearInterp3D(const LinearInterp3D & interp);
 
 /** Destructor of LinearInterp3D. */
  ~LinearInterp3D();
  
  LinearInterp1D *x3terp;
  LinearInterp2D ** interpsxy;
  LinearInterp1D * interpz;
  double ***yloc;
  
};

/**
 * \class PolyInterp3D 
 * \brief Class for three dimensional polinomial interpolation.
 * \sa BaseInterp3D
 */
class PolyInterp3D : public BaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a PolyInterp3D. */
  virtual void Destroy();
  
  /** Creates a PolyInterp3D. */
  static PolyInterp3D * Create();
 
  /** Copies a PolyInterp3D */
  static PolyInterp3D * Copy(PolyInterp3D *interp);
  
  virtual PolyInterp3D * Clone() const; 
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateJacobian(double *x, double **y);  
  
  /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDivergence(double *x, double *y);  
 
  /** Sets the interpolation order for each dimension, i.e. the order of the polinomials used to interpolate. Default: m1 = m2 = 2 (linear interpolation) */
  void SetInterpolationOrder(int m1, int m2, int m3);
  
  /** Initializes the PolyInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of PolyInterp3D. */
  PolyInterp3D();
 
 /** Copy constructor of PolyInterp3D. */
  PolyInterp3D(const PolyInterp3D & interp);
 
 /** Destructor of PolyInterp3D. */
  ~PolyInterp3D();

  int mm1, mm2, mm3;
  PolyInterp1D *x3terp;
  PolyInterp2D ** interpsxy;
  PolyInterp1D * interpz;
  double ***yloc;
  
};

/**
 * \class SplineInterp3D 
 * \brief Class for three dimensional spline interpolation.
 * \sa BaseInterp3D
 */
class SplineInterp3D : public BaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a SplineInterp3D. */
  virtual void Destroy();
  
  /** Creates a SplineInterp3D. */
  static SplineInterp3D * Create();
 
  /** Copies a SplineInterp3D */
  static SplineInterp3D * Copy(SplineInterp3D *interp);
  
  virtual SplineInterp3D * Clone() const; 
  
  /** Sets the dimension of the interpolation windows (how many points you want to use to perform the interpolation) */
  void SetWindow(int w);
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateJacobian(double *x, double **y);  
  
  /** Evaluates (if nd = nf) the divergence of the vector field in point x. The result is stored in y (scalar). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDivergence(double *x, double *y);  
 
  /** Initializes the SplineInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of SplineInterp3D. */
  SplineInterp3D();
 
 /** Copy constructor of SplineInterp3D. */
  SplineInterp3D(const SplineInterp3D & interp);
 
 /** Destructor of SplineInterp3D. */
  ~SplineInterp3D();

  PolyInterp1D *x3terp;
  SplineInterp2D ** interpsxy;
  SplineInterp1D * interpz;
  int window;
  double ***yloc;
};

/**
 * \class SmoothingSplineInterp3D 
 * \brief Class for three dimensional approximation with smoothing splines.
 * \sa BaseInterp3D
 */
class SmoothingSplineInterp3D : public BaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a SmoothingSplineInterp3D. */
  virtual void Destroy();
  
  /** Creates a SmoothingSplineInterp3D. */
  static SmoothingSplineInterp3D * Create();
 
  /** Copies a SmoothingSplineInterp3D */
  static SmoothingSplineInterp3D * Copy(SmoothingSplineInterp3D *interp);
  
  virtual SmoothingSplineInterp3D * Clone() const; 
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
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
  
  /** Initializes the SmoothingSplineInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of SmoothingSplineInterp3D. */
  SmoothingSplineInterp3D();
 
 /** Copy constructor of SmoothingSplineInterp3D. */
  SmoothingSplineInterp3D(const SmoothingSplineInterp3D & interp);
 
 /** Destructor of SmoothingSplineInterp3D. */
  ~SmoothingSplineInterp3D();

  PolyInterp1D *x3terp;
  SmoothingSplineInterp2D ** interpsxy;
  SmoothingSplineInterp1D * interpz;
  int window;
  double ***yloc;
  double SS;
};

} // namespace bal

#endif


