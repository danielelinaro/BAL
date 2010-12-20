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
 *  \brief Definition of classes balBaseInterp3D balLinearInterp3D 
 */

#ifndef _BALINTERP3D_
#define _BALINTERP3D_

#include <cmath>
#include <cstdio>
#include "balObject.h"
#include "balInterp2D.h"
using namespace std;

/**
 * \class balBaseInterp3D 
 * \brief Base class for three dimensional interpolation.
 * \sa balInterpolator
 */
class balBaseInterp3D : public balInterpolator {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balBaseInterp3D. */
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
 
 /** Constructor of balBaseInterp3D. */
  balBaseInterp3D();
 
 /** Copy constructor of balBaseInterp3D. */
  balBaseInterp3D(const balBaseInterp3D & interp);
 
 /** Destructor of balBaseInterp3D. */
  ~balBaseInterp3D();
  
  int nnx1, nnx2, nnx3;
  double *xx1, *xx2, *xx3, **yy;
};

/**
 * \class balLinearInterp3D 
 * \brief Class for three dimensional linear interpolation.
 * \sa balBaseInterp3D
 */
class balLinearInterp3D : public balBaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balLinearInterp3D. */
  virtual void Destroy();
  
  /** Creates a balLinearInterp3D. */
  static balLinearInterp3D * Create();
 
  /** Copies a balLinearInterp3D */
  static balLinearInterp3D * Copy(balLinearInterp3D *interp);
  
  virtual balLinearInterp3D * Clone() const; 
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
  
  /** Initializes the balLinearInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balLinearInterp3D. */
  balLinearInterp3D();
 
 /** Copy constructor of balLinearInterp3D. */
  balLinearInterp3D(const balLinearInterp3D & interp);
 
 /** Destructor of balLinearInterp3D. */
  ~balLinearInterp3D();
  
  balLinearInterp1D *x3terp;
  balLinearInterp2D ** interpsxy;
  balLinearInterp1D * interpz;
  double ***yloc;
  
};

/**
 * \class balPolyInterp3D 
 * \brief Class for three dimensional polinomial interpolation.
 * \sa balBaseInterp3D
 */
class balPolyInterp3D : public balBaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balPolyInterp3D. */
  virtual void Destroy();
  
  /** Creates a balPolyInterp3D. */
  static balPolyInterp3D * Create();
 
  /** Copies a balPolyInterp3D */
  static balPolyInterp3D * Copy(balPolyInterp3D *interp);
  
  virtual balPolyInterp3D * Clone() const; 
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
 
  /** Sets the interpolation order for each dimension, i.e. the order of the polinomials used to interpolate. Default: m1 = m2 = 2 (linear interpolation) */
  void SetInterpolationOrder(int m1, int m2, int m3);
  
  /** Initializes the balPolyInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balPolyInterp3D. */
  balPolyInterp3D();
 
 /** Copy constructor of balPolyInterp3D. */
  balPolyInterp3D(const balPolyInterp3D & interp);
 
 /** Destructor of balPolyInterp3D. */
  ~balPolyInterp3D();

  int mm1, mm2, mm3;
  balPolyInterp1D *x3terp;
  balPolyInterp2D ** interpsxy;
  balPolyInterp1D * interpz;
  double ***yloc;
  
};

/**
 * \class balSplineInterp3D 
 * \brief Class for three dimensional spline interpolation.
 * \sa balBaseInterp3D
 */
class balSplineInterp3D : public balBaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balSplineInterp3D. */
  virtual void Destroy();
  
  /** Creates a balSplineInterp3D. */
  static balSplineInterp3D * Create();
 
  /** Copies a balSplineInterp3D */
  static balSplineInterp3D * Copy(balSplineInterp3D *interp);
  
  virtual balSplineInterp3D * Clone() const; 
  
  /** Sets the dimension of the interpolation windows (how many points you want to use to perform the interpolation) */
  void SetWindow(int w);
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
 
  /** Initializes the balSplineInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balSplineInterp3D. */
  balSplineInterp3D();
 
 /** Copy constructor of balSplineInterp3D. */
  balSplineInterp3D(const balSplineInterp3D & interp);
 
 /** Destructor of balSplineInterp3D. */
  ~balSplineInterp3D();

  balPolyInterp1D *x3terp;
  balSplineInterp2D ** interpsxy;
  balSplineInterp1D * interpz;
  int window;
  double ***yloc;
};

/**
 * \class balSmoothingSplineInterp3D 
 * \brief Class for three dimensional approximation with smoothing splines.
 * \sa balBaseInterp3D
 */
class balSmoothingSplineInterp3D : public balBaseInterp3D {
  
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a balSmoothingSplineInterp3D. */
  virtual void Destroy();
  
  /** Creates a balSmoothingSplineInterp3D. */
  static balSmoothingSplineInterp3D * Create();
 
  /** Copies a balSmoothingSplineInterp3D */
  static balSmoothingSplineInterp3D * Copy(balSmoothingSplineInterp3D *interp);
  
  virtual balSmoothingSplineInterp3D * Clone() const; 
  
  /** Evaluates the function in point x. The result is stored in array y. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int Evaluate(double *x, double *y);
  
  /** Evaluates the Jacobian matrix of the function in point x. The result is stored in matrix y (of dimensions nf x 3). 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  virtual int EvaluateDerivative(double *x, double **y);  
 
  /** Sets the smoothing factor. If S = 0 an interpolation is performed.
  * By increasing S the function becomes smoother and smoother.
  * Default: S = 0 (interpolation).
  */
  void SetSmoothingParameter(double S);

  /** Sets the dimension of the interpolation window (how many points you want to use to perform the interpolation) */
  void SetWindow(int w);
  
  /** Initializes the balSmoothingSplineInterp3D. You have to perform this operation before evaluating the function. 
   *  If an error has occurred the return value is -1, otherwise it is 0. */
  int Init();

 protected:
 
 /** Constructor of balSmoothingSplineInterp3D. */
  balSmoothingSplineInterp3D();
 
 /** Copy constructor of balSmoothingSplineInterp3D. */
  balSmoothingSplineInterp3D(const balSmoothingSplineInterp3D & interp);
 
 /** Destructor of balSmoothingSplineInterp3D. */
  ~balSmoothingSplineInterp3D();

  balPolyInterp1D *x3terp;
  balSmoothingSplineInterp2D ** interpsxy;
  balSmoothingSplineInterp1D * interpz;
  int window;
  double ***yloc;
  double SS;
};
#endif


