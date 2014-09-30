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
 *  \brief Definition of classes BaseInterp2D and derivates.
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
 * \brief Base class for two dimensional interpolation of vector functions \f$f: R^{2} \rightarrow R^{m}\f$.
 * \sa Interpolator
 */
class BaseInterp2D : public Interpolator {
  
 public:
 /** Constructor of BaseInterp2D. */
  BaseInterp2D();
 
 /** Copy constructor of BaseInterp2D. */
  BaseInterp2D(const BaseInterp2D & interp);
 
 /** Destructor of BaseInterp2D. */
  virtual ~BaseInterp2D();

 	/** Sets the interpolation points.
	 * \param xi1 First components of the interpolation points.
	 * \param xi2 Second components of the interpolation points.
	 * \param yi  Matrix containing the different values of the function in correspondence to each input point.
	 *			  It must have dimensions \f$ m \times (n_{x1} \times n_{x2})\f$ such as \f$y_i[k][i+n_{x1}*j] = f_k(x_i^1[i],x_i^2[j])\f$.
	 * \param nx1 Number of interpolation points along the first dimension (dimension of array xi1).
	 * \param nx2 Number of interpolation points along the second dimension (dimension of array xi2).
	 * \param nf Number of function values corresponding to each input point (dimension \f$m\f$ of the function codomain).
	 */
  virtual void SetInterpolationPoints(double * xi1, double * xi2, double **yi, int nx1, int nx2, int nf);
  
 protected:
  
  int nnx1, nnx2;
  double *xx1, *xx2, **yy;
};

/**
 * \class LinearInterp2D 
 * \brief Class for two dimensional linear interpolation of vector functions \f$f: R^2 \rightarrow R^{m}\f$.
 * \sa BaseInterp2D
 */
class LinearInterp2D : public BaseInterp2D {
  
 public:
 /** Constructor of LinearInterp2D. */
  LinearInterp2D();
 
 /** Copy constructor of LinearInterp2D. */
  LinearInterp2D(const LinearInterp2D & interp);
 
 /** Destructor of LinearInterp2D. */
  ~LinearInterp2D();

  virtual int Evaluate(double *x, double *y);
  
  virtual int EvaluateJacobian(double *x, double **y);  

  virtual int EvaluateDivergence(double *x, double *y);  

  int Init();

 protected:
  
  LinearInterp1D *x1terp, *x2terp;
  
};

/**
 * \class PolyInterp2D 
 * \brief Class for two dimensional polinomial interpolation of vector functions \f$f: R^2 \rightarrow R^{m}\f$.
 * \sa BaseInterp2D
 */
class PolyInterp2D : public BaseInterp2D {
  
 public:

 /** Constructor of PolyInterp2D. */
  PolyInterp2D();
 
 /** Copy constructor of PolyInterp2D. */
  PolyInterp2D(const PolyInterp2D & interp);
 
 /** Destructor of PolyInterp2D. */
  ~PolyInterp2D();

  virtual int Evaluate(double *x, double *y);
  
  virtual int EvaluateJacobian(double *x, double **y);  
  
  virtual int EvaluateDivergence(double *x, double *y);  
 
  /** Sets the interpolation order for each dimension, i.e.\ the order of the polinomials used to interpolate. Default: m1 = m2 = 2 (linear interpolation) */
  void SetInterpolationOrder(int m1, int m2);

  int Init();

 protected:
 
  int mm1, mm2;
  PolyInterp1D ** interpsx, *interpy;
  PolyInterp1D *x2terp;
  double ***yloc;
  
};

/**
 * \class SplineInterp2D 
 * \brief Class for two dimensional interpolation of vector functions \f$f: R^2 \rightarrow R^{m}\f$ using splines.
 * \sa BaseInterp2D
 */
class SplineInterp2D : public BaseInterp2D {
  
 public:
 
 /** Constructor of SplineInterp2D. */
  SplineInterp2D();
 
 /** Copy constructor of SplineInterp2D. */
  SplineInterp2D(const SplineInterp2D & interp);
 
 /** Destructor of SplineInterp2D. */
  ~SplineInterp2D();

  virtual int Evaluate(double *x, double *y);

  virtual int EvaluateJacobian(double *x, double **y);  

  virtual int EvaluateDivergence(double *x, double *y);  
 
  int Init();

 protected:
 
  PolyInterp1D *x1terp, *x2terp;
  int wt[16][16];
  double ***y1d, ***y2d, ***y12d, ***c;
};

/**
 * \class SmoothingSplineInterp2D 
 * \brief Class for two dimensional approximation of vector functions \f$f: R^2 \rightarrow R^{m}\f$ using smoothing splines.
 * \sa BaseInterp2D
 */
class SmoothingSplineInterp2D : public BaseInterp2D {
  
 public:
 
 /** Constructor of SmoothingSplineInterp2D. */
  SmoothingSplineInterp2D();
 
 /** Copy constructor of SmoothingSplineInterp2D. */
  SmoothingSplineInterp2D(const SmoothingSplineInterp2D & interp);
 
 /** Destructor of SmoothingSplineInterp2D. */
  ~SmoothingSplineInterp2D();

  virtual int Evaluate(double *x, double *y);
  
  virtual int EvaluateJacobian(double *x, double **y);  

  virtual int EvaluateDivergence(double *x, double *y);  
  
  /** Sets the smoothing factor. If S = 0 an interpolation is performed.
  * By increasing S the function becomes smoother and smoother.
  * Default: S = 0 (interpolation).
  */
  void SetSmoothingParameter(double S);

  /** Sets the dimension of the interpolation window.
   * \param w Number of points the user wants to use to perform the interpolation. */
  void SetWindow(int w);
  
  int Init();

 protected:
 
  int window;
  double SS;
  SmoothingSplineInterp1D ** interpsx, *interpy;
  PolyInterp1D *x2terp;
  double ***yloc;
  
};

} // namespace bal

#endif


