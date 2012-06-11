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
 *  \brief Definition of classes BaseInterp1D, LinearInterp1D, PolyInterp1D and SplineInterp1D. 
 */

#ifndef _BALINTERP1D_
#define _BALINTERP1D_

#include <cmath>
#include <iostream>
#include "balInterpolator.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define NATURAL_SPLINE 1.0e99
#define FINITE_DIFFERENCES_STEP 1.0e-6

namespace bal {

/**
 * \class BaseInterp1D 
 * \brief Base class for one dimensional interpolation of vector functions \f$f: R \rightarrow R^{m}\f$.
 *
 * \example interp.cpp 
 *
 * \sa Interpolator
 */
class BaseInterp1D : public Interpolator {
		
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;

  /** Copies a BaseInterp1D. */
  static BaseInterp1D * Copy(const BaseInterp1D & interp);
  

  /** Destroys a BaseInterp1D. */
  virtual void Destroy();

 /**
  * Sets interpolation points.
  * \param xi Vector containing the abscissas of the points.
  * \param yi Matrix containing the different ordinates corresponding to each ascissa. \f$y_i\f$ must have dimension \f$m \times length\f$. 
  * \param length Number of interpolation points.
  * \param nf Number of ordinates for each abscissa (dimension \f$m\f$ of the function codomain).
  */
  virtual void SetInterpolationPoints(double * xi, double **yi, int length, int nf);
  
  /**
   * \param x Generic domain point.
   * \return Returns a value \f$j\f$ such that \f$x\f$ is (insofar as possible) centered in the subrange 
   * \f$xx[j,\dots,j+mm-1]\f$, where \f$xx\f$ is the vector of stored points and \f$mm\f$ is the interpolation order. 
   * The returned value is not less than 0, nor greater than n-1.
   * 
   * Before calling this function, the values in \f$xx\f$ must be already monotonic, either increasing or decreasing.
   *
   * Taken from <a href="http://www.nr.com">Numerical Recipes: 
   * The art of scientific computing" [Third Edition]</a>, page 115.
   */
  int Locate(const double x);
  
  /**
   * \param x Generic domain point.
   * \return Returns a value \f$j\f$ such that \f$x\f$ is (insofar as possible) centered in the subrange 
   * \f$xx[j,\dots,j+mm-1]\f$, where \f$xx\f$ is the vector of stored points and \f$mm\f$ is the interpolation order. 
   * The returned value is not less than 0, nor greater than n-1.
   * 
   * Before calling this function, the values in \f$xx\f$ must be already monotonic, either increasing or decreasing.
   *
   * Taken from <a href="http://www.nr.com/">Numerical Recipes: 
   * The art of scientific computing" [Third Edition]</a>, page 116.
   */
  int Hunt(const double x);

  /** Returns true if it is better you use Hunt() instead of Locate() */
  bool nextHunt();
  
 protected:

 /** Constructor of BaseInterp1D. */
  BaseInterp1D();

 /** Copy constructor of BaseInterp1D. */
  BaseInterp1D(const BaseInterp1D & interp);
 
 /** Destructor of BaseInterp1D. */
  ~BaseInterp1D();
  
  int n, mm, jsav, cor, dj; 
  double *xx, **yy; 
};


/**
 * \class LinearInterp1D 
 * \brief Class for one dimensional linear interpolation of vector functions \f$f: R \rightarrow R^{m}\f$.
 *
 * \sa BaseInterp1D
 */
class LinearInterp1D : public BaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a LinearInterp1D. */
  virtual void Destroy();
  
  /** Creates a LinearInterp1D */
  static LinearInterp1D * Create();
  
  /** Copies a LinearInterp1D */
  static LinearInterp1D * Copy(LinearInterp1D *interp);
  virtual LinearInterp1D * Clone() const;
  
	virtual int Evaluate(double *x, double *y);
  
  virtual int EvaluateJacobian(double *x, double **y);

  virtual int EvaluateDivergence(double *x, double *y);  

 protected:
 
  /** Constructor of LinearInterp1D. */
  LinearInterp1D();
  
  /** Copy constructor of LinearInterp1D. */
  LinearInterp1D(const LinearInterp1D & interp);

  /** Destructor of LinearInterp1D. */
  ~LinearInterp1D();
  
};

/**
 * \class PolyInterp1D 
 * \brief Class for one dimensional polynomial interpolation of vector functions \f$f: R \rightarrow R^{m}\f$.
 *
 * \sa BaseInterp1D
 */
class PolyInterp1D : public BaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a PolyInterp1D. */
  virtual void Destroy();
  
  /** Creates a PolyInterp1D */
  static PolyInterp1D * Create();
  
  /** Copies a PolyInterp1D */
  static PolyInterp1D * Copy(PolyInterp1D *interp);
  virtual PolyInterp1D * Clone() const;

  /** Sets the interpolation order, i.e. the order of the polinomial used to interpolate. Default: 2 (linear interpolation) */
  void SetInterpolationOrder(int m);

  virtual int Evaluate(double *x, double *y);

  virtual int EvaluateJacobian(double *x, double **y);

  virtual int EvaluateDivergence(double *x, double *y);  
  
 protected:
  
  /** Constructor of PolyInterp1D. */
  PolyInterp1D();

  /** Copy constructor of PolyInterp1D. */
  PolyInterp1D(const PolyInterp1D & interp);
  
  /** Destructor of PolyInterp1D. */
  ~PolyInterp1D();
  
 private:
  double dy;
};

/**
 * \class SplineInterp1D 
 * \brief Class for one dimensional interpolation of vector functions \f$f: R \rightarrow R^{m}\f$ using splines.
 *
 * \sa BaseInterp1D
 */
class SplineInterp1D : public BaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a SplineInterp1D. */
  virtual void Destroy();
  
  /** Creates a SplineInterp1D */
  static SplineInterp1D * Create();
  
  /** Copies a SplineInterp1D */
  static SplineInterp1D * Copy(SplineInterp1D *interp);
  virtual SplineInterp1D * Clone() const;
  
  /** Sets the value of the first derivative at the boundaries of the domain. There is a single value for each component of the codomain.
   *  Default value: natural splines (second derivative equal to 0 in the boundaries) */
  void SetBoundaryConditions(double yp1, double ypn);

  int Init();
  
  virtual int Evaluate(double *x, double *y);

  virtual int EvaluateJacobian(double *x, double **y);

  virtual int EvaluateDivergence(double *x, double *y); 

  
 protected:
  
  /** Constructor of SplineInterp1D. */
  SplineInterp1D();
  
  /** Copy constructor of SplineInterp1D. */
  SplineInterp1D(const SplineInterp1D & interp);
  
  /** Destructor of SplineInterp1D. */
  ~SplineInterp1D();
 
  /** Method used by Init(). */
  void Sety2();

 private:
  double **y2,yyp1,yypn;
};

/**
 * \class SmoothingSplineInterp1D 
 * \brief Class for one dimensional approximation of vector functions \f$f: R \rightarrow R^{m}\f$ using smoothing Reinsch splines.
 * Notice that an interpolation is not performed.
 *
 * \sa BaseInterp1D
 */
class SmoothingSplineInterp1D : public BaseInterp1D {
 public:
  
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  
  /** Destroys a SmoothingSplineInterp1D. */
  virtual void Destroy();

  /** Creates a SmoothingSplineInterp1D. */
  static SmoothingSplineInterp1D * Create();
  virtual SmoothingSplineInterp1D * Clone() const;
  
  /** Copies a SmoothingSplineInterp1D. */
  static SmoothingSplineInterp1D * Copy(SmoothingSplineInterp1D *interp);
  
  /** Sets the smoothing factors. The smoothing spline is computed by finding a third order polinomial which minimizes the
   *  integral of the square of the second derivative under these constraints:
   * \f[
   * \sum_i{( \frac{f(x_i) - y_i}{dy_i})^2} \leq S
   * \f]
   * If \f$S=0\f$ or \f$dy_i=0\f$ for each \f$i\f$, an interpolation is performed.
   * By increasing \f$S\f$ or \f$dy_i\f$ values the function becomes smoother and smoother.
   */
  void SetSmoothingParameters(double *dy, double S);
   
  /** Sets the smoothing factors. The smoothing spline is computed by finding a third order polinomial which minimizes the
   *  integral of the square of the second derivative under these constraints:
   * \f[
   * \sum_i{(f(x_i) - y_i)^2} \leq smooth^2
   * \f]
   * If \f$smooth = 0\f$ an interpolation is performed.
   * By increasing \f$smooth\f$ the function becomes smoother and smoother.
   */
  void SetSmoothingParameters(double smooth);

  int Init();

  virtual int Evaluate(double *x, double *y);

  virtual int EvaluateJacobian(double *x, double **y);

  virtual int EvaluateDivergence(double *x, double *y);  
  
 protected:
  
  /** Constructor of SmoothingSplineInterp1D. */
  SmoothingSplineInterp1D();
  
  /** Copy constructor of SmoothingSplineInterp1D. */
  SmoothingSplineInterp1D(const SmoothingSplineInterp1D & interp);
  
  /** Destructor of SmoothingSplineInterp1D. */
  ~SmoothingSplineInterp1D();
  
  /** Method used by Init(). */
  int ComputeCoefficients(int i);
  
 private:
  double *ddy,SS,**a,**b,**c,**d;
  int _DEALLOC_DDY;
};

} // namespace bal

#endif

