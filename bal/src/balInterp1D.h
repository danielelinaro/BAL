/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterp1D.h
 *
 *   Copyright (C) 2009 Daniele Linaro
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
 *  \brief Classes for one dimensional interpolation.
 *  
 */

#ifndef _BALINTERP1D_
#define _BALINTERP1D_

#include <cmath>
#include <cstdio>
#include "balObject.h"
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

class balBaseInterp1D : public balObject {
		
	public:
		/** Returns the name of the class. */
		virtual const char * GetClassName() const { return "balBaseInterp1D" ; }
		/** Destroys a balBaseInterp1D. */
		virtual void Destroy() { this->~balBaseInterp1D(); }

		double interp(double x) {
			int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo,x);
		}

		/**
		 * Given a value x, return a value j such that x is (insofar as possible) centered in the subrange 
		 * xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either 
		 * increasing or decreasing. The returned value is not less than 0, nor greater than n-1. 
		 *
		 * Taken from "Numerical Recipes: The art of scientific computing". Third
		 * Edition, page 115.
		 */
		int locate(const double x);
		
		/**
		 *	Given a value x, return a value j such that x is (insofar as possible) centered in the subrange 
		 *	xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either 
		 *	increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
		 *
		 * Taken from "Numerical Recipes: The art of scientific computing". Third
		 * Edition, page 116.
		 */
		int hunt(const double x);

	protected:
		balBaseInterp1D(double * x, const double *y, int length, int m)
			: n(length), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) { 
				dj = MIN(1, (int)pow((double)n,0.25));
		}
		~balBaseInterp1D() {}

		virtual double rawinterp(int jlo, double x) = 0;

	public:
		int n, mm, jsav, cor, dj; 
		const double *xx, *yy; 
};


class balLinearInterp1D : public balBaseInterp1D {
	public:
		/** Returns the name of the class. */
		virtual const char * GetClassName() const { return "balLinearInterp1D" ; }
		/** Destroys a balLinearInterp1D. */
		virtual void Destroy() { this->~balLinearInterp1D(); }
		/** Creates a balLinearInterp1D */
		static balLinearInterp1D * Create(double * xv, double * yv, int length) {
			return new balLinearInterp1D(xv,yv,length);
		}

	protected:
		balLinearInterp1D(double * xv, double * yv, int length) :
			balBaseInterp1D(xv,yv,length,2) {}
		~balLinearInterp1D() {}

		virtual double rawinterp(int j, double x) {
			if (xx[j]==xx[j+1]) 
				return yy[j];
			return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
		}
};

class balPolyInterp1D : public balBaseInterp1D {
	public:
		/** Returns the name of the class. */
		virtual const char * GetClassName() const { return "balPolyInterp1D" ; }
		/** Destroys a balPolyInterp1D. */
		virtual void Destroy() { this->~balPolyInterp1D(); }
		/** Creates a balLinearInterp1D */
		static balPolyInterp1D * Create(double * xv, double * yv, int length, int m) {
			return new balPolyInterp1D(xv,yv,length,m);
		}

	protected:
		balPolyInterp1D(double * xv, double * yv, int length, int m) :
			balBaseInterp1D(xv,yv,length,m), dy(0.0) {}
		~balPolyInterp1D() {}

		virtual double rawinterp(int j, double x);

	private:
		double dy;
};

class balSplineInterp1D : public balBaseInterp1D {
	public:
		/** Returns the name of the class. */
		virtual const char * GetClassName() const { return "balSplineInterp1D" ; }
		/** Destroys a balSplineInterp1D. */
		virtual void Destroy() { this->~balSplineInterp1D(); }
		/** Creates a balLinearInterp1D */
		static balSplineInterp1D * Create(double * xv, double * yv, int length, double yp1=1.e99, double ypn=1.e99) {
			return new balSplineInterp1D(xv,yv,length,yp1,ypn);
		}

	protected:
		balSplineInterp1D(double * xv, double * yv, int length, double yp1, double ypn) :
			balBaseInterp1D(xv,yv,length,2) {
			y2 = new double[length];
			sety2(xv,yv,yp1,ypn);
		}
		~balSplineInterp1D() { delete y2; }

		void sety2(double *xv, double *yv, double yp1, double ypn);
		virtual double rawinterp(int j, double x);

	private:
		double *y2;
};

#endif


