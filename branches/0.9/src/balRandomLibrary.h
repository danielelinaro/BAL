/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balRandomLibrary.h
 *
 *   Copyright (C) 2010 Daniele Linaro
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
 * \file balRandomLibrary.h
 * \brief Definition of the classes UniformRandom and NormalBM
 */

#ifndef _BALRANDLIB_
#define _BALRANDLIB_

#include <cmath>

typedef unsigned long long ullong;
typedef unsigned int uint;

namespace bal {

namespace random {

/**
 * \class UniformRandom
 * \brief Random number generator uniformly distributed in the interval [0,1]. Taken from Numerical Recipes in C++
 */
class UniformRandom { 
 private:
  ullong u,v,w; 
 public:
  UniformRandom(ullong j) : v(4101842887655102017LL), w(1) { 
    u = j ^ v; int64(); 
    v = u; int64(); 
    w = v; int64(); 
  } 
  
  inline ullong int64() { 
    u = u * 2862933555777941757LL + 7046029254386353087LL; 
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8; 
    w = 4294957665U*(w & 0xffffffff) + (w >> 32); 
    ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4; 
    return (x + v) ^ w;
  }
  inline double doub() {
    return 5.42101086242752217E-20 * int64();
  }
  inline uint int32() {
    return (uint) int64();
  }
};

/**
 * \class NormalBM
 * \brief Gaussian random number generator with given mean and variance. Taken from Numerical Recipes in C++, it uses
 * the Box-Muller algorithm.
 */
class NormalBM : UniformRandom {
 private:
  double mu,sig;
  double storedval;
 public:
 NormalBM(double mmu, double ssig, ullong i)
   : UniformRandom(i), mu(mmu), sig(ssig), storedval(0.) {}
  
  double deviate() {
    double v1,v2,rsq,fac;
    if (storedval == 0.) {
      do { 
	v1 = 2.0*doub()-1.0;
	v2 = 2.0*doub()-1.0; 
	rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0); 
      fac = sqrt(-2.0*log(rsq)/rsq);
      storedval = v1*fac; 
      return mu + sig*v2*fac; 
    } else {
      fac = storedval; 
      storedval = 0.; 
      return mu + sig*fac; 
    } 
  } 
};

} // namespace random

} // namespace bal

#endif

