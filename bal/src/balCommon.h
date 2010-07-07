/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balCommon.h
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
 * \file balCommon.h
 * \brief Header file with common definitions for other bal classes.
 */

/**
 * \mainpage A Bifurcation Analysis Library (BAL)
 * \author Daniele Linaro <daniele.linaro@unige.it>
 * \version 0.9.1
 * \date 2009 - 2010
 *
 * This library aims to provide an easy way to compute brute-force
 * bifurcation diagrams of (autonomous) dynamical systems described 
 * by sets of ordinary differential equations (ODE's).
 */

#ifndef _BALCOMMON_
#define _BALCOMMON_

#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#ifdef CVODE26
#include <sundials/sundials_direct.h>
#endif
#include <nvector/nvector_serial.h>

/** Ith numbers components 0..NEQ-1 */
#define Ith(v,i)    NV_Ith_S(v,i)
/** IJth numbers rows,cols 0..NEQ-1 */
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

/* colors */
#define ESC ''
#define RED "[31m"
#define GREEN "[32m"
#define YELLOW "[33m"
#define BLUE "[34m"
#define MAGENTA "[35m"
#define CYAN "[36m"
#define NORMAL "[00m"

#define PI				(3.141592653589793)

#define	LIST_MAX_SIZE	1000

#include <exception>

/**
 * \class balException
 * \brief Class that implements the exceptions thrown by the BAL library
 */
class balException : public std::exception {
 public:
 balException(const char* description = NULL) : errorDescription(description) {}
  virtual const char* what() const throw() {
    if(errorDescription == NULL)
      return "balException";
    return errorDescription;
  }
 private:
  const char* errorDescription;
};

#endif

