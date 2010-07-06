/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balBifurcationParameters.h
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
 * \file balBifurcationParameters.h
 * \brief Class for storing the tuples of parameters required in the
 * computation of a bifurcation diagram.
 */

#ifndef _BALBIFURCATIONPARAMETERS_
#define _BALBIFURCATIONPARAMETERS_

#include "balParameters.h"
#include "balCommon.h"

/**
 * \class balBifurcationParameters
 * \brief Class for storing the tuples of parameters required in the
 * computation of a bifurcation diagram.
 *
 * This class is used by balBifurcationDiagram to iterate over all the
 * possible tuples of parameters needed to compute a bifurcation diagram.
 * Suppose that the user wants to integrate a system that depends on four
 * parameters: if they want to compute a bidimensional diagram, then it is
 * sufficient to specify the upper and lower bounds of the 4 parameters.
 * Those that do not change will have the same upper and lower bounds and
 * a number of steps equal to one. The 2 parameters that have to change
 * will have distinct values for the upper and lower bounds and a number of
 * steps greater than one. If the number of steps along the two parameters
 * is n1 and n2 respectively, then the total number of tuples given by
 * balBifurcationParameters will be n1*n2.
 * An instance of balBifurcationParameters stores internally the current
 * tuple of parameters, so that a balDynamicalSystem can use an instance of
 * balBifurcationParameters as if it were of type balParameters.
 *
 * \see balParameters balBifurcationDiagram balDynamicalSystem
 */
class balBifurcationParameters : public balParameters {
 public:
  virtual const char * GetClassName () const { return "balBifurcationParameters"; }
  static balBifurcationParameters * Create () { return new balBifurcationParameters; }
  virtual void Destroy () { delete this; }
  
  virtual void SetNumber(int n);
  
  bool SetIthParameterLowerBound(int i, double p);
  bool SetIthParameter(int i, double p);
  bool SetIthParameterUpperBound(int i, double p);
  void SetParameterBounds(balParameters * lower, balParameters * upper);
  double GetIthParameterUpperBound(int i) throw(balException);
  double GetIthParameter(int i) throw(balException);
  double GetIthParameterLowerBound(int i) throw(balException);
  
  bool SetNumberOfSteps(int i, int s);
  void SetNumberOfSteps(const int * s);
  int GetNumberOfSteps(int i) const;
  
  inline void Reset() { Setup(); }
  int GetTotalNumberOfTuples() const { return total; }
  
  bool Next();
  bool HasTuples() const;
  bool HasNext() const;
  bool IsFirst() const;
  bool IsLast() const;
  
 protected:
  balBifurcationParameters();
  ~balBifurcationParameters();
  void Setup();
  
 private:
  /** The lower bounds of the parameters. */
  balParameters * plower;
  /** The upper bounds of the parameters. */
  balParameters * pupper;
  /** The parameter steps for every parameter. */
  double * steps;
  /** The number of steps associated with every parameter. */
  int * nsteps;
  /** The current steps associated with every parameter. */
  int * isteps;
  /** The total number of steps, i.e. the product of the values contained
   * in nsteps. */
  int total;
  /** The current parameters' tuple. */
  int count;

  /** Tells whether memory has been allocated or not */
  bool _dealloc;
};

#endif

