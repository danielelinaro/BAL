/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balBifurcationParameters.h
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

#ifndef _BALBIFURCATIONPARAMETERS_
#define _BALBIFURCATIONPARAMETERS_

#include "balParameters.h"
#include "balCommon.h"

/** 
 * \file balBifurcationParameters.h
 * \brief Definition of the class BifurcationParameters
 */

namespace bal {

/**
 * \class BifurcationParameters
 * \brief Class for storing the tuples of parameters required in the
 * computation of a bifurcation diagram.
 *
 * This class is used by BifurcationDiagram to iterate over all the
 * possible tuples of parameters needed to compute a bifurcation diagram.
 *
 * Suppose that the user wants to integrate a system that depends on four
 * parameters: if they want to compute a bidimensional diagram, then it is
 * sufficient to specify the upper and lower bounds of the 4 parameters.
 * Those that do not change will have the same upper and lower bounds and
 * a number of steps equal to one. The 2 parameters that have to change
 * will have distinct values for the upper and lower bounds and a number of
 * steps greater than one. If the number of steps along the two parameters
 * is \f$n_1\f$ and \f$n_2\f$ respectively, then the total number of tuples given by
 * BifurcationParameters will be \f$n_1\f$*\f$n_2\f$.
 *
 * An instance of BifurcationParameters stores internally the current
 * tuple of parameters, so that a DynamicalSystem can use an instance of
 * BifurcationParameters as if it were of type Parameters.
 *
 * Here is an example of use:
 * it defines a vector of parameters \f$[p_0,p_1,p_2,p_3]\f$ where \f$[p_1,p_2,p_3]\f$ are set to
 * a fixed value (won't change during bifurcation analysis), while for \f$[p_0]\f$
 * is set the range of variation \f$[p_0^{min},p_0^{max}]\f$ and the number of steps
 * in which the interval is divided \f$p_0^{steps}\f$.
 * 
 * \snippet bifdiag.cpp Defining space of parameters
 *
 * \example bifparams.cpp
 * \sa Parameters BifurcationDiagram DynamicalSystem
 */
class BifurcationParameters : public Parameters {
 public:
  virtual const char * GetClassName () const;
  static BifurcationParameters * Create ();
  virtual void Destroy ();
  
  virtual void SetNumber(int n);
	
	bool SetIthParameter(int i, double p);
  bool SetIthParameterLowerBound(int i, double p);
  bool SetIthParameterUpperBound(int i, double p);
	/** \remark Parameters that don't change during analysis have the same lower and upper bound. */
  void SetParameterBounds(Parameters * lower, Parameters * upper);
  
  double GetIthParameterUpperBound(int i) throw(Exception);
  double GetIthParameter(int i) throw(Exception);
  double GetIthParameterLowerBound(int i) throw(Exception);
  
  bool SetNumberOfSteps(int i, int s);
  void SetNumberOfSteps(const int * s);
  int GetNumberOfSteps(int i) const;
  
  void Reset();
  int GetTotalNumberOfTuples() const;
  
  bool Next();
  bool HasTuples() const;
  bool HasNext() const;
  bool IsFirst() const;
  bool IsLast() const;
  
 protected:
  BifurcationParameters();
  ~BifurcationParameters();
  void Setup();
  
 private:
  /** The lower bounds of the parameters. */
  Parameters * plower;
  /** The upper bounds of the parameters. */
  Parameters * pupper;
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

} // namespace bal

#endif

