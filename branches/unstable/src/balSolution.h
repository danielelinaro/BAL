/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balSolution.h
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
 * \file balSolution.h
 * \brief Definition of the class Solution
 */

#ifndef _BALSOLUTION_
#define _BALSOLUTION_

#include "balObject.h"
#include "balCommon.h"
#include "balParameters.h"

namespace bal {

/**
 * \class Solution
 * \brief Class that contains the result of an integration of a dynamical
 * system
 * \sa DynamicalSystem Parameters ODESolver
 */
class Solution : public Object {
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  /** Creates a new Solution. */
  static Solution * Create();
  /** Copies a Solution */
  static Solution * Copy(Solution * solution);
  /** Destroys a Solution. */
  virtual void Destroy();
  
  int GetRows() const;
  int GetColumns() const;
  void GetSize(int * r, int * c) const;
  void SetSize(int r, int c);
  
  Parameters * GetParameters() const;
  void SetParameters(Parameters * p);
  
  realtype * GetData() const;
  void SetData(int r, int c, realtype * data);
  
  int GetNumberOfTurns() const;
  void SetNumberOfTurns(int _nturns);

	realtype * GetLyapunovExponents() const;
  void SetLyapunovExponents(int n, realtype * lp);
	
  int GetID() const;
  void SetID(int id);
	
	bool IsLyapunovMode() const;

 protected:
  /* Protected destructor of the class. */
  virtual ~Solution();
  Solution();
  Solution(const Solution & solution);
  
 private:
  Parameters * parameters;
  realtype * buffer;
	realtype * lyapunov_exponents;
	int spectrum_dimension;
  int rows, columns;
  int nturns;
  int ID;
	bool lyapunov_mode;
};

bool CompareBalSolutions(Solution *sol1, Solution *sol2);

} // namespace bal

#endif

