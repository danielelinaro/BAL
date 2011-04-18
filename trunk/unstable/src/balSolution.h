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

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include "balObject.h"
#include "balCommon.h"
#include "balParameters.h"

namespace bal {

/**
 * \class Solution
 * \brief Class that contains the result of an integration of a dynamical
 * system.
 * \sa DynamicalSystem Parameters ODESolver
 */
class Solution : public Object {
public:
  Solution(int r, int c, realtype *buf);
  Solution(const Solution& solution);
  virtual ~Solution();

  virtual Object* Clone() const;
  std::string ToString() const;

  int GetRows() const;
  int GetColumns() const;
  void GetSize(int *r, int *c) const;
  void SetSize(int r, int c);
  
  Parameters* GetParameters() const;
  void SetParameters(const Parameters *p);
  
  realtype* GetData() const;
  
  int GetNumberOfTurns() const;
  void SetNumberOfTurns(int nturns);

  realtype* GetLyapunovExponents() const;
  void SetLyapunovExponents(const realtype *lp);
  bool IsLyapunovMode() const;
	
  int GetID() const;
  void SetID(int id);
  
  bool operator< (const Solution& sol) const;

 private:
  boost::shared_ptr<Parameters> parameters;
  boost::shared_array<realtype> buffer;
  boost::shared_array<realtype> lyapunov_exponents;
  int spectrum_dimension;
  int rows, columns;
  int nturns;
  int ID;
  bool lyapunov_mode;
};

bool CompareSolutions(Solution *sol1, Solution *sol2);

} // namespace bal

#endif


