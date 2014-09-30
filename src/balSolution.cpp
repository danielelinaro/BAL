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
 * \file balSolution.cpp
 * \brief Implementation of the class Solution
 */

#include <cstring>
#include "balSolution.h"

namespace bal {

Solution::Solution(int r, int c, realtype *buf)
  : rows(r), columns(c),
    nturns(0), spectrum_dimension(c-2), ID(0),
    lyapunov_mode(false),
    buffer(new realtype[r*c]), lyapunov_exponents(new realtype[c-2]) {
#ifdef DEBUG
  std::cout << "Solution constructor.\n";
#endif
  memcpy(buffer.get(), buf, r*c*sizeof(realtype));
}

Solution::Solution(const Solution& solution) 
  : rows(solution.rows), columns(solution.columns),
    nturns(solution.nturns), spectrum_dimension(solution.spectrum_dimension),
    ID(solution.ID), lyapunov_mode(solution.lyapunov_mode),
    buffer(new realtype[solution.rows*solution.columns]),
    lyapunov_exponents(new realtype[solution.spectrum_dimension]),
    parameters(dynamic_cast<Parameters*>(solution.parameters->Clone())) {
#ifdef DEBUG
  std::cout << "Solution copy constructor.\n";
#endif
  memcpy(buffer.get(), solution.buffer.get(), rows*columns*sizeof(realtype));
  memcpy(lyapunov_exponents.get(), solution.lyapunov_exponents.get(), spectrum_dimension*sizeof(realtype));
}

Solution::~Solution() {
#ifdef DEBUG
  std::cout << "Solution destructor.\n";
#endif
}

int Solution::GetRows() const {
  return rows;
}

int Solution::GetColumns() const {
  return columns;
}

void Solution::GetSize(int *r, int *c) const {
  *r = rows;
  *c = columns;
}

realtype* Solution::GetData() const {
  return buffer.get();
}	

void Solution::SetParameters(const Parameters *p) {
  parameters = boost::shared_ptr<Parameters>(new Parameters(*p));
}

Parameters* Solution::GetParameters() const {
  return parameters.get();
}

void Solution::SetNumberOfTurns(int nturns_) {
  nturns = nturns_;
}
  
int Solution::GetNumberOfTurns() const {
  return nturns;
}

void Solution::SetLyapunovExponents(const realtype *lp) {
  memcpy(lyapunov_exponents.get(), lp, spectrum_dimension*sizeof(realtype));
}

realtype* Solution::GetLyapunovExponents() const {
  return lyapunov_exponents.get();
}

void Solution::SetID(int id) {
  ID = id;
}

int Solution::GetID() const {
  return ID;
}

bool Solution::IsLyapunovMode() const {
  return lyapunov_mode;
}

bool Solution::operator< (const Solution& sol) const {
  return GetID() < sol.GetID();
}

bool CompareSolutions(Solution *sol1, Solution *sol2) {
  return sol1->GetID() < sol2->GetID();
}

} // namespace bal

