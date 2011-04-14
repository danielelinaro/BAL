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

#include "balSolution.h"

namespace bal {

Solution::Solution() {
  rows = columns = 0;
  buffer = NULL;
	lyapunov_exponents = NULL;
	spectrum_dimension = 0;
  parameters = NULL;
  ID = 0;
	lyapunov_mode = false;
}

Solution::Solution(const Solution& solution) {
  Solution();
  solution.GetSize(&rows,&columns);
  SetSize(rows,columns);
  memcpy(buffer,solution.buffer,rows*columns*sizeof(realtype));
	nturns = solution.nturns;
  ID = solution.ID;
	if (solution.lyapunov_mode)
		SetLyapunovExponents(solution.spectrum_dimension,solution.lyapunov_exponents);
}

Solution::~Solution() {
  if(buffer != NULL) delete [] buffer;
	if(lyapunov_exponents !=NULL) delete [] lyapunov_exponents;
  if(parameters != NULL) parameters->Destroy();
}

const char * Solution::GetClassName() const {
  return "Solution";
}

Solution * Solution::Create() {
  return new Solution;
}

Solution * Solution::Copy(Solution * solution) {
  return new Solution(*solution);
}

void Solution::Destroy() {
  delete this;
}

int Solution::GetRows() const {
  return rows;
}

int Solution::GetColumns() const {
  return columns;
}

void Solution::SetSize(int r, int c) {
  if(buffer != NULL) delete [] buffer;
  rows = r;
  columns = c;
  buffer = new realtype[rows*columns];
}

void Solution::GetSize(int * r, int * c) const {
  *r = rows;
  *c= columns;
}

void Solution::SetData(int r, int c, realtype * data) {
  SetSize(r,c);
  memcpy(buffer,data,rows*columns*sizeof(realtype));
}

realtype * Solution::GetData() const {
  return buffer;
}	

void Solution::SetParameters(Parameters * p) {
  if(parameters != NULL) parameters->Destroy();
  parameters = Parameters::Copy(p);
}

Parameters * Solution::GetParameters() const {
  return parameters;
}

void Solution::SetNumberOfTurns(int _nturns) {
  nturns = _nturns;
}
  
int Solution::GetNumberOfTurns() const {
  return nturns;
}

void Solution::SetLyapunovExponents(int n, realtype * lp) {
	if (lyapunov_exponents == NULL || n != spectrum_dimension)
		lyapunov_exponents = new realtype[n];
	spectrum_dimension = n;
	lyapunov_mode = true;
	for (int i=0; i<spectrum_dimension ; i++)
		lyapunov_exponents[i] = lp[i];
}

realtype * Solution::GetLyapunovExponents() const {
	return lyapunov_exponents;
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

bool CompareBalSolutions(Solution *sol1, Solution *sol2) {
  return sol1->GetID() < sol2->GetID();
}

} // namespace bal

