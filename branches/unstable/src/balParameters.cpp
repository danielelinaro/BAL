/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    Parameters.cpp
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
 * \file Parameters.cpp
 * \brief Implementation of the class Parameters
 */

#include "balParameters.h"
#include <sstream>

namespace bal {

Parameters::Parameters(int np) : p(np), pars(new double[np]) {
  for(int i=0; i<p; i++)
    pars[i] = 0.0;
}

Parameters::Parameters(const Parameters& params) : p(params.p), pars(new double[params.p]) {
  for(int i=0; i<p; i++)
    pars[i] = params.pars[i];
}

Parameters::~Parameters () {
}

Object* Parameters::Clone() const {
  return new Parameters(*this);
}

int Parameters::GetNumber () const {
  return p;
}

double& Parameters::operator[] (int k) {
  return pars[k];
}

void Parameters::operator= (const Parameters& params) {
  if(p != params.p)
    throw Exception("Wrong number of parameters in DynamicalSystem::SetParameters");
  for(int i=0; i<p; i++)
    pars[i] = params.pars[i];
}

double Parameters::At (int k) const {
  return pars[k];
}

double* Parameters::GetParameters () const {
  return pars.get();
}

std::string Parameters::ToString() const {
  std::stringstream ss;
  ss << "(";
  if(p > 0) {
    for(int i=0; i<p-1; i++)
      ss << pars[i] << ",";
    ss << pars[p-1];
  }
  ss << ")";
  return ss.str();
}

} // namespace bal
