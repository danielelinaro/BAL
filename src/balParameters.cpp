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
#include <iostream>
#include <sstream>

namespace bal {

Parameters::Parameters(int np) : p(np) {
#ifdef DEBUG
  std::cout << "Parameters constructor.\n";
#endif
  if (p) {
    pars = new double[p];
    for(int i=0; i<p; i++)
      pars[i] = 0.0;
  }
}

Parameters::Parameters(const Parameters& params) : p(params.p) {
#ifdef DEBUG
  std::cout << "Parameters copy constructor.\n";
#endif
  if (p) {
    pars = new double[p];
    for(int i=0; i<p; i++)
      pars[i] = params.pars[i];
  }
}

Parameters::~Parameters () {
#ifdef DEBUG
  std::cout << "Parameters destructor.\n";
#endif
  if (p)
    delete pars;
}

int Parameters::GetNumber () const {
  return p;
}

double& Parameters::operator[] (int k) {
  if (k>=0 && k<p)
    return pars[k];
  throw "Parameters::operator[] - index out of bounds";
}

void Parameters::operator= (const Parameters& params) {
  if(p != params.p)
    throw Exception("Parameters::operator= - wrong number of parameters");
  for(int i=0; i<p; i++)
    pars[i] = params.pars[i];
}

double& Parameters::At (int k) {
  if (k>=0 && k<p)
    return pars[k];
  throw "Parameters::At - index out of bounds";
}

double* Parameters::GetParameters () const {
  return pars;
}

Parameters* Parameters::Clone() const {
  return new Parameters(*this);
}

void Parameters::CopyValues(const Parameters* params){
  if(p != params->p)
    throw Exception("Parameters::CopyValues - wrong number of parameters");
  for (int i=0; i<p; i++)
    pars[i] = params->pars[i];
}

std::ostream& operator<< (std::ostream& os, Parameters& pars) {
  os << "[";
  for (int i=0; i<pars.GetNumber()-1; i++)
    os << pars[i] << ",";
  os << pars[pars.GetNumber()-1] << "]";
  return os;
}

} // namespace bal

