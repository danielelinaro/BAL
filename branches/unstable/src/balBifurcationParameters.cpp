/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balBifurcationParameters.cpp
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
 * \file balBifurcationParameters.cpp
 * \brief Implementation of the class BifurcationParameters
 */

#include "balBifurcationParameters.h"

namespace bal {

BifurcationParameters::BifurcationParameters(int np) : Parameters(np), plower(np), pupper(np),
						       steps(new double[np]),
						       nsteps(new int[np]),
						       isteps(new int[np]) {
  for(int i=0; i<np; i++) {
    steps[i] = 0.0;
    nsteps[i] = 1;
    isteps[i] = 0;
  }
  std::cout << "BifurcationParameters constructor.\n";
}

BifurcationParameters::BifurcationParameters(const BifurcationParameters& bp) : Parameters(bp.p),
										plower(bp.p),
										pupper(bp.p),
										steps(new double[bp.p]),
										nsteps(new int[bp.p]),
										isteps(new int[bp.p]) {
  for(int i=0; i<p; i++) {
    steps[i] = bp.steps[i];
    nsteps[i] = bp.nsteps[i];
    isteps[i] = bp.isteps[i];
  }
  std::cout << "BifurcationParameters copy constructor.\n";
}

BifurcationParameters::~BifurcationParameters() {
  std::cout << "BifurcationParameters destructor.\n";
}

Object* BifurcationParameters::Clone() const {
  return new BifurcationParameters(*this);
}

void BifurcationParameters::SetParameterBounds(const Parameters& lower, const Parameters& upper) {
  plower = lower;
  pupper = upper;
  Setup();
}

void BifurcationParameters::SetIthParameterLowerBound(int i, double p) {
  plower[i] = p;
  Setup();
}

void BifurcationParameters::SetIthParameter(int i, double p) {
  plower[i] = p;
  pupper[i] = p;
  nsteps[i] = 1;
  Setup();
}

void BifurcationParameters::SetIthParameterUpperBound(int i, double p) {
  pupper[i] = p;
  Setup();
}

double BifurcationParameters::GetIthParameterLowerBound(int i) {
  return plower[i];
}

double BifurcationParameters::GetIthParameter(int i) {
  return At(i);
}

double BifurcationParameters::GetIthParameterUpperBound(int i) {
  return pupper[i];
}

bool BifurcationParameters::SetNumberOfSteps(int i, int s) {
  if(i>=0 && i<p && s>0)
    nsteps[i] = s;
  else
    return false;
  Setup();
  return true;
}

void BifurcationParameters::SetNumberOfSteps(const int *s) {
  for(int i=0; i<p; i++)
    nsteps[i] = s[i];
  Setup();
}

int BifurcationParameters::GetNumberOfSteps(int i) const {
  if(i>=0 && i<p)
    return nsteps[i];
  return -1;
}

int BifurcationParameters::GetTotalNumberOfTuples() const {
  return total;
}

void BifurcationParameters::Reset() {
  Setup();
}

void BifurcationParameters::Setup() {
  total = 1;
  for(int i=0; i<p; i++) {
    (*this)[i] = plower[i];
    if(nsteps[i] == 1) {
      steps[i] = (pupper[i] - plower[i]);
    }
    else {
      steps[i] = (pupper[i] - plower[i]) / (nsteps[i]-1);
    }
    total *= nsteps[i];
    isteps[i] = 1;
  }
  count = 1;
}

bool BifurcationParameters::Next() {
  count++;
  if(count <= total) {
    for(int i=0; i<p; i++) {
      (*this)[i] += steps[i];
      isteps[i]++;
      if(isteps[i] > nsteps[i]) {
	(*this)[i] = plower[i];
	isteps[i] = 1;
      }
      else {
	break;
      }
    }
    return true;
  }
  return false;
}

bool BifurcationParameters::HasTuples() const {
  return count <= total;
}

bool BifurcationParameters::HasNext() const {
  return count < total;
}

bool BifurcationParameters::IsFirst() const {
  return count == 1;
}

bool BifurcationParameters::IsLast() const {
  return count == total;
}

} // namespace bal

