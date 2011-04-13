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

BifurcationParameters::BifurcationParameters() {
  plower = Parameters::Create();
  pupper = Parameters::Create();
  nsteps = NULL;
  isteps = NULL;
  steps = NULL;
  _dealloc = false;
}

BifurcationParameters::~BifurcationParameters() {
  plower->Destroy();
  pupper->Destroy();
  if(_dealloc) {
    delete nsteps;
    delete isteps;
    delete steps;
  }
}

BifurcationParameters * BifurcationParameters::Create () {
  return new BifurcationParameters;
}

void BifurcationParameters::Destroy () {
  delete this;
}

const char * BifurcationParameters::GetClassName () const {
  return "BifurcationParameters";
}

void BifurcationParameters::SetNumber(int n) {
  if(n > 0) {
    Parameters::SetNumber(n);
    plower->SetNumber(n);
    pupper->SetNumber(n);
    if(_dealloc) {
      delete nsteps;
      delete isteps;
      delete steps;
    }
    nsteps = new int[n];
    isteps = new int[n];
    steps = new double[n];
    _dealloc = true;
    for(int i=0; i<n; i++) {
      plower->At(i) = 0.0;
      pupper->At(i) = 0.0;
      nsteps[i] = 1;
    }
    Setup();
  }
}

void BifurcationParameters::SetParameterBounds(Parameters * lower, Parameters * upper) {
  if(lower->GetNumber() == upper->GetNumber()) {
    SetNumber(lower->GetNumber());
    for(int i=0; i<lower->GetNumber(); i++) {
      plower->At(i) = lower->At(i);
      pupper->At(i) = upper->At(i);
    }
    Setup();
  }
}

bool BifurcationParameters::SetIthParameterLowerBound(int i, double p) {
  if(i<0 || i>=plower->GetNumber())
    return false;
  plower->At(i) = p;
  Setup();
  return true;
}

bool BifurcationParameters::SetIthParameter(int i, double p) {
  if(i<0 || i>=pupper->GetNumber())
    return false;
  plower->At(i) = p;
  pupper->At(i) = p;
  nsteps[i] = 1;
  Setup();
  return true;
}

bool BifurcationParameters::SetIthParameterUpperBound(int i, double p) {
  if(i<0 || i>=pupper->GetNumber())
    return false;
  pupper->At(i) = p;
  Setup();
  return true;
}

double BifurcationParameters::GetIthParameterLowerBound(int i) throw(Exception) {
  if(i<0 || i>=plower->GetNumber())
    throw Exception("Index out of range");
  return plower->At(i);
}

double BifurcationParameters::GetIthParameter(int i) throw(Exception) {
  if(i<0 || i>=GetNumber())
    throw Exception("Index out of range");
  return At(i);
}

double BifurcationParameters::GetIthParameterUpperBound(int i) throw(Exception) {
  if(i<0 || i>=pupper->GetNumber())
    throw Exception("Index out of range");
  return pupper->At(i);
}

bool BifurcationParameters::SetNumberOfSteps(int i, int s) {
  if(i>=0 && i<plower->GetNumber() && s>0)
    nsteps[i] = s;
  else
    return false;
  Setup();
  return true;
}

void BifurcationParameters::SetNumberOfSteps(const int * s) {
  for(int i=0; i<plower->GetNumber(); i++)
    nsteps[i] = s[i];
  Setup();
}

int BifurcationParameters::GetNumberOfSteps(int i) const {
  if(i>=0 && i<plower->GetNumber())
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
  for(int i=0; i<plower->GetNumber(); i++) {
    At(i) = plower->At(i);
    if(nsteps[i] == 1) {
      steps[i] = (pupper->At(i) - plower->At(i));
    }
    else {
      steps[i] = (pupper->At(i) - plower->At(i)) / (nsteps[i]-1);
    }
    total *= nsteps[i];
    isteps[i] = 1;
  }
  count = 1;
}

bool BifurcationParameters::Next() {
  count++;
  if(count <= total) {
    for(int i=0; i<plower->GetNumber(); i++) {
      At(i) = At(i) + steps[i];
      isteps[i]++;
      if(isteps[i] > nsteps[i]) {
	At(i) = plower->At(i);
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

