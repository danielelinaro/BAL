/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balInterpSystem.cpp
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
 * \file balInterpSystem.cpp
 * \brief Implementation of the class balInterpSystem
 */

#include "balInterpSystem.h"

balDynamicalSystem* balInterpSystemFactory() {
  return balInterpSystem::Create();
}

balInterpSystem::balInterpSystem() {
  //SetDimension(2);
  SetNumberOfParameters(0);
  //xderiv = N_VNew_Serial(GetDimension());
  interpolator = NULL;
}

balInterpSystem::balInterpSystem(const balInterpSystem& interpsystem) : balDynamicalSystem( interpsystem ) {
  if (interpsystem.interpolator != NULL)
    interpolator = interpsystem.interpolator->Clone();
  else
    interpolator = NULL;
  //xderiv = N_VNew_Serial(eye.GetDimension());
  /*for(int i = 0; i < eye.GetDimension(); i++)
    Ith(xderiv,i)=Ith(eye.xderiv,i);
  _dealloc = eye._dealloc;
  */
}

balInterpSystem::~balInterpSystem() {
  //N_VDestroy_Serial(xderiv);
}

balInterpSystem * balInterpSystem::Create () {
  return new balInterpSystem;
}

void balInterpSystem::Destroy () {
  delete this;
}

balDynamicalSystem * balInterpSystem::Copy() {
  return new balInterpSystem(*this);
}
  
const char * balInterpSystem::GetClassName () const { 
  return "balInterpSystem";
}

bool balInterpSystem::HasEvents() const {
  return true;
}

int balInterpSystem::SetInterpolator(balInterpolator *interp) {
  if (interp->GetDomainDimensions() != interp->GetCodomainDimensions()) {
    cerr<<"balInterpSystem::SetInterpolator() - Invalid interpolator\n";
    return -1;
  }
  interpolator = interp;
  SetDimension(interpolator->GetDomainDimensions());
  SetNumberOfEvents(GetDimension());
  return 0;
}

int balInterpSystem::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) { 
  int i;
  int n = GetDimension();
  int res = interpolator->Evaluate(NV_DATA_S(x),NV_DATA_S(xdot));
  if (res==-1)
    return ! CV_SUCCESS;
  //for (i=0; i<n; i++)
  //  Ith(xdot,i) = y[i];
  return CV_SUCCESS;
}

int balInterpSystem::Events (realtype t, N_Vector x, realtype * event, void * data) {
  /*RHS(t,x,xderiv,data);
  for(int i=0; i<GetNumberOfEvents(); i++)
    event[i] = Ith(xderiv,i);
  */return CV_SUCCESS;
}
