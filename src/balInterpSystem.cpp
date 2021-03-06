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
 * \brief Implementation of the class InterpSystem
 */

#include <iostream>
#include "balInterpSystem.h"

bal::DynamicalSystem* InterpSystemFactory() {
  return bal::InterpSystem::Create();
}

namespace bal {

InterpSystem::InterpSystem() {
  //SetDimension(2);
  SetNumberOfParameters(0);
  //xderiv = N_VNew_Serial(GetDimension());
  interpolator = NULL;
  backward = false;
  arclength = false;
  _dealloc = false;
}

InterpSystem::InterpSystem(const InterpSystem& interpsystem) : DynamicalSystem( interpsystem ) {
  if (interpsystem.interpolator != NULL) {
    interpolator = interpsystem.interpolator->Clone();
    _dealloc = true;
  }
  else { 
    interpolator = NULL;
    _dealloc = false;
  }
  backward = interpsystem.backward;
  arclength = interpsystem.arclength;
  //xderiv = N_VNew_Serial(eye.GetDimension());
  /*for(int i = 0; i < eye.GetDimension(); i++)
    Ith(xderiv,i)=Ith(eye.xderiv,i);
  _dealloc = eye._dealloc;
  */
}

InterpSystem::~InterpSystem() {  
  if (_dealloc)
     interpolator->Destroy();
  //N_VDestroy_Serial(xderiv);
}

InterpSystem * InterpSystem::Create () {
  return new InterpSystem;
}

void InterpSystem::Destroy () {
  delete this;
}

InterpSystem * InterpSystem::Copy(InterpSystem *is) {
  return new InterpSystem(*is);
}

DynamicalSystem * InterpSystem::Clone() const {
  return new InterpSystem(*this);
}
  
const char * InterpSystem::GetClassName () const { 
  return "InterpSystem";
}

bool InterpSystem::HasEvents() const {
  return false;
}

int InterpSystem::SetInterpolator(Interpolator *interp) {
  if (interp->GetDomainDimensions() != interp->GetCodomainDimensions()) {
          std::cerr << "InterpSystem::SetInterpolator() - Invalid interpolator\n";
    return -1;
  }
  interpolator = interp;
  SetDimension(interpolator->GetDomainDimensions());
  SetNumberOfEvents(GetDimension());
  return 0;
}

int InterpSystem::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) { 
  int i;
  int n = GetDimension();
  // NV_DATA_S gets a pointer to data conteined in N_Vector object
  int res = interpolator->Evaluate(NV_DATA_S(x),NV_DATA_S(xdot));
  if (res==-1)
    return ! CV_SUCCESS;
  //for (i=0; i<n; i++)
  //  Ith(xdot,i) = y[i];
  if (backward) {    
    for (i=0; i<n; i++) {
      Ith(xdot,i) = -Ith(xdot,i);
    }
  }
  if (arclength) {
    double norm = 0.0;
    for (i=0; i<n; i++) {
      norm += Ith(xdot,i)*Ith(xdot,i);
    }
    //norm += 1.0;
    //printf("-%e\n",norm);
    norm = sqrt(norm);
    //norm = 1.0e-6;
    //printf("%e\t",norm);
    norm *= 1.0e6;
    //printf("%e\n",norm);
    for (i=0; i<n; i++) {
      //printf("%d - %e\t",i,Ith(xdot,i));
      Ith(xdot,i) /= norm;
      //printf("%e\n",Ith(xdot,i));
    }
  }
  return CV_SUCCESS;
}

int InterpSystem::Events (realtype t, N_Vector x, realtype * event, void * data) {
  /*RHS(t,x,xderiv,data);
  for(int i=0; i<GetNumberOfEvents(); i++)
    event[i] = Ith(xderiv,i);
  */return CV_SUCCESS;
}

bool InterpSystem::SpecialOptions(const void *opt) {
  bool *b = (bool*) opt;
  backward = b[0];
  arclength = b[1];
  return true;
}

} // namespace bal

