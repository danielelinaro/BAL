/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balDynasty.cpp
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
 * \file balDynasty.cpp
 * \brief Implementation of the class Dynasty
 */

#include "balDynasty.h"
#include <cmath>

bal::DynamicalSystem* DynastyFactory() {
  return new bal::Dynasty;
}

namespace bal {

Dynasty::Dynasty() : DynamicalSystem() {
  SetDimension(3);
  SetNumberOfParameters(7);
  //SetNumberOfEvents(GetDimension());
  SetNumberOfEvents(1);
  constraint_type = MINIMA;
  xderiv = N_VNew_Serial(GetDimension());
}

Dynasty::Dynasty(const Dynasty& dyn) : DynamicalSystem(dyn) {
  constraint_type = dyn.constraint_type;
  xderiv = N_VNew_Serial(dyn.GetDimension());
  for(int i = 0; i < dyn.GetDimension(); i++)
    Ith(xderiv,i)=Ith(dyn.xderiv,i);
}

Dynasty::~Dynasty() {
  N_VDestroy_Serial(xderiv);
}

int Dynasty::RHS (realtype t, N_Vector x, N_Vector xdot, void * data) {
  realtype x1, x2, x3;
  realtype r, e, b, d, g, h, q;
  Parameters *parameters;
  
  parameters = (Parameters*) data;
  r = exp(parameters->At(0));
  e = parameters->At(1);
  b = parameters->At(2);
  d = parameters->At(3);
  g = parameters->At(4);
  h = parameters->At(5);
  q = parameters->At(6);

  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);

  if(x1 < 0.0 || x2 < 0.0 || x3 < 0.0) {
    // small negative components are trimmed to zero
#ifdef NEGATIVES
    fprintf(stderr, "x = [%.12f,%.12f,%.12f]\n", x1, x2, x3);
#endif
    if(x1 > -EPS)
      x1 = 0.0;
    else
      return CV_ERR_FAILURE;
    if(x2 > -EPS)
      x2 = 0.0;
    else
      return CV_ERR_FAILURE;
    if(x3 > -EPS)
      x3 = 0.0;
    else
      return CV_ERR_FAILURE;
  }

  Ith (xdot, 0) = x1*(1-x1-x2/(b+x1)-h*x3);
  Ith (xdot, 1) = q*x2*(e*x1/(b+x1)-1-x3/(d+x2));
  Ith (xdot, 2) = r*(x1*x2/(b+x1)-g*x3);

#ifdef DEBUG
  fprintf(stderr, "   x = [%f,%f,%f]\n",Ith(x,0),Ith(x,1),Ith(x,2));
  fprintf(stderr, "xdot = [%f,%f,%f]\n",Ith(xdot,0),Ith(xdot,1),Ith(xdot,2));
#endif

  return CV_SUCCESS;
}

#ifdef CVODE25
int Dynasty::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
			  void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int Dynasty::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
			  void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  realtype x1, x2, x3;
  realtype r, e, b, d, g, h, q;
  realtype denx, deny, denxsquare, denysquare;
  Parameters *parameters;
  
  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
 
  parameters = (Parameters*) jac_data;
  r = exp(parameters->At(0));
  e = parameters->At(1);
  b = parameters->At(2);
  d = parameters->At(3);
  g = parameters->At(4);
  h = parameters->At(5);
  q = parameters->At(6);

  denx = 1/(b+x1);
  deny = 1/(d+x2);
  denxsquare = denx*denx;
  denysquare = deny*deny;

  IJth (J, 0, 0) = 1 - 2*x1 - b*x2*denxsquare - h*x3;
  IJth (J, 0, 1) = - x1*denx;
  IJth (J, 0, 2) = - h*x1;
  IJth (J, 1, 0) = q*e*b*x2*denxsquare;
  IJth (J, 1, 1) = q*e*x1*denx - q*d*x3*denysquare - q;
  IJth (J, 1, 2) = -q*x2*deny;
  IJth (J, 2, 0) = r*b*x2*denxsquare;
  IJth (J, 2, 1) = r*x1*denx;
  IJth (J, 2, 2) = -r*g;

  return CV_SUCCESS;
}

int Dynasty::Events (realtype t, N_Vector x, realtype * event, void * data) {
  RHS(t,x,xderiv,data);
  for(int i=0; i<GetNumberOfEvents(); i++)
    event[i] = Ith(xderiv,i);
  return CV_SUCCESS;
}

void Dynasty::EventsConstraints (realtype t, N_Vector x, int * constraints, void * data) {
  realtype x1, x2, x3;
  realtype r, e, b, d, g, h, q;
  realtype denx, deny, denxsquare, denysquare;
  realtype ris[3], xdot[3];
  Parameters *parameters;
  
  if(constraint_type != MINIMA && constraint_type != MAXIMA) {
    for(int i=0; i<GetNumberOfEvents(); i++)
      constraints[i] = 1;
    return;
  }

  x1 = Ith (x, 0);
  x2 = Ith (x, 1);
  x3 = Ith (x, 2);
 
  parameters = (Parameters*) data;
  r = exp(parameters->At(0));
  e = parameters->At(1);
  b = parameters->At(2);
  d = parameters->At(3);
  g = parameters->At(4);
  h = parameters->At(5);
  q = parameters->At(6);

  denx = 1/(b+x1);
  deny = 1/(d+x2);
  denxsquare = denx*denx;
  denysquare = deny*deny;

  RHS(t,x,xderiv,data);
  for(int i=0; i<GetDimension(); i++)
    xdot[i] = Ith(xderiv,i);

  /*
  ris[0] = xdot[0] - 2*x1*xdot[0] - h*x3*xdot[0] - h*x1*xdot[2] - b*x2*denxsquare*xdot[0] - x1*denx*xdot[1];
  ris[1] = q*e*b*x2*denxsquare*xdot[0] + q*e*x1*denx*xdot[1] - q*d*x3*denysquare*xdot[1] - q*x2*deny*xdot[2];
  ris[2] = r*b*x2*denxsquare*xdot[0] + r*x1*denx*xdot[1] - r*g*xdot[2];
  */

  ris[0] = (1-x1-x2*denx-h*x3)*xdot[0]+x1*(-xdot[0]+x2*xdot[0]*denxsquare-xdot[1]*denx-h*xdot[2]);
  ris[1] = q*(-1+e*x1*denx-x3*deny)*xdot[1]+q*xdot[1]*(-e*x1*xdot[0]*denxsquare+e*xdot[0]*denx+x3*xdot[1]*denysquare-xdot[2]*deny);
  ris[2] = r*(-x1*x2*xdot[0]*denxsquare+x2*xdot[0]*denx+x1*xdot[1]*denx-g*xdot[2]);
	
  switch(constraint_type) {
  case MINIMA:
    for(int i=0; i<GetNumberOfEvents(); i++)
      constraints[i] = (ris[i] > 0 ? 1 : 0);
    break;
  case MAXIMA:
    for(int i=0; i<GetNumberOfEvents(); i++)
      constraints[i] = (ris[i] < 0 ? 1 : 0);
    break;
  default:
    for(int i=0; i<GetNumberOfEvents(); i++)
      constraints[i] = 1;
  }
}

bool Dynasty::SpecialOptions(const void *opt) {
  char *s = (char *) opt;
  if(strncasecmp(s,"minima",6) == 0) {
    constraint_type = MINIMA;
    return true;
  }
  if(strncasecmp(s,"maxima",6) == 0) {
    constraint_type = MAXIMA;
    return true;
  }
  if(strncasecmp(s,"any",3) == 0) {
    constraint_type = ANY;
    return true;
  }
  return false;
}

bool Dynasty::HasJacobian() const {
  return true;
}
 
bool Dynasty::HasEvents() const {
  return true;
}
 
bool Dynasty::HasEventsConstraints() const {
  return true;
}

DynamicalSystem* Dynasty::Clone() const {
  return new Dynasty(*this);
}

} // namespace bal

