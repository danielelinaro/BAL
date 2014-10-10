/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balHeartNeuron.cpp
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
 * \file balHeartNeuron.cpp
 * \brief Implementation of the class HeartNeuron
 */

#include <cmath>
#include "balHeartNeuron.h"

bal::DynamicalSystem* HeartNeuronFactory() {
  return new bal::HeartNeuron;
}

namespace bal {

HeartNeuron::HeartNeuron() : DynamicalSystem(3,3,1,false),
  C(0.5),gK2(30.),EK(-0.07),ENa(0.045),gNa(200.),E1(-0.046),g1(8.),tauNa(0.0405022) {
  A[0] = -150.;
  A[1] = 500.;
  A[2] = -83.;
  B[0] = 0.0305;
  B[1] = 0.0333;
  B[2] = 0.018;
  xderiv = N_VNew_Serial(GetDimension());
}

HeartNeuron::HeartNeuron(const HeartNeuron& hn) : DynamicalSystem(hn),
  C(hn.C),gK2(hn.gK2),EK(hn.EK),ENa(hn.ENa),gNa(hn.gNa),E1(hn.E1),g1(hn.g1),tauNa(hn.tauNa) {
  A[0] = hn.A[0];
  A[1] = hn.A[1];
  A[2] = hn.A[2];
  B[0] = hn.B[0];
  B[1] = hn.B[1];
  B[2] = hn.B[2];
  xderiv = N_VNew_Serial(GetDimension());
  for(int i = 0; i < GetDimension(); i++)
    Ith(xderiv,i) = Ith(hn.xderiv,i);
}

HeartNeuron::~HeartNeuron(){
  N_VDestroy_Serial(xderiv);
}

bool HeartNeuron::HasJacobian() const {
  return false;
}

bool HeartNeuron::HasEvents() const {
  return true;
}

bool HeartNeuron::HasEventsConstraints() const {
  return true; 
}
  
int HeartNeuron::RHS(realtype t, N_Vector x, N_Vector xdot, void *sys) {
  realtype VK2shift,Iapp,tauK2;
  realtype V,hNa,mK2;
  realtype f;
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  Parameters *parameters = ds->GetParameters();

  VK2shift = parameters->At(0);
  Iapp = parameters->At(1);
  tauK2 = parameters->At(2);
  
  V = Ith(x,0);
  hNa = Ith(x,1);
  mK2 = Ith(x,2);
  
  f = BoltzmannF(A[0],B[0],V);
  
  Ith(xdot,0) = -1.0/C * (gK2*mK2*mK2*(V-EK) + g1*(V-E1) + gNa*f*f*f*hNa*(V-ENa) - Iapp);
  Ith(xdot,1) = (BoltzmannF(A[1],B[1],V) - hNa) / tauNa;
  Ith(xdot,2) = (BoltzmannF(A[2],B[2]+VK2shift,V) - mK2) / tauK2; 
  
  return CV_SUCCESS;
}

int HeartNeuron::Events(realtype t, N_Vector x, realtype * event, void * data){
  realtype Iapp;
  realtype V,hNa,mK2;
  realtype f;
  
  Parameters *param;
  param = (Parameters *) data;
  Iapp = param->At(1);
  
  V = Ith(x,0);
  hNa = Ith(x,1);
  mK2 = Ith(x,2);
  
  f = BoltzmannF(A[0],B[0],V);
  
  event[0] = -1.0/C * (gK2*mK2*mK2*(V-EK) + g1*(V-E1) + gNa*f*f*f*hNa*(V-ENa) - Iapp);
  return CV_SUCCESS;
}

#ifdef CVODE25
int HeartNeuron::Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
			      void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
#ifdef CVODE26
int HeartNeuron::Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
			      void *sys, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
#endif
  realtype VK2shift,Iapp,tauK2;
  realtype V,hNa,mK2;
  realtype f;
  DynamicalSystem *ds = static_cast<DynamicalSystem*>(sys);
  Parameters *parameters = ds->GetParameters();
  
  VK2shift = parameters->At(0);
  Iapp = parameters->At(1);	
  tauK2 = parameters->At(2);	
  
  V = Ith(x,0);
  hNa = Ith(x,1);
  mK2 = Ith(x,2);
  
  f = BoltzmannF(A[0],B[0],V);
  
  IJth(J,0,0) = -(1.0/C) * (gK2*mK2*mK2 + g1 + gNa*f*f*f*hNa + 3.0*gNa*f*f*hNa*(V-ENa)*BoltzmannDFDV(A[0],B[0],V));
  IJth(J,0,1) = -(1.0/C) * (gNa*f*f*f*(V-ENa));
  IJth(J,0,2) = -(2.0/C) * (gK2*mK2*(V-EK));
  IJth(J,1,0) = BoltzmannDFDV(A[1],B[1],V)/tauNa;
  IJth(J,1,1) = -1.0/tauNa;
  IJth(J,1,2) = 0.0;
  IJth(J,2,0) = BoltzmannDFDV(A[2],B[2]+VK2shift,V)/tauK2;
  IJth(J,2,1) = 0.0;
  IJth(J,2,2) = -1.0/tauK2;

  return CV_SUCCESS;	
}

void HeartNeuron::EventsConstraints (realtype t, N_Vector x, int * constraints, void * data){
  realtype VK2shift,Iapp,tauK2;
  realtype V,hNa,mK2;
  
  realtype xdot[3];
  realtype ris;
  realtype f,dfdv;
  
  Parameters *param;
  param = (Parameters *) data;
  VK2shift = param->At(0);
  Iapp = param->At(1);
  tauK2 = param->At(1);
  
  V = Ith(x,0);
  hNa = Ith(x,1);
  mK2 = Ith(x,2);
  
  RHS(t,x,xderiv,data);       
  for(int i=0; i<GetDimension(); i++) xdot[i] = Ith(xderiv,i);   
  
  f = BoltzmannF(A[0],B[0],V);
  dfdv = BoltzmannDFDV(A[0],B[0],V);
  
  ris = -1.0/C * (2*gK2*mK2*xdot[2]*(V-EK) + gK2*mK2*mK2*xdot[0] + g1*xdot[0] +
		  gNa*f*f*f*xdot[1]*(V-ENa) + gNa*f*f*f*hNa*xdot[0] + 3*gNa*f*f*dfdv*hNa*xdot[0]*(V-ENa));
  
  if (ris < 0) 
    constraints[0] = 1;
  else
    constraints[0] = 0;
}

DynamicalSystem* HeartNeuron::Clone() const {
  return new HeartNeuron(*this);
}

} // namespace bal

