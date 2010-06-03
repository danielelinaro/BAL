/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    dynsys.cpp
 *
 *   Copyright (C) 2009 Daniele Linaro
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

#include <iostream>
#include <nvector/nvector_serial.h>
#include "balObject.h"
#include "balDynamicalSystem.h"
#include "balHindmarshRose.h"
#include "balParameters.h"
using namespace std;

// TEST balDynamicalSystem and balHindmarshRose
int main(int argc, char *argv[]) {
	
	// parameters
	balParameters * pars = balParameters::Create();
	pars->SetNumber(4);
	pars->At(0) = 3.0;
	pars->At(1) = 5.0;
	pars->At(2) = 0.01;
	pars->At(3) = 4.0;

	// balDynamicalSystem
	balDynamicalSystem * dynsys = balDynamicalSystem::Create();
	cout << dynsys->GetClassName() << endl;
	dynsys->Destroy();

	// balHindmarshRose
	int i;
	balDynamicalSystem * hr = balHindmarshRose::Create();
	N_Vector x = N_VNew_Serial(hr->GetDimension());
	N_Vector xdot = N_VNew_Serial(hr->GetDimension());
	for(i=0; i<hr->GetDimension(); i++) NV_Ith_S(x,i) = 0.0;
	cout << hr->GetClassName() << endl;
	hr->SetParameters(pars);
	balDynamicalSystem::RHSWrapper(0,x,xdot,hr);

	cout << "xdot = (";
	for(i=0; i<hr->GetDimension()-1; i++) cout << NV_Ith_S(xdot,i) << ",";
	cout << NV_Ith_S(xdot,i) << ")" << endl;

	N_VDestroy_Serial(x);
	N_VDestroy_Serial(xdot);
	hr->Destroy();

	return 0;
}

