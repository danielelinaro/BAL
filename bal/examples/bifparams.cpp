/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    bifparams.cpp
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
#include "balParameters.h"
#include "balBifurcationParameters.h"
using namespace std;

// TEST balBifurcationParameters
int main(int argc, char *argv[]) {

	// parameters
	balParameters * pars = balParameters::Create();
	pars->SetNumber(4);
	pars->At(0) = 3.0;
	pars->At(1) = 5.0;
	pars->At(2) = 0.01;
	pars->At(3) = 4.0;

	balBifurcationParameters * bp = balBifurcationParameters::Create();
	balParameters * parupper = balParameters::Create();
	parupper->SetNumber(4);
	parupper->At(0) = pars->At(0)+1;
	parupper->At(1) = pars->At(1)+1;
	parupper->At(2) = pars->At(2);
	parupper->At(3) = pars->At(3);
	bp->SetParameterBounds(pars, parupper);
	int steps[4] = {6,6,1,1};
	bp->SetNumberOfSteps(steps);

	cout << "par lower: " << *pars << endl;
	cout << "par upper: " << *parupper << endl;

	// print all tuples
	while(bp->HasTuples()) {
		cout << *bp << endl;
		bp->Next();
	}
	cout << endl;
	bp->Reset();
	// leave last tuple out
	while(bp->HasNext()) {
		cout << *bp << endl;
		bp->Next();
	}
	cout << endl;
	pars->Destroy();
	parupper->Destroy();
	bp->Reset();

	return 0;
}

