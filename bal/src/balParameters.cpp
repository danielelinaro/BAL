/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balParameters.cpp
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

#include "balParameters.h"

balParameters::balParameters(const balParameters& params) {
	p = 0;
	pars = NULL;
	SetNumber(params.GetNumber());
	for(int i=0; i<p; i++)
		pars[i] = params.pars[i];
}

void balParameters::SetNumber (int numpars) {
	if(numpars > 0) {
		p = numpars;
		if(pars != NULL)
			delete pars;
		pars = new double[p];
		for(int i=0; i<p; i++) pars[i] = 0.0;
	}
}

int balParameters::GetNumber () const {
	return p;
}

double & balParameters::At (int k) {
	return pars[k];
}

double * balParameters::GetParameters () const {
	return pars;
}

void balParameters::CopyValues(balParameters* _par){
	for (int i = 0; i < _par->GetNumber(); i++)
		pars[i] = _par->At(i);
}

