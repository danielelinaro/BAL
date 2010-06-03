/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    main.cpp
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

#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <nvector/nvector_serial.h>
#include "balCommon.h"
#include "balObject.h"
#include "balDynamicalSystem.h"
#include "balHindmarshRose.h"
#include "balParameters.h"
#include "balSolution.h"
#include "balLogger.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"
#include "balBifurcationParameters.h"
#include "balHeartNeuron.h"
#include "balPLL.h"
#include "balInterp1D.h"

using namespace std;

void testInterp();
void bifd(double p1l, double p1u, double p2l, double p2u, int steps[2], realtype x0[3], const char * filename);

int main(int argc, char *argv[]) {
	// TEST balObject
	/*
	balObject *obj = balObject::Create();
	cout << obj->GetClassName() << endl;
	cout << obj->IsA("balObject") << endl;
	obj->Destroy();
	*/

	/*
	// TEST balParameters
	balParameters * pars = balParameters::Create();
	cout << pars->GetClassName() << endl;
	pars->SetNumber(4);
	pars->At(0) = 3.0;
	pars->At(1) = 5.0;
	pars->At(2) = 0.01;
	pars->At(3) = 4.0;
	cout << *pars << endl;
	balParameters * parsCopy = balParameters::Copy(pars);
	cout << *parsCopy << endl;
	*/

	/*
	double buffer[10] = {0, 0, 0, 0, -2, 1, 1, 1, 1, -1}; 
	// TEST balH5Logger
	balLogger * logger = balH5Logger::Create();
	logger->SetFilename("test.1.h5");
	logger->SetNumberOfColumns(5);
	logger->SetParameters(pars);
	logger->SaveBuffer(buffer, 2);
	logger->SetFilename("test.2.h5");
	logger->SaveBuffer(buffer, 2);
	logger->Destroy();
	// TEST balASCIILogger
	logger = balASCIILogger::Create();
	logger->SetFilename("test.1.dat");
	logger->SetNumberOfColumns(5);
	logger->SetParameters(pars);
	logger->SaveBuffer(buffer, 2);
	logger->SetFilename("test.2.dat");
	logger->SaveBuffer(buffer, 2);
	logger->Destroy();
	*/
	
	// TEST balDynamicalSystem e balHindmarshRose
	/*
	int i;
	balDynamicalSystem * dynsys = balDynamicalSystem::Create();
	balDynamicalSystem * hr = balHindmarshRose::Create();
	N_Vector x = N_VNew_Serial(hr->GetDimension());
	N_Vector xdot = N_VNew_Serial(hr->GetDimension());
	for(i=0; i<hr->GetDimension(); i++) NV_Ith_S(x,i) = 0.0;
	cout << dynsys->GetClassName() << endl;
	cout << hr->GetClassName() << endl;
	dynsys->SetParameters(pars);
	hr->SetParameters(pars);
	balDynamicalSystem::RHSWrapper(0,NULL,NULL,dynsys);
	balDynamicalSystem::RHSWrapper(0,x,xdot,hr);

	cout << "xdot = (";
	for(i=0; i<hr->GetDimension()-1; i++) cout << NV_Ith_S(xdot,i) << ",";
	cout << NV_Ith_S(xdot,i) << ")" << endl;

	N_VDestroy_Serial(x);
	N_VDestroy_Serial(xdot);
	dynsys->Destroy();
	*/

	// TEST balBifurcationParameters
	/*
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

	cout << endl;
	while(bp->HasNext()) {
		cout << *bp << endl;
		bp->Next();
	}
	cout << endl;
	bp->Reset();
	while(bp->HasNext()) {
		cout << *bp << endl;
		bp->Next();
	}
	cout << endl;
	pars->Destroy();
	parupper->Destroy();
	bp->Reset();
	*/

	// TEST balODESolver
	/*
	balODESolver * solver = balODESolver::Create();
	hr->SetParameters(bp);
	solver->SetDynamicalSystem(hr);
	solver->SetTransientDuration(0.0);
	solver->HaltAtEquilibrium(true);
	solver->SetIntegrationMode(balTRAJ);
	solver->Solve();
	solver->SetIntegrationMode(balEVENTS);
	solver->Solve();
	solver->SetIntegrationMode(balBOTH);
	solver->Solve();
	pars->At(1) = 1.0;
	solver->Solve();
	solver->Destroy();
	*/

	/**/
	// TEST balBifurcationDiagram
	int steps[4] = {1,503,1,1};
	realtype x0[3] = {0.5,0.5,0.5};
	balBifurcationParameters * bp = balBifurcationParameters::Create();
	bp->SetNumber(4);
	bp->SetIthParameter(0,3.0);
	bp->SetIthParameter(2,0.01);
	bp->SetIthParameter(3,4.0);
	bp->SetIthParameterLowerBound(1,2.5);
	bp->SetIthParameterUpperBound(1,4.5);
	bp->SetNumberOfSteps(steps);
	balHindmarshRose * hr = balHindmarshRose::Create();
	hr->SetParameters(bp);
	balBifurcationDiagram * bifd = balBifurcationDiagram::Create();
	bifd->SetDynamicalSystem(hr);
	bifd->SetFilename(argv[1]);
	bifd->GetODESolver()->SetIntegrationMode(balBOTH);
	bifd->GetODESolver()->SetTransientDuration(0);
	bifd->GetODESolver()->HaltAtEquilibrium(true);
	bifd->GetODESolver()->SetFinalTime(1e6);
	bifd->GetODESolver()->SetMaxNumberOfIntersections(100);
	bifd->GetODESolver()->SetX0(x0);

	cout << "Computing the bifurcation diagram...\n";
	flush(cout);
//	clock_t start = clock();
	bifd->SetNumberOfThreads(4);
	
	bool m;
	char s;
	printf("MULTITHREADED? (y or else)\t");
	scanf("%c",&s);
	
	m = ((s =='y')?true:false);
	
	bifd->ComputeDiagram(m);
//	clock_t end = clock();
	cout << "... done!\n";
	bifd->Destroy();
	hr->Destroy();
	bp->Destroy();
	//printf("Time elapsed: %f\n", ((double) end - start) / CLOCKS_PER_SEC);
	
	
	
	
	
	/**/

	// TEST balHeartNeuron
	/*
	int steps[2] = {1,11};
	realtype x0[3] = {0.,.4,.4};
	balBifurcationParameters * bp = balBifurcationParameters::Create();
	bp->SetNumber(3);
	bp->SetIthParameterLowerBound(0,-0.022);
	bp->SetIthParameterUpperBound(0,-0.022);
	bp->SetIthParameterLowerBound(1,0.016);
	bp->SetIthParameterUpperBound(1,0.018);
	bp->SetIthParameterLowerBound(2,0.25);
	bp->SetIthParameterUpperBound(2,0.25);
	bp->SetNumberOfSteps(steps);
	balHeartNeuron * hn = balHeartNeuron::Create();
	hn->SetParameters(bp);
	balBifurcationDiagram * bifd = balBifurcationDiagram::Create();
	bifd->RestartFromX0(false);
	bifd->SetDynamicalSystem(hn);
	bifd->SetFilename("HeartNeuron2D.h5");
	bifd->GetODESolver()->SetIntegrationMode(balBOTH);
	bifd->GetODESolver()->SetTransientDuration(250);
	bifd->GetODESolver()->HaltAtEquilibrium(true);
	bifd->GetODESolver()->SetX0(x0);
	//bifd->GetODESolver()->HaltAtCycle(true);
	bifd->GetODESolver()->SetFinalTime(25000);
	bifd->GetODESolver()->SetTimeStep(5e-3);
	bifd->GetODESolver()->SetMaxNumberOfIntersections(1000);
	bifd->ComputeDiagram();
	bifd->Destroy();
	hn->Destroy();
	bp->Destroy();
	*/
	/*
	int steps[2] = {200,200};
	realtype x0[3] = {0.,.4,.4};
	if(fork()) {
		if(fork()) {
			bifd(-0.025,-0.01,0.,0.02,steps,x0,"HeartNeuron2Da.h5");
		}
		else {
			bifd(-0.01,0.005,0.,0.02,steps,x0,"HeartNeuron2Db.h5");
		}
	}
	else {
		if(fork()) {
			bifd(-0.023,-0.021,0.015,0.017,steps,x0,"HeartNeuron2Dc.h5");
		}
		else {
			bifd(-0.023,-0.021,0.017,0.019,steps,x0,"HeartNeuron2Dd.h5");
		}
	}
	*/

	//testInterp();

	return EXIT_SUCCESS;
}

void testInterp() {
	balBaseInterp1D * interp;
	double x[10], y[10];
	for(int i=0; i<10; i++) {
		x[i] = i;
		y[i] = i*i;
	}
	//interp = balLinearInterp1D::Create(x,y,10);
	//interp = balPolyInterp1D::Create(x,y,10,3);
	interp = balSplineInterp1D::Create(x,y,10);
	for(double xx=0.1; xx<=9.9; xx+=0.1) {
		fprintf(stdout, "%e %e\n", xx, interp->interp(xx));
	}
	interp->Destroy();
}

void bifd(double p1l, double p1u, double p2l, double p2u, int steps[2], realtype x0[3], const char * filename) {
	balBifurcationParameters * bp = balBifurcationParameters::Create();
	bp->SetNumber(3);
	bp->SetIthParameterLowerBound(0,p1l);
	bp->SetIthParameterUpperBound(0,p1u);
	bp->SetIthParameterLowerBound(1,p2l);
	bp->SetIthParameterUpperBound(1,p2u);
	bp->SetIthParameterLowerBound(2,0.25);
	bp->SetIthParameterUpperBound(2,0.25);
	bp->SetNumberOfSteps(steps);
	balHeartNeuron * hn = balHeartNeuron::Create();
	hn->SetParameters(bp);
	balBifurcationDiagram * bifd = balBifurcationDiagram::Create();
	bifd->SetDynamicalSystem(hn);
	bifd->SetFilename(filename);
	bifd->GetODESolver()->SetIntegrationMode(balEVENTS);
	bifd->GetODESolver()->SetTransientDuration(200);
	bifd->GetODESolver()->HaltAtEquilibrium(true);
	//bifd->GetODESolver()->HaltAtCycle(true);
	bifd->GetODESolver()->SetX0(x0);
	bifd->GetODESolver()->SetFinalTime(10e6);
	bifd->GetODESolver()->SetMaxNumberOfIntersections(200);
	bifd->ComputeDiagram();
	bifd->Destroy();
	hn->Destroy();
	bp->Destroy();
}


