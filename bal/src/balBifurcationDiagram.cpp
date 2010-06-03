/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balBifurcationDiagram.cpp
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

#include "balBifurcationDiagram.h"

void ResetColours(int d) {
	printf("%c%s", ESC, RED);
	printf("Aborting...\n");
	printf("%c%s", ESC, NORMAL);
	// !!!!! NEED TO ADD MEMORY MANAGEMENT (I.E. FREEING USED MEMORY BEFORE
	// QUITTING) !!!!!
	exit(1);
}

balBifurcationDiagram::balBifurcationDiagram() {
	logger = balH5Logger::Create();
	destroy_logger = true;
	solver = balODESolver::Create();
	destroy_solver = true;
	restart_from_x0 = true;
	nthreads = 2;
	classification = NULL;
	SetFilename("balBifurcationDiagram.h5");
	signal(SIGINT, ResetColours);
}

bool balBifurcationDiagram::SetNumberOfThreads(int _nthreads) {
	if(_nthreads <= 0) {
		return false;
	}
	nthreads = _nthreads;
	return true;
}

balBifurcationDiagram::~balBifurcationDiagram() {
	if(destroy_logger)
		logger->Destroy();
	if(destroy_solver)
		solver->Destroy();
	if(classification != NULL)
		delete classification;
}

void balBifurcationDiagram::SetDynamicalSystem(balDynamicalSystem * sys) {
	system = sys;
	solver->SetDynamicalSystem(sys);
	parameters = system->GetParameters();
	npar = parameters->GetNumber();
}

balDynamicalSystem * balBifurcationDiagram::GetDynamicalSystem() const {
	return system;
}

void balBifurcationDiagram::SetLogger(balLogger * log) {
	if(destroy_logger) {
		logger->Destroy();
		destroy_logger = false;
	}
	logger = log;
}

balLogger * balBifurcationDiagram::GetLogger() const {
	return logger;
}

void balBifurcationDiagram::SetODESolver(balODESolver * sol) {
	if(destroy_solver) {
		solver->Destroy();
		destroy_solver = false;
	}
	solver = sol;
}

balODESolver * balBifurcationDiagram::GetODESolver() const {
	return solver;
}

void balBifurcationDiagram::SetFilename(const char * filename) {
	logger->SetFilename(filename);
}

const char * balBifurcationDiagram::GetFilename() {
	return logger->GetFilename();
}


bool balBifurcationDiagram::SaveClassificationData(const char *filename) const {
	if(classification == NULL)
		return false;

	//classification->sort(balDoubleArrayComparer());

	int i;
	list<double *>::iterator it;
	FILE *fid; 

	fid = fopen(filename,"w"); 
	if(fid == NULL)
		return false;

	for(it=classification->begin(); it!=classification->end(); it++) {
		for(i=0; i<npar; i++)
			fprintf(fid, "%e ", (*it)[i]);
		fprintf(fid,"%d\n", (int) (*it)[npar]);
	}

	fclose(fid);
	return true;
}

double** balBifurcationDiagram::GetClassificationData() const {
	if(classification == NULL)
		return NULL;

	//classification->sort(balDoubleArrayComparer());

	int i, j;
	list<double *>::iterator it;
	double **data;
	data = new double*[classification->size()];
	for(it=classification->begin(), i=0; it!=classification->end(); it++, i++) {
		data[i] = new double[npar+1];
		for(j=0; j<npar+1; j++)
			data[i][j] = (*it)[j];
	}
}

void balBifurcationDiagram::ComputeDiagram() {
	if(classification != NULL)
		delete classification;
	classification = new list<double *>;
	if(nthreads == 1)
		ComputeDiagramSingleThreaded();
	else
		ComputeDiagramMultiThreaded();
}

double* balBifurcationDiagram::BuildClassificationEntry(balSolution *sol) {
	balParameters *p = sol->GetParameters();
	double *entry = new double[npar+1];
	for(int i=0; i<npar; i++)
		entry[i] = p->At(i);
	entry[npar] = sol->GetNumberOfTurns();
	return entry;
}

void balBifurcationDiagram::ComputeDiagramSingleThreaded() {
	balBifurcationParameters * pars = (balBifurcationParameters *) parameters;
	pars->Reset();
	int total = pars->GetTotalNumberOfTuples();
	balSolution * solution;
	for(int cnt=0; cnt<total; cnt++, pars->Next()) {
		printf("%c%s", ESC, GREEN);
		printf("[%05d/%05d]\r", cnt+1, total); fflush(stdout);
		printf("%c%s", ESC, NORMAL);
		solver->Solve();
		solution = solver->GetSolution();
		classification->push_back(BuildClassificationEntry(solution));
		logger->SaveSolution(solution);
		solution->Destroy();
		if(! restart_from_x0) {
			solver->SetX0(solver->GetXEnd());
		}
	}
	printf("\n"); fflush(stdout);
}


void balBifurcationDiagram::ComputeDiagramMultiThreaded() {
	int i, cnt;
	balBifurcationParameters * pars = (balBifurcationParameters *) parameters;
	pars->Reset();
	
	int total = pars->GetTotalNumberOfTuples();
	int prologue = total%nthreads;
	boost::thread * threads[nthreads];
	balODESolver** threaded_solvers = new balODESolver*[nthreads];
	balParameters** threaded_parameters = new balParameters*[nthreads];
	
	/* creating #nthreads indipendent solvers with indipendent parameters array */
	for(i = 0; i < nthreads; i++) {
		threaded_solvers[i] = balODESolver::Copy(solver);
		threaded_parameters[i] = balParameters::Copy(parameters);
		threaded_solvers[i]->SetDynamicalSystemParameters(threaded_parameters[i]);
	}

	/* launching writing routine threaded. it activates only when list size >= LIST_MAX_SIZE */
	logger_thread = new boost::thread(&balLogger::SaveSolutionThreaded,logger,&solution_list,&list_mutex,&q_empty,&q_full);

	/* calculating the first #total%nthreads solutions unthreaded */
	for(cnt = 0; cnt < prologue; cnt++, pars->Next()) {
		IntegrateAndEnqueue(solver);
		printf("%c%s", ESC, GREEN);
		printf("[%05d/%05d]\r", cnt+1, total); fflush(stdout);
		printf("%c%s", ESC, NORMAL);
	}
	
	for (cnt=0; cnt < total/nthreads; cnt++) {
		/* setting to each threaded solver the correct parameters array taken from pars->Next() */
		for(i = 0; i < nthreads; i++) {
			threaded_parameters[i]->CopyValues(pars);
			pars->Next();
		}
		
		/* launching the nthread solvers */
		for (i = 0; i < nthreads; i++)
			threads[i] = new boost::thread(&balBifurcationDiagram::IntegrateAndEnqueue,this,threaded_solvers[i]);
					
		for (i = 0; i < nthreads; i++) {
			(*(threads[i])).join();
			printf("%c%s", ESC, GREEN);
			printf("[%05d/%05d]\r", (prologue + cnt*nthreads + i + 1), total); fflush(stdout);
			printf("%c%s", ESC, NORMAL);
		}
	}
	
	printf("\n"); fflush(stdout);
	
	for(i = 0; i < nthreads; i++) {
		threaded_solvers[i]->Destroy();
		threaded_parameters[i]->Destroy();
	}
	
	/* interrupting logger_thread */
	(*logger_thread).interrupt();
	
	/* waiting logger thread to writing the remaining solution in the queue and exit */
	(*logger_thread).join();
}

void balBifurcationDiagram::IntegrateAndEnqueue(balODESolver * sol) {
	sol->Solve();
	balSolution *solution = sol->GetSolution();
	classification->push_back(BuildClassificationEntry(solution));
	
	{
		boost::mutex::scoped_lock lock(list_mutex);
		while (solution_list.size() >= LIST_MAX_SIZE) 
		{  
			/* notifies the logger_thread (waiting on q_full) that the solution list in now full */
		   q_full.notify_one();
			/* the thread goes on wait on q_empty cond var, unlocking list_mutex until logger_thread has emptied the list */
		   q_empty.wait(lock);
		}
		solution_list.push_back(solution);
	}	
	
	if(! restart_from_x0) {
		sol->SetX0(sol->GetXEnd());
	}
}

