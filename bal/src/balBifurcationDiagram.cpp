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

/** 
 * \file balBifurcationDiagram.cpp
 * \brief Implementation of the class balBifurcationDiagram
 */

#include "balBifurcationDiagram.h"

void ResetColours(int d) {
  printf("%c%s", ESC, RED);
  printf("Aborting...\n");
  printf("%c%s", ESC, NORMAL);
  // !!!!! NEED TO ADD MEMORY MANAGEMENT (I.E. FREEING USED MEMORY BEFORE QUITTING) !!!!!
  exit(1);
}

bool CompareBalClassificationEntries(balClassificationEntry *entry1, balClassificationEntry *entry2) {
  return entry1->GetID() < entry2->GetID();
}

/***** balClassificationEntry *****/

balClassificationEntry::balClassificationEntry(balSolution *sol) {
  n = sol->GetParameters()->GetNumber() + 1;
  data = new double[n];
  id = sol->GetID();
  for(int i=0; i<n-1; i++)
    data[i] = sol->GetParameters()->At(i);
  data[n-1] = sol->GetNumberOfTurns();
}

int balClassificationEntry::GetN() const {
  return n;
}

int balClassificationEntry::GetID() const {
  return id;
}

double* balClassificationEntry::GetData() const {
  return data;
}

/***** balBifurcationDiagram *****/

const char * balBifurcationDiagram::GetClassName() const {
  return "balBifurcationDiagram";
}

balBifurcationDiagram * balBifurcationDiagram::Create() {
  return new balBifurcationDiagram;
}

void balBifurcationDiagram::Destroy() {
  delete this;
}

balBifurcationDiagram::balBifurcationDiagram() {
  logger = balH5Logger::Create();
  destroy_logger = true;
  solver = balODESolver::Create();
  destroy_solver = true;
  restart_from_x0 = true;
  nthreads = 2;
  mode = balPARAMS;
  classification = NULL;
  destroy_classification = false;
  nX0 = 0;
  X0 = NULL;
  SetFilename("balBifurcationDiagram.h5");
  signal(SIGINT, ResetColours);
}

void balBifurcationDiagram::SetNumberOfThreads(int _nthreads) {
  if(_nthreads > 0)
    nthreads = _nthreads;
}

int balBifurcationDiagram::GetNumberOfThreads() const {
  return nthreads;
}

balBifurcationDiagram::~balBifurcationDiagram() {
  if(destroy_logger)
    logger->Destroy();
  if(destroy_solver)
    solver->Destroy();
  if(destroy_classification)
    delete classification;
}

void balBifurcationDiagram::SetDynamicalSystem(balDynamicalSystem * sys) {
  system = sys;
  ndim = system->GetDimension();
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
  if(!destroy_classification)
    return false;
  
  classification->sort(CompareBalClassificationEntries);
  
  int i;
  double *entry;
  list<balClassificationEntry *>::iterator it;
  FILE *fid;
  
  fid = fopen(filename,"w"); 
  if(fid == NULL)
    return false;
  
  for(it=classification->begin(); it!=classification->end(); it++) {
    entry = (*it)->GetData();
    for(i=0; i<(*it)->GetN()-1; i++)
      fprintf(fid, "%e ", entry[i]);
    fprintf(fid,"%d\n", (int) entry[(*it)->GetN()-1]);
  }
  
  fclose(fid);
  return true;
}

double** balBifurcationDiagram::GetClassificationData() const {
  if(!destroy_classification)
    return NULL;
  
  classification->sort(CompareBalClassificationEntries);
  
  int i, j;
  double **data, *entry;
  list<balClassificationEntry *>::iterator it;

  data = new double*[classification->size()];
  for(it=classification->begin(), i=0; it!=classification->end(); it++, i++) {
    entry = (*it)->GetData();
    data[i] = new double[(*it)->GetN()];
    for(j=0; j<(*it)->GetN(); j++)
      data[i][j] = entry[j];
  }
  return data;
}

void balBifurcationDiagram::SetInitialConditions(int nx0, double **x0) {
  if(nx0 > 0) {
    nX0 = nx0;
    X0 = x0;
  }
}

bool balBifurcationDiagram::RestartsFromX0() const {
  return restart_from_x0;
}

void balBifurcationDiagram::RestartFromX0(bool restart) {
  restart_from_x0 = restart;
}

bool balBifurcationDiagram::SetMode(int _mode) {
  switch(_mode) {
  case balPARAMS:
    break;
  case balIC:
    break;
  default:
    return false;
  }
  mode = _mode;
}

void balBifurcationDiagram::ComputeDiagram() {
  if(destroy_classification)
    delete classification;
  classification = new list<balClassificationEntry *>;
  destroy_classification = true;
  ComputeDiagramMultiThread();
}

void balBifurcationDiagram::ComputeDiagramSingleThread() {
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
    classification->push_back(new balClassificationEntry(solution));
    logger->SaveSolution(solution);
    solution->Destroy();
    if(! restart_from_x0) {
      solver->SetX0(solver->GetXEnd());
    }
  }
  printf("\n"); fflush(stdout);
}

void balBifurcationDiagram::ComputeDiagramMultiThread() {
  int i, idx, cnt, solutionId;
  // the total number of integrations
  int total;
  // the number of integrations that must be performed outside of the main loop
  int prologue;
  // the number of for loops that will be performed after the prologue
  int nloops;

  balBifurcationParameters *pars;

  switch(mode) {
  case balPARAMS:
    pars = (balBifurcationParameters *) parameters;
    pars->Reset();
    total = pars->GetTotalNumberOfTuples();
    break;
  case balIC:
    RestartFromX0(true);
    idx = 0;
    total = nX0;
    break;
  }
  
  nloops = total / nthreads;
  prologue = total % nthreads;

  /**
   * launching the thread of the writing routine: it activates only when list size >= LIST_MAX_SIZE
   **/
  logger_thread = new boost::thread(&balLogger::SaveSolutionThreaded,logger,&solution_list,&list_mutex,&q_empty,&q_full);
  
  /**
   * calculating the first total % nthreads solutions in a serial fashion
   **/
  for(cnt = 0, solutionId = 1; cnt < prologue; cnt++, solutionId++) {
    IntegrateAndEnqueue(solver,solutionId);
    printf("%c%s", ESC, GREEN);
    printf("[%05d/%05d]\r", cnt+1, total);
    printf("%c%s", ESC, NORMAL); fflush(stdout);
    switch(mode) {
    case balPARAMS:
      pars->Next();
      break;
    case balIC:
      idx++;
      break;
    }
  }
  
  /*
   * we make nthreads copies of both the solvers and the parameters of the dynamical system
   */
  balODESolver** lsol = new balODESolver*[nthreads];
  balParameters** lpar = new balParameters*[nthreads];
  for(i = 0; i < nthreads; i++) {
    lpar[i] = balParameters::Copy(parameters);
    lsol[i] = balODESolver::Copy(solver);
    lsol[i]->SetDynamicalSystemParameters(lpar[i]);
  }

  // the array of pointers to the threads
  boost::thread * threads[nthreads];

  for (cnt = 0; cnt < nloops; cnt++) {
    /* we set to each solver a different tuple of parameters */
    switch(mode) {
    case balPARAMS:
      for(i = 0; i < nthreads; i++, pars->Next())
	lpar[i]->CopyValues(pars);
      break;
    case balIC:
      for(i = 0; i < nthreads; i++, idx++)
	lsol[i]->SetX0(X0[idx]);
      break;
    }
    
    /* launch the nthread solvers */
    for (i = 0; i < nthreads; i++, solutionId++)
      threads[i] = new boost::thread(&balBifurcationDiagram::IntegrateAndEnqueue,this,lsol[i],solutionId);
    
    for (i = 0; i < nthreads; i++) {
      threads[i]->join();
      printf("%c%s", ESC, GREEN);
      printf("[%05d/%05d]\r", (prologue + cnt*nthreads + i + 1), total);
      printf("%c%s", ESC, NORMAL); fflush(stdout);
    }
  }
  
  printf("\n"); fflush(stdout);
  
  /**
   * we destroy the solvers and the parameters that were allocated previously
   **/
  for(i = 0; i < nthreads; i++) {
    lsol[i]->Destroy();
    lpar[i]->Destroy();
  }
  
  /* interrupting logger_thread */
  logger_thread->interrupt();
  
  /* waiting logger thread to writing the remaining solution in the queue and exit */
  logger_thread->join();
}

void balBifurcationDiagram::IntegrateAndEnqueue(balODESolver * sol, int solutionId) {
  sol->Solve();
  balSolution *solution = sol->GetSolution();
  solution->SetID(solutionId);
  classification->push_back(new balClassificationEntry(solution));
  
  {
    boost::mutex::scoped_lock lock(list_mutex);
    while (solution_list.size() >= LIST_MAX_SIZE) {  
      
      /**
       * we notify the thread that will write (which is waiting on q_full)
       * that the solution list in now full
       **/
      q_full.notify_one();

      /**
       * this thread goes on wait on q_empty conditional variable,
       * thus unlocking list_mutex until the thread that writes has emptied the list
       **/
      q_empty.wait(lock);
    }
    // insert a new solution into the list
    solution_list.push_back(solution);
  }	
  
  if(! restart_from_x0) {
    sol->SetX0(sol->GetXEnd());
  }
}
