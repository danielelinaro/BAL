/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balBifurcationDiagram.cpp
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
 * \file BifurcationDiagram.cpp
 * \brief Implementation of the class BifurcationDiagram
 */

#include "balBifurcationDiagram.h"

namespace bal {

void ResetColours(int d) {
  printf("%c%s", ESC, RED);
  printf("Aborting...\n");
  printf("%c%s", ESC, NORMAL);
  // !!!!! NEED TO ADD MEMORY MANAGEMENT (I.E. FREEING USED MEMORY BEFORE QUITTING) !!!!!
  exit(1);
}

bool CompareSummaryEntries(SummaryEntry *entry1, SummaryEntry *entry2) {
  return entry1->GetID() < entry2->GetID();
}

/***** SummaryEntry *****/

SummaryEntry::SummaryEntry(Solution *sol, diagram_mode mode) {
  int r, c, np, nx;
  realtype *buffer;
  np = sol->GetParameters()->GetNumber();
  n = np; // number of parameters
  if(mode == IC) {
    sol->GetSize(&r,&c);
    nx = c-2;
    n += 2*nx; // initial and final condition
  }
  n++; // classification field
  data = new double[n];
  id = sol->GetID();
  for(int i=0; i<np; i++)
    data[i] = sol->GetParameters()->At(i);
  switch(mode) {
  case MLE:
    data[n-1] = sol->GetLyapunovExponents()[0]; //saving Maximal Lyapunov exponent (MLE)
    break;
  case IC:
    buffer = sol->GetData();
    for(int i=0; i<nx; i++) {
      data[np+i] = buffer[1+i];
      data[np+nx+i] = buffer[(r-1)*c+1+i];
    }
    data[n-1] = sol->GetNumberOfTurns();
    break;
  default:
    data[n-1] = sol->GetNumberOfTurns();
  }
}

SummaryEntry::~SummaryEntry() {
  delete data;
}

int SummaryEntry::GetN() const {
  return n;
}

int SummaryEntry::GetID() const {
  return id;
}

double* SummaryEntry::GetData() const {
  return data;
}

/***** BifurcationDiagram *****/

BifurcationDiagram::BifurcationDiagram() : logger(new H5Logger()), solver(new ODESolver()) {
  restart_from_x0 = true;
  nthreads = 2;
  mode = PARAMS;
  nX0 = 0;
  X0 = NULL;
  signal(SIGINT, ResetColours);
}

BifurcationDiagram::~BifurcationDiagram() {
  delete logger;
  delete solver;
}

/*
Object* BifurcationDiagram::Clone() const {
  return new BifurcationDiagram(*this);
}
*/

void BifurcationDiagram::SetNumberOfThreads(int _nthreads) {
  if(_nthreads > 0)
    nthreads = _nthreads;
}

int BifurcationDiagram::GetNumberOfThreads() const {
  return nthreads;
}

void BifurcationDiagram::SetDynamicalSystem(DynamicalSystem *sys) {
  system = sys;
  ndim = system->GetDimension();
  solver->SetDynamicalSystem(system);
}

void BifurcationDiagram::SetDynamicalSystem(DynamicalSystem& sys) {
  SetDynamicalSystem(&sys);
}

DynamicalSystem* BifurcationDiagram::GetDynamicalSystem() const {
  return system;
}

void BifurcationDiagram::SetBifurcationParameters(BifurcationParameters *pars) {
  parameters = pars;
  npar = parameters->GetNumber();
  system->SetParameters(parameters);
}

void BifurcationDiagram::SetBifurcationParameters(BifurcationParameters& pars) {
  SetBifurcationParameters(&pars);
}

BifurcationParameters* BifurcationDiagram::GetBifurcationParameters() const {
  return parameters;
}

void BifurcationDiagram::SetLogger(Logger *log) {
  logger = log;
}

Logger* BifurcationDiagram::GetLogger() const {
  return logger;
}

void BifurcationDiagram::SetODESolver(ODESolver *sol) {
  solver = sol;
}

ODESolver* BifurcationDiagram::GetODESolver() const {
  return solver;
}

void BifurcationDiagram::SetFilename(const char * filename, bool compress) {
  logger->Open(filename,compress);
}

std::string BifurcationDiagram::GetFilename() {
  return logger->GetFilename();
}

bool BifurcationDiagram::SaveSummaryData(const char *filename) {
  if (!summary.size())
    return false;
  
  // CHECK
  //summary.sort(CompareSummaryEntries);
  
  int i;
  double *entry;
  std::list<SummaryEntry *>::iterator it;
  FILE *fid;
  
  fid = fopen(filename,"w"); 
  if(fid == NULL)
    return false;
  
  for(it=summary.begin(); it!=summary.end(); it++) {
    entry = (*it)->GetData();
    for(i=0; i<(*it)->GetN()-1; i++)
      fprintf(fid, "%e ", entry[i]);
    if (mode == MLE)
      fprintf(fid,"%e\n", (double) entry[(*it)->GetN()-1]);
    else 
      fprintf(fid,"%d\n", (int) entry[(*it)->GetN()-1]);
  }
  
  fclose(fid);
  return true;
}

double** BifurcationDiagram::GetSummaryData(int *size) {
  if(!summary.size())
    return NULL;
  
  summary.sort(CompareSummaryEntries);
  
  int i, j;
  double **data, *entry;
  std::list<SummaryEntry *>::iterator it;

  data = new double*[summary.size()];
  it = summary.begin();
  if(size != NULL) {
    size[0] = summary.size();
    size[1] = (*it)->GetN();
  }
  for(i=0; it!=summary.end(); it++, i++) {
    entry = (*it)->GetData();
    data[i] = new double[(*it)->GetN()];
    for(j=0; j<(*it)->GetN(); j++)
      data[i][j] = entry[j];
  }
  return data;
}

void BifurcationDiagram::SetInitialConditions(int nx0, double **x0) {
  if(nx0 > 0) {
    nX0 = nx0;
    X0 = x0;
  }
}

bool BifurcationDiagram::RestartsFromX0() const {
  return restart_from_x0;
}

void BifurcationDiagram::RestartFromX0(bool restart) {
  restart_from_x0 = restart;
}

bool BifurcationDiagram::SetMode(diagram_mode _mode) {
  if (mode != PARAMS && mode != IC && mode != MLE)
    return false;
  mode = _mode;
  return true;
}

int BifurcationDiagram::GetMode() const {
  return mode;
}

void BifurcationDiagram::ComputeDiagram() {
  solutions.clear();
  summary.clear();
  ComputeDiagramMultiThread();
}

void BifurcationDiagram::ComputeDiagramMultiThread() {
  int i, idx, cnt, solutionId;
  // the total number of integrations
  int total;
  // the number of integrations that must be performed outside of the main loop
  int prologue;
  // the number of for loops that will be performed after the prologue
  int nloops;

  if (solver->GetIntegrationMode() == LYAP)
    restart_from_x0 = true;
  
  switch(mode) {
    case PARAMS:
      parameters->Reset();
      total = parameters->GetTotalNumberOfTuples();
      break;
    case IC:
      RestartFromX0(true);
      idx = 0;
      solver->SetX0(X0[idx]);
      total = nX0;
      break;
    default:
      fprintf(stderr, "Unknown mode of operation in BifurcationDiagram::ComputeDiagramMultiThread.\n");
      return;
  }
  
  nloops = total / nthreads;
  prologue = total % nthreads;

  /**
   * launching the thread of the writing routine: it activates only when list size >= LIST_MAX_SIZE
   **/
  if (solver->GetIntegrationMode() != LYAP) {
    logger_thread = new boost::thread(&LoggerThread,logger,&solutions,&list_mutex,&q_empty,&q_full);
  }
  /**
   * calculating the first total % nthreads solutions in a serial fashion
   **/
  for(cnt = 0, solutionId = 1; cnt < prologue; cnt++, solutionId++) {
    IntegrateAndEnqueue(solver,solutionId);
    printf("%c%s", ESC, GREEN);
#ifdef DEBUG
    if(mode == IC)
      printf("[%05d/%05d]\r", idx+1, total);
    else
      printf("[%05d/%05d]\r", cnt+1, total);
#else
    printf("[%05d/%05d]\r", cnt+1, total);
#endif
    printf("%c%s", ESC, NORMAL); fflush(stdout);
    switch(mode) {
    case PARAMS:
      parameters->Next();
      break;
    case IC:
      idx++;
      if (idx<total)
        solver->SetX0(X0[idx]);
      break;
    }
  }

  /*
   * we make nthreads copies of both the solvers and the parameters of the dynamical system
   */
  ODESolver *lsol[nthreads];
  BifurcationParameters *lpar[nthreads];
  for(i = 0; i < nthreads; i++) {
    lpar[i] = parameters->Clone();
    lsol[i] = solver->Clone();
    // CHECK: not sure the following two lines are correct
    lsol[i]->SetDynamicalSystem(solver->GetDynamicalSystem()->Clone());
    lsol[i]->GetDynamicalSystem()->SetParameters(lpar[i]);
  }

  // the array of pointers to the threads
  boost::thread *threads[nthreads];

  for (cnt = 0; cnt < nloops; cnt++) {
    /* we set to each solver a different tuple of parameters */
    switch(mode) {
    case PARAMS:
      for(i = 0; i < nthreads; i++, parameters->Next())
	lpar[i]->CopyValues(parameters);
      break;
    case IC:
      for(i = 0; i < nthreads; i++, idx++)
	lsol[i]->SetX0(X0[idx]);
      break;
    }
    /* launch the nthread solvers */
    for (i = 0; i < nthreads; i++, solutionId++)
      threads[i] = new boost::thread(&BifurcationDiagram::IntegrateAndEnqueue,this,lsol[i],solutionId);
		
    
    for (i = 0; i < nthreads; i++) {
      threads[i]->join();
      delete threads[i];
      printf("%c%s", ESC, GREEN);
#ifdef DEBUG
      if(mode == IC)
	printf("[%05d/%05d]\r", idx-nthreads+i+1, total);
      else
	printf("[%05d/%05d]\r", (prologue + cnt*nthreads + i + 1), total);
#else
      printf("[%05d/%05d]\r", (prologue + cnt*nthreads + i + 1), total);
#endif
      printf("%c%s", ESC, NORMAL); fflush(stdout);
    }
  }
  
  printf("\n"); fflush(stdout);
  
  /**
   * we destroy the solvers and the parameters that were allocated previously
   **/
  for (i=0; i<nthreads; i++) {
    DynamicalSystem *ds = lsol[i]->GetDynamicalSystem();
    delete lsol[i];
    delete ds;
    delete lpar[i];
  }

  if (solver->GetIntegrationMode() != LYAP) {
    /* interrupting logger_thread */
    logger_thread->interrupt();
    /* waiting logger thread to writing the remaining solution in the queue and exit */
    logger_thread->join();
    delete logger_thread;
  }
}

void BifurcationDiagram::IntegrateAndEnqueue(ODESolver * sol, int solutionId) {
  sol->Solve();
  Solution *solution = sol->GetSolution();
  solution->SetID(solutionId);
  
  {
    boost::mutex::scoped_lock lock(list_mutex);
    while (solutions.size() >= LIST_MAX_SIZE) {  
      
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
    if (solver->GetIntegrationMode() != LYAP)
      solutions.push_back(solution);
    // insert a new summary into the list
    summary.push_back(new SummaryEntry(solution,mode));
  }	
  
  if (solver->GetIntegrationMode() == LYAP)
    delete solution;
  
  if(! restart_from_x0)
    sol->SetX0(sol->GetXEnd());

}

} // namespace bal

