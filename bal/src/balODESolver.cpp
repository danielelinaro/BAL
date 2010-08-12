/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balODESolver.cpp
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
 * \file balODESolver.cpp
 * \brief Implementation of the class balODESolver
 */

#include "balODESolver.h"

const char * balODESolver::GetClassName() const {
  return "balODESolver";
}

balODESolver * balODESolver::Create() {
  return new balODESolver;
}

balODESolver * balODESolver::Copy (balODESolver * solver) {
  return new balODESolver(*solver);
}

void balODESolver::Destroy() {
  delete this;
}

balODESolver::balODESolver() {
  neq = npar = nev = 0;
  buffer = NULL;
  delete_buffer = false;
  events = NULL;
  events_constraints = NULL;
  dynsys = NULL;
  params = NULL;
  reltol = RTOL;
  abstol = ATOL;
  mode = balTRAJ;
  stiff = true;
  t = 0.0;
  tstep = STEP;
  ttran = T_TRAN;
  tfinal = T_END;
  lyap_tstep = 100*STEP;
  setup = false;
  max_intersections = DEFAULT_INTERSECTIONS;
  bufsize = 0;
  rows = 0;
  cols = 0;
  cvode_mem = NULL;
  errfp = fopen (ERROR_FILE, "w");
  halt_at_equilibrium = false;
  halt_at_cycle = false;
  nturns = -1;
  equilibrium_tolerance = EQUIL_TOL;
  cycle_tolerance = CYCLE_TOL;
  x = NULL;
  xdot = NULL;
  x0 = NULL;
  x_inters = NULL;
  lyapunov_exponents = NULL;
	delete_lyapunov_exponents = false;
  nvectors_allocated = false;
  class_event = 1;
}

balODESolver::balODESolver(const balODESolver& solver) {
  int i;
	buffer = NULL;
  delete_buffer = false;
  rows = 0;
  nvectors_allocated = false;
  SetDynamicalSystem((solver.GetDynamicalSystem())->Copy());
  for(i=0; i<dynsys->GetDimension(); i++)
    Ith(x0,i) = Ith(solver.x0,i);
  
  reltol = solver.GetRelativeTolerance();
  abstol = solver.GetAbsoluteTolerance();
  mode = solver.GetIntegrationMode();
  stiff = solver.stiff;
  t = 0.0;
  tstep = solver.GetTimeStep();
  ttran = solver.GetTransientDuration();
  tfinal = solver.GetFinalTime();
  lyap_tstep = solver.GetLyapunovTimeStep();
	delete_lyapunov_exponents = false;
	if (solver.delete_lyapunov_exponents) {
		lyapunov_exponents = new realtype[dynsys->GetOriginalDimension()];
		delete_lyapunov_exponents = true;
		for(i=0; i<dynsys->GetOriginalDimension(); i++)
			lyapunov_exponents[i] = solver.GetLyapunovExponents()[i];
	}
  setup = false;
  max_intersections = solver.GetMaxNumberOfIntersections();
  bufsize = 0;
  cvode_mem = NULL;
  errfp = fopen (ERROR_FILE, "a");
  halt_at_equilibrium = solver.HaltsAtEquilibrium();
  halt_at_cycle = solver.HaltsAtCycle();
  nturns = solver.GetNumberOfTurns();
  equilibrium_tolerance = solver.GetEquilibriumTolerance();
  cycle_tolerance = solver.GetCycleTolerance();
  class_event = solver.GetClassificationEvent();
}

balODESolver::~balODESolver() {
  if(delete_buffer || buffer != NULL) {
    delete buffer;
  }
  if(events != NULL) {
    delete events;
  }
  if(events_constraints != NULL) {
    delete events_constraints;
  }
  if(nvectors_allocated) {
    N_VDestroy_Serial(x);
    N_VDestroy_Serial(xdot);
    N_VDestroy_Serial(x0);
    N_VDestroy_Serial(x_inters);
  }
  if(cvode_mem != NULL) {
    CVodeFree(&cvode_mem);
  }
  fclose(errfp);

  if(delete_lyapunov_exponents){
    delete lyapunov_exponents;
  }
}

double * balODESolver::GetBuffer() const {
  return buffer;
}

int balODESolver::GetBufferSize() const {
  return bufsize;
}

void balODESolver::GetBufferSize(int * r, int * c) const {
  *r = rows;
  *c = cols;
}

void balODESolver::SetDynamicalSystem(balDynamicalSystem * ds) {
  if(ds != NULL) {
    dynsys = ds;
    neq	= dynsys->GetDimension();
    if(dynsys->HasEvents()) {
      nev = dynsys->GetNumberOfEvents();
      events = new int[nev];
      if(dynsys->HasEventsConstraints()) {
		// if constraints are presents, their number must be equal to the
		// number of events. if, for some reason, the i-th constraint makes
		// no sense, then it is sufficient that the dynamical system class
		// always return events_constraints[i] = 1.
				events_constraints = new int[nev];
      }
    }
    params = dynsys->GetParameters();
    npar = params->GetNumber();
    // the number of columns of the buffer is equal to the number of dimensions of the system
    // plus 2, i.e. the time instant and a label that describes the type of
    // record
    cols = neq + 2;
    if(nvectors_allocated) {
      N_VDestroy_Serial(x);
      N_VDestroy_Serial(xdot);
      N_VDestroy_Serial(x0);
      N_VDestroy_Serial(x_inters);
      nvectors_allocated = false;
    }
    x = N_VNew_Serial(neq);
    xdot = N_VNew_Serial(neq);
    x0 = N_VNew_Serial(neq);
    x_inters = N_VNew_Serial(neq);
    nvectors_allocated = true;
    // perform setup!
    setup = false;
  }
}

balDynamicalSystem * balODESolver::GetDynamicalSystem() const {
  return dynsys;
}

void balODESolver::SetDynamicalSystemParameters(balParameters * par){
  dynsys->SetParameters(par);
  params = par;
}

balSolution * balODESolver::GetSolution() const {
  if(rows == 0)
    return NULL;
  balSolution * solution = balSolution::Create();
  solution->SetData(rows,cols,buffer);
  solution->SetParameters(params);
	if (mode == balLYAP)
		solution->SetLyapunovExponents(dynsys->GetOriginalDimension(),lyapunov_exponents);
  else 
		solution->SetNumberOfTurns(nturns);
  return solution;
}

realtype balODESolver::GetTransientDuration () const {
  return ttran;
}

void balODESolver::SetTransientDuration (realtype tran) {
  if (tran >= 0) {
    ttran = tran;
    if(buffer != NULL) {
      delete_buffer = true;
    }
  }
}

realtype balODESolver::GetFinalTime () const {
  return tfinal;
}

void balODESolver::SetFinalTime (realtype final) {
  if (final >= 0) {
    tfinal = final;
    if(buffer != NULL) {
      delete_buffer = true;
    }
  }
}

realtype balODESolver::GetTimeStep () const {
  return tstep;
}

void balODESolver::SetTimeStep (realtype step) {
  if (step > 0) {
    tstep = step;
    if(buffer != NULL) {
      delete_buffer = true;
    }
  }
}

realtype balODESolver::GetLyapunovTimeStep () const {
  return lyap_tstep;
}

void balODESolver::SetLyapunovTimeStep (realtype tstep) {
  if (tstep > 0)
    lyap_tstep = tstep;
}

realtype balODESolver::GetRelativeTolerance () const {
  return reltol;
}

void balODESolver::SetRelativeTolerance (realtype rtol) {
  if (rtol > 0)
    reltol = rtol;
}

realtype balODESolver::GetAbsoluteTolerance () const {
  return abstol;
}

void balODESolver::SetAbsoluteTolerance (realtype atol) {
  if (atol > 0)
    abstol = atol;
}

integration_mode balODESolver::GetIntegrationMode () const {
  return mode;
}

void balODESolver::SetIntegrationMode (integration_mode m) {
  mode = m;
  if(buffer != NULL) {
    delete_buffer = true;
  }
}

void balODESolver::IsStiff(bool stiffness) {
  stiff = stiffness;
}

int balODESolver::GetMaxNumberOfIntersections() const {
  return max_intersections;
}

void balODESolver::SetMaxNumberOfIntersections(int intersections) {
  if(intersections > 0) {
    max_intersections = intersections;
    if(buffer != NULL) {
      delete_buffer = true;
    }
  }
}

void balODESolver::HaltAtEquilibrium(bool halt) {
  halt_at_equilibrium = halt;
}

bool balODESolver::HaltsAtEquilibrium() const {
  return halt_at_equilibrium;
}

void balODESolver::HaltAtCycle(bool halt) {
  halt_at_cycle = halt;
}

bool balODESolver::HaltsAtCycle() const {
  return halt_at_cycle;
}

void balODESolver::SetClassificationEvent(int ce) {
  if(ce > 0)
    class_event = ce;
}

int balODESolver::GetClassificationEvent() const {
  return class_event;
}

int balODESolver::GetNumberOfTurns() const {
  return nturns;
}

realtype balODESolver::GetEquilibriumTolerance() const {
  return equilibrium_tolerance;
}

void balODESolver::SetEquilibriumTolerance(realtype tol) {
  if(tol > 0)
    equilibrium_tolerance = tol;
}

realtype balODESolver::GetCycleTolerance() const {
  return cycle_tolerance;
}

void balODESolver::SetCycleTolerance(realtype tol) {
  if(tol > 0)
    cycle_tolerance = tol;
}

N_Vector balODESolver::GetX() const {
  return x;
}

N_Vector balODESolver::GetXdot() const {
  return xdot;
}

N_Vector balODESolver::GetX0() const {
  return x0;
}

realtype * balODESolver::GetXEnd() const {
  if(rows == 0)
    return NULL;
  return buffer + (rows-1)*cols + 1;
}

void balODESolver::SetX0(N_Vector X0, int n) {
  if(X0 != NULL && dynsys != NULL) {
    int stop;
    if(n == -1)
      stop = dynsys->GetDimension();
    else
      stop = n;
    for(int i=0; i<stop; i++)
      Ith(x0,i) = Ith(X0,i);
  }
}

void balODESolver::SetX0(realtype * X0, int n) {
  if(X0 != NULL && dynsys != NULL) {
    int stop;
    if(n == -1)
      stop = dynsys->GetDimension();
    else
      stop = n;
    for(int i=0; i<stop; i++)
      Ith(x0,i) = X0[i];
  }
}

inline void balODESolver::SetOrthonormalBaseIC() {
  int length = dynsys->GetOriginalDimension();
  for(int i=0; i<length; i++) {
    for(int j=0; j<length; j++)
      Ith(x0,length+i*length+j) = (i==j ? 1.0 : 0.0) ;
  }
}

realtype * balODESolver::GetLyapunovExponents() const {
  return lyapunov_exponents;
}

bool balODESolver::Setup() {
  int flag;
  
  if (cvode_mem != NULL) {
    CVodeFree (&cvode_mem);
  }
  
  if(stiff)
    cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
  else
    cvode_mem = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);

  if (cvode_mem == NULL) {
    fprintf (stderr, "error on CVodeCreate.\n");
    return false;
  }
  
#ifdef CVODE25
  flag = CVodeMalloc (cvode_mem, balDynamicalSystem::RHSWrapper, 0.0, x, CV_SS, reltol, &abstol);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeMalloc.\n");
    return false;
  }
#endif
#ifdef CVODE26
  flag = CVodeInit (cvode_mem, balDynamicalSystem::RHSWrapper, 0.0, x);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeInit.\n");
    return false;
  }
  flag = CVodeSStolerances (cvode_mem, reltol, abstol);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeSStolerances.\n");
    return false;
  }
#endif
  
  flag = CVodeSetErrFile (cvode_mem, errfp);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeSetErrFile.\n");
    return false;
  }

#ifdef CVODE25
  flag = CVodeSetFdata (cvode_mem, dynsys);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeSetFdata.\n");
    return false;
  }
#endif
#ifdef CVODE26
  flag = CVodeSetUserData (cvode_mem, dynsys);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeSetUserData.\n");
    return false;
  }
#endif

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense (cvode_mem, neq);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVDense.\n");
    return false;
  }

  if (dynsys->HasJacobian()) {
    /* Set the Jacobian routine to Jac (user-supplied) */
#ifdef CVODE25
    flag = CVDenseSetJacFn (cvode_mem, balDynamicalSystem::JacobianWrapper, dynsys);
    if (flag != CV_SUCCESS) {
      fprintf (stderr, "Error on CVDenseSetJacFn.\n");
      return false;
    }
#endif
#ifdef CVODE26
    flag = CVDlsSetDenseJacFn (cvode_mem, balDynamicalSystem::JacobianWrapper);
    if (flag != CV_SUCCESS) {
      fprintf (stderr, "Error on CVDlsSetDenseJacFn.\n");
      return false;
    }
#endif
  }

  CVodeSetMaxNumSteps (cvode_mem, MAX_NUM_STEPS);
  CVodeSetMaxErrTestFails (cvode_mem, 10);
  
  // Allocate memory for the integration buffer
  AllocateSolutionBuffer();
  
  setup = true;
  return true;
}

bool balODESolver::AllocateSolutionBuffer() {
  if(delete_buffer || buffer == NULL) {
    int lrows;
    if(delete_buffer) {
      delete buffer;
      delete_buffer = false;
    }
    switch(mode) {
    case balTRAJ:
		case balLYAP:
      lrows = (int) ceil ((tfinal - ttran) / tstep);
      if(lrows < 0) lrows = 0;
      lrows += 2;
      break;
    case balEVENTS:
      lrows = max_intersections + 2;
      break;
    case balBOTH:
      lrows = (int) ceil ((tfinal - ttran) / tstep);
      // the number of detected intersections is increased only after
      // transient evolution.
      if(lrows < 0) 
				lrows = 2;
      else
				lrows += max_intersections + 2;
    }
    bufsize = lrows * cols;
    try {
      buffer = new realtype[bufsize];
    } catch (bad_alloc&) {
      fprintf (stderr, "Not enough memory to allocate for the solution buffer...\n");
      return false;
    }
  }
  return true;
}

void balODESolver::ResetInitialCondition() {
  t = 0.0;
  for(int i=0; i<neq; i++) {
    Ith(x,i) = Ith(x0,i);
  }
  ResetPositionInBuffer();
  StoreRecordInBuffer(balSTART);
}

void balODESolver::ResetPositionInBuffer() {
  rows = 0;
}

void balODESolver::StoreRecordInBuffer(int lbl) {
  int actual_pos = rows*cols;
  buffer[actual_pos] = t;
  for(int i=0; i<neq; i++)
    buffer[actual_pos+i+1] = Ith(x,i);
  buffer[actual_pos+neq+1] = (realtype) lbl;
  rows++;
}

void balODESolver::ChangeCurrentLabel(int lbl) {
  buffer[rows*cols-1] = (realtype) lbl;
}

void balODESolver::SkipTransient(bool *equilibrium, bool *error) {
  int flag;
  *equilibrium = false;
  *error = false;

  while (t < ttran) {
    flag = CVode (cvode_mem, ttran, x, &t, CV_NORMAL);
    if (flag < 0 && flag != CV_TOO_MUCH_WORK && (flag != CV_ILL_INPUT || mode == balTRAJ)) {
      *error = true;
      break;
    }
    flag = CheckEquilibrium();
    if(flag != EQUIL_FALSE)
      *equilibrium = true;
    if(flag == EQUIL_BREAK)
      break;
  }
  // if it's an equilibrium point, the label is changed at the end...
  StoreRecordInBuffer(balTRAN_END);
}

int balODESolver::CheckEquilibrium() {
  dynsys->RHS (t, x, xdot, params);
  if(EuclideanDistance(neq,xdot) < equilibrium_tolerance) {
    nturns = 0;
    if(halt_at_equilibrium) {
      ChangeCurrentLabel(balEQUIL);
      return EQUIL_BREAK;
    }
    return EQUIL_TRUE;
  }
  return EQUIL_FALSE;
}

int balODESolver::CheckCycle(int guess) {
  if(EuclideanDistance(neq,x,x_inters) < cycle_tolerance) {
    if(guess > 0) {
      nturns = guess;
      if(halt_at_cycle)
	return CYCLE_BREAK;
    }
    return CYCLE_TRUE;
  }
  return CYCLE_FALSE;
}

void balODESolver::SetSolutionLength(int length) {
  if (length > 0) {
    rows = length / cols;
  }
}

realtype balODESolver::EuclideanDistance(int length, N_Vector x, N_Vector y) const {
  realtype dst = 0.0;
  if(y != NULL) {
    for(int i=0; i<length; i++)
      dst += (Ith(x,i) - Ith(y,i)) * (Ith(x,i) - Ith(y,i));
  }
  else {
    for(int i=0; i<length; i++)
      dst += Ith(x,i) * Ith(x,i);
  }
  return sqrt(dst);
}

bool balODESolver::GramSchmidtOrthonorm(realtype * x, realtype * xnorm, realtype * znorm) const {
  int n = dynsys->GetOriginalDimension();
  int i,j,k;
	realtype dp;
  realtype * tmp = new realtype[n];
 
  znorm[0] = Norm(n,x);
  for (k=0; k<n; k++)
    xnorm[k] = x[k]/znorm[0];
  
  for (i=1; i<n ; i++) {
    
    for (k=0; k<n; k++)
      tmp[k] = x[n*i+k]; 
    
    for (j=0; j<i; j++) {
			dp = DotProduct(n,x+(n*i),xnorm+(n*j));
      for (k=0; k<n; k++)
				tmp[k] = tmp[k] - dp * xnorm[n*j+k];
    }
    
    znorm[i] = Norm(n,tmp);
    
    for (k=0; k<n; k++)
      xnorm[n*i+k] = tmp[k]/znorm[i];
  }

  //fprintf(stderr, "norms: %f %f %f\n", Norm(n,xnorm), Norm(n,xnorm+n), Norm(n,xnorm+2*n));
  //fprintf(stderr, "dots: %f %f %f\n", DotProduct(n,xnorm,xnorm+n), DotProduct(n,xnorm,xnorm+2*n), DotProduct(n,xnorm+n,xnorm+2*n));

  /*
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++)
      fprintf(stderr, "%f ", xnorm[i*n+j]);
    fprintf(stderr, "\n");
  }
  */

  delete tmp;
  return true;
} 

inline realtype balODESolver::Norm(int length, realtype * x) const {
  realtype norm = 0.0;
  for(int i = 0; i<length; i++)
    norm += x[i]*x[i];
  return sqrt(norm);
}  

inline realtype  balODESolver::DotProduct(int length, realtype * x, realtype * y) const {
  realtype res = 0.0;
  for(int i=0; i<length; i++)
    res += x[i]*y[i];
  return res;
}
 

bool balODESolver::ResetCVode() {
  int flag;

  switch(mode) {
  case balTRAJ:
  case balLYAP:
#ifdef CVODE25
    flag = CVodeRootInit (cvode_mem, 0, NULL, dynsys);
#endif
#ifdef CVODE26
    flag = CVodeRootInit (cvode_mem, 0, NULL);
#endif
    break;
  case balEVENTS:
  case balBOTH:
#ifdef CVODE25
    flag = CVodeRootInit (cvode_mem, nev, balDynamicalSystem::EventsWrapper, dynsys);
#endif
#ifdef CVODE26
    flag = CVodeRootInit (cvode_mem, nev, balDynamicalSystem::EventsWrapper);
#endif
  }
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeRootInit: flag = %d.\n", flag);
    return false;
  }

#ifdef CVODE25
  flag = CVodeReInit (cvode_mem, balDynamicalSystem::RHSWrapper, 0.0, x, CV_SS, reltol, &abstol);
#endif
#ifdef CVODE26:
  flag = CVodeReInit (cvode_mem, 0.0, x);
#endif
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return false;
  }
  return true;
}

bool balODESolver::Solve() {
  if(! setup)
    Setup();
  dynsys->Reset();
  switch (mode) {
  case balTRAJ:
    return SolveWithoutEvents();
  case balEVENTS:
  case balBOTH:
    return SolveWithEvents();
  case balLYAP:
    return SolveLyapunov();
  }
  return false;
}

bool balODESolver::SolveLyapunov() {
	
	int i;
	realtype * temp_x0 = new realtype[dynsys->GetOriginalDimension()];
	for(i=0; i<dynsys->GetOriginalDimension(); i++)
		temp_x0[i] = Ith(x0,i);
	realtype temp_ttran = ttran;
	realtype tend = tfinal;
  realtype t_ = 0.0;
	tfinal = ttran;
	delete_buffer = true;
	//fprintf(stderr,"ci: %lf\t%lf\t%lf\t%lf\n", Ith(x0,0),Ith(x0,1),Ith(x0,2),dynsys->GetParameters()->At(0));
	SolveWithoutEvents();
	//fprintf(stderr,"ttran esaurito stato: %lf\t%lf\t%lf\n", Ith(x,0),Ith(x,1),Ith(x,2));
	
	realtype * new_x0 = new realtype [dynsys->GetDimension()];
	for (i=0; i<dynsys->GetDimension(); i++)
		new_x0[i] = Ith(x,i);
		
	dynsys->Extend(true);
	SetDynamicalSystem(dynsys);
  
	int N = dynsys->GetDimension();
  int n = dynsys->GetOriginalDimension();
  realtype * x_ = new realtype[N];
  realtype * xnorm = new realtype[n*n];
  realtype * znorm = new realtype[n];
  realtype * cum = new realtype[n];
	if (!delete_lyapunov_exponents){
		lyapunov_exponents = new realtype[n];
		delete_lyapunov_exponents = true;
	}
  for (i=0; i<n; i++)
    cum[i] = 0.0;
	
	SetX0(new_x0,n);
	delete new_x0;
  SetOrthonormalBaseIC();
  
	//fprintf(stderr,"ci_ext: ");
//	for (i=0;i<N; i++)
//		fprintf(stderr,"%lf ",Ith(x0,i));
//	fprintf(stderr,"\n");
  
  tfinal = lyap_tstep;
  ttran = 0.0;
	
	delete_buffer = true;
	Setup();
	dynsys->Reset();
  while(t_ < tend){
    SolveWithoutEvents();
    t_ += lyap_tstep;
    for(i=0; i<N; i++)
      x_[i] = Ith(x,i);
    GramSchmidtOrthonorm(x_+n,xnorm,znorm);
    for(i=0; i<n*n; i++)
      x_[i+n] = xnorm[i];
    for(i=0; i<n; i++)
      cum[i] += log(znorm[i]);///log(2.0);
    SetX0(x_);
		
		/*if((int)t_ % 100 == 0){
			for(i=0; i<n; i++){
				lyapunov_exponents[i] = cum[i]/t_;
				printf("%e ",lyapunov_exponents[i]);
				printf(" %lf ",x_[i]);
			}
			printf("\n");
      }*/
	}
  
  for(i=0; i<n; i++){
		lyapunov_exponents[i] = cum[i]/t_;
		//printf("%e ",lyapunov_exponents[i]);
	}
	//printf("\n");
	
	ttran = temp_ttran;
	tfinal = tend;
	dynsys->Extend(false);
	SetDynamicalSystem(dynsys);
	delete_buffer = true;
	Setup();
	SetX0(temp_x0);
  
	delete temp_x0;
  delete x_;
  delete xnorm;
  delete znorm;
  delete cum;
}

bool balODESolver::SolveWithoutEvents() {
  int flag;
  realtype tout;
  bool eq, err;
  
  // check whether the buffer is allocated or it needs to be
  // reallocated
  AllocateSolutionBuffer();
  // set initial conditions (this automatically saves the i.c.
  // in the first row of the solution buffer)
  ResetInitialCondition();
  // set the number of turns of the solution to ''unknown''
  nturns = -1;

  if(! ResetCVode())
    return false;
	
  /************* TRANSIENT *****************/
  SkipTransient(&eq,&err);
  
  /************* TRAJECTORY *****************/
  if(!(eq && halt_at_equilibrium) && !err) {
    tout = t+tstep;
    while (tout < tfinal+tstep) {
      flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
      StoreRecordInBuffer(balSTEP);
      
      if (flag < 0 && flag != CV_TOO_MUCH_WORK) {
	/*
	 * this is because the flag CV_ILL_INPUT is returned when two
	 * events are found at a very small time interval, which is the
	 * case if at least one of the Poincare' sections corresponds to
	 * an extremum of a state variables, which is a condition always
	 * verified when the trajectory is on an equilibrium
	 */
				if(flag == CV_ILL_INPUT) {
					if(CheckEquilibrium() == EQUIL_BREAK)
						break;
				}
	/* this is done if the error is not CV_ILL_INPUT */
				err = true;
				break;
			}

      if(CheckEquilibrium() == EQUIL_BREAK)
				break;

      tout += tstep;
    }
  }

  /* if the integrator stopped because of an error, we change the label of the last row in
   * the integration buffer and stop the integration procedure */
  if (err) {
    ChangeCurrentLabel(balERROR);
    return false;
  }
  /* if the point is an equilibrium, we change the label of the last row in
   * the integration buffer and stop the integration procedure */
  if (eq || nturns == 0)
    ChangeCurrentLabel(balEQUIL);

  return true;
}

bool balODESolver::SolveWithEvents() {
  int flag, eq_flag, i, j;
  int intersections; // the total number of intersections with Poincare' sections
  int class_inters; // the number of intersections with the Poincare' section used to detect cycles
  int nturns_guess;
  realtype tout;
  bool eq, err, cycle, restart;
  
  // check whether the buffer is allocated or it needs to be
  // reallocated
  AllocateSolutionBuffer();
  // set initial conditions
  ResetInitialCondition();
  // set the number of turns of the solution to ''unknown''
  nturns = -1;
  
  if(! ResetCVode())
    return false;

  /************* TRANSIENT *****************/
  SkipTransient(&eq,&err);

  /************* INTERSECTIONS ************/
  restart = false;
  if(!(eq && halt_at_equilibrium) && !err) {
    if(mode == balEVENTS)
      tout = tfinal;
    else if(mode == balBOTH)
      tout = t+tstep;
    intersections = 0;
    class_inters = 0;
    nturns_guess = -1;
    cycle = false;
    while (intersections < max_intersections && t < tfinal) {
      flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
      if (flag < 0 && flag != CV_TOO_MUCH_WORK && flag != CV_ILL_INPUT) {
	/*
	 * this is because the flag CV_ILL_INPUT is returned when two
	 * events are found at a very small time interval, which is the
	 * case if at least one of the Poincare' sections corresponds to an
	 * extremum of a state variables, which is a condition always
	 * verified when the trajectory is on an equilibrium
	 */
	StoreRecordInBuffer(balERROR);
	err = true;
	break;
      }

      eq_flag = CheckEquilibrium();
      if(eq_flag == EQUIL_BREAK)
	break;
      else if(eq_flag == EQUIL_TRUE)
	eq = true;

      if (flag == CV_SUCCESS && mode == balBOTH) {
	StoreRecordInBuffer(balSTEP);
	tout += tstep;
      }
      else if (flag == CV_ROOT_RETURN && !eq) {
	CVodeGetRootInfo (cvode_mem, events);
	if(dynsys->HasEventsConstraints()) {
	  dynsys->EventsConstraints(t,x,events_constraints,params);
	  /* 
	   * call ManageEvents so that the dynamical system can (optionally)
	   * change something in its internal structure. See for example
	   * balPLL, a switch system.
	   */
	  dynsys->ManageEvents(t,x,events,events_constraints);
	}
	else {
	  dynsys->ManageEvents(t,x,events);
	}
	
	for (i = 0; i < nev; i++) {
#ifdef CVODE25
	  if (events[i] == 1) {
#endif
#ifdef CVODE26
	  // this is because the return values of CVodeGetRootInfo have
	  // changed and are not only {0,1}, but {-1,0,1}, according to the
	  // decrease or increase of the test function
	  if (events[i] != 0) {
#endif
	    if (! dynsys->HasEventsConstraints() || events_constraints[i] == 1) {
	      intersections++;
	      StoreRecordInBuffer(i+1);
	      
	      /*** start cycle detection ***/
	      if((i+1) == class_event) {
		class_inters++;
		if(class_inters == 1) {
		  // save the first intersection
		  for(j=0; j<neq; j++)
		    Ith(x_inters,j) = Ith(x,j);
		}
		else {
		  if(nturns_guess == -1) {
		    if(CheckCycle(nturns_guess) == CYCLE_TRUE) {
		      // guess the number of turns
		      nturns_guess = class_inters-1;
		    }
		    else if(intersections >= round(max_intersections/2) && !restart) {
		      // restart to keep track of intersections if we reach
		      // half of the allowed intersections and we haven't
		      // made a guess yet
		      class_inters = 0;
		      restart = true;
		    }
		  }
		  else if(class_inters == 2*nturns_guess + 1) {
		    flag = CheckCycle(nturns_guess);
		    if(flag == CYCLE_FALSE) {
		      // restart to keep track of intersections and see if
		      // we can detect a new cycle
		      nturns_guess = -1;
		      class_inters = 0;
		    }
		    else if(flag == CYCLE_BREAK) {
		      cycle = true;
		      break;
		    }
		  }
		}
	      }
	      /*** end cycle detection ***/
	    }
	  }
	}
	if(cycle) break;
      }
    }
  }
  /* if the integrator stopped because of an error, we change the label of the last row in
   * the integration buffer and stop the integration procedure */
  if (err) {
    ChangeCurrentLabel(balERROR);
    return false;
  }
  
  /* if we weren't able to detect a periodicity, there's a good chance 
   * that the system be chaotic */
  if (nturns == -1)
    nturns = max_intersections;
  
  /* if the point is an equilibrium, we change the label of the last row in
   * the integration buffer */
  if (eq || nturns == 0)
    ChangeCurrentLabel(balEQUIL);
    
  return true;
}

