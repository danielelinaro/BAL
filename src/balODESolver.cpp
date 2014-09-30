/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balODESolver.cpp
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
 * \file balODESolver.cpp
 * \brief Implementation of the class ODESolver
 */

#include <fstream>
#include "balODESolver.h"

namespace bal {

ODESolver::ODESolver() {
#ifdef DEBUG
  std::cout << "ODESolver constructor.\n";
#endif
  neq = npar = nev = 0;
  reltol = RTOL;
  abstol = ATOL;
  mode = TRAJ;
  stiff = true;
  t0 = 0.0;
  t = t0;
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
  nvectors_allocated = false;
  class_event = 1;
}

ODESolver::ODESolver(const ODESolver& solver) {
#ifdef DEBUG
  std::cout << "ODESolver copy constructor.\n";
#endif
  rows = 0;
  nvectors_allocated = false;
  SetDynamicalSystem(dynamic_cast<DynamicalSystem*>(solver.dynsys->Clone()));
  for(int i=0; i<dynsys->GetDimension(); i++)
    Ith(x0,i) = Ith(solver.x0,i);
  reltol = solver.reltol;
  abstol = solver.abstol;
  mode = solver.mode;
  stiff = solver.stiff;
  t0 = solver.t0;
  t = t0;
  tstep = solver.tstep;
  ttran = solver.ttran;
  tfinal = solver.tfinal;
  lyap_tstep = solver.lyap_tstep;
  lyapunov_exponents = boost::shared_array<realtype>(new realtype[dynsys->GetOriginalDimension()]);
  for(int i=0; i<dynsys->GetOriginalDimension(); i++)
    lyapunov_exponents[i] = solver.lyapunov_exponents[i];
  max_intersections = solver.max_intersections;
  bufsize = 0;
  errfp = fopen(ERROR_FILE, "a");
  halt_at_equilibrium = solver.halt_at_equilibrium;
  halt_at_cycle = solver.halt_at_cycle;
  nturns = solver.nturns;
  equilibrium_tolerance = solver.equilibrium_tolerance;
  cycle_tolerance = solver.cycle_tolerance;
  class_event = solver.class_event;
  cvode_mem = NULL;

  // we force cvode_mem to be re-allocated.
  setup = false;
}

ODESolver::~ODESolver() {
#ifdef DEBUG
  std::cout << "ODESolver destructor.\n";
#endif
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
}

boost::shared_array<realtype> ODESolver::GetBuffer() const {
  return buffer;
}

int ODESolver::GetBufferSize() const {
  return bufsize;
}

void ODESolver::GetBufferSize(int *r, int *c) const {
  *r = rows;
  *c = cols;
}

void ODESolver::SetDynamicalSystem(DynamicalSystem *ds) {
  if(nvectors_allocated) {
    N_VDestroy_Serial(x);
    N_VDestroy_Serial(xdot);
    N_VDestroy_Serial(x0);
    if (nev)
      N_VDestroy_Serial(x_inters);
    nvectors_allocated = false;
  }

  dynsys = boost::shared_ptr<DynamicalSystem>(dynamic_cast<DynamicalSystem*>(ds->Clone()));
  neq = dynsys->GetDimension();
  nev = dynsys->GetNumberOfEvents();
  if(nev) {
    events = boost::shared_array<int>(new int[nev]);
    if(dynsys->HasEventsConstraints()) {
	// if constraints are presents, their number must be equal to the
	// number of events. if, for some reason, the i-th constraint makes
	// no sense, then it is sufficient that the dynamical system class
	// always return events_constraints[i] = 1.
      events_constraints = boost::shared_array<int>(new int[nev]);
    }
  }
  //~~
  //params = dynsys->GetParameters();
  //npar = params->GetNumber();
  //~~
  
  // the number of columns of the buffer is equal to the number of dimensions of the system
  // plus 2, i.e. the time instant and a label that describes the type of
  // record
  cols = neq + 2;
  printf("%d\n", neq);
  x = N_VNew_Serial(neq);
  xdot = N_VNew_Serial(neq);
  x0 = N_VNew_Serial(neq);
  x_inters = N_VNew_Serial(neq);
  nvectors_allocated = true;
  // perform setup!
  setup = false;
}

boost::shared_ptr<DynamicalSystem> ODESolver::GetDynamicalSystem() const {
  return dynsys;
}

void ODESolver::SetDynamicalSystemParameters(boost::shared_ptr<Parameters>& par){
  dynsys->SetParameters(par);
  params = par;
}

Solution* ODESolver::GetSolution() const {
  if(rows == 0)
    return NULL;
  Solution *solution = new Solution(rows,cols,buffer.get());
  solution->SetParameters(dynsys->GetParameters().get());
  if (mode == LYAP)
    solution->SetLyapunovExponents(lyapunov_exponents.get());
  else 
    solution->SetNumberOfTurns(nturns);
  return solution;
}

realtype ODESolver::GetInitialTime() const {
  return t0;
}

void ODESolver::SetInitialTime(realtype T0) {
  if(T0 >= 0)
    t0 = T0;
}

realtype ODESolver::GetTransientDuration () const {
  return ttran;
}

void ODESolver::SetTransientDuration (realtype tran) {
  if (tran >= 0)
    ttran = tran;
}

realtype ODESolver::GetFinalTime () const {
  return tfinal;
}

void ODESolver::SetFinalTime (realtype final) {
  if (final >= 0)
    tfinal = final;
}

realtype ODESolver::GetTimeStep () const {
  return tstep;
}

void ODESolver::SetTimeStep (realtype step) {
  if (step > 0)
    tstep = step;
}

realtype ODESolver::GetLyapunovTimeStep () const {
  return lyap_tstep;
}

void ODESolver::SetLyapunovTimeStep (realtype tstep) {
  if (tstep > 0)
    lyap_tstep = tstep;
}

realtype ODESolver::GetRelativeTolerance () const {
  return reltol;
}

void ODESolver::SetRelativeTolerance (realtype rtol) {
  if (rtol > 0)
    reltol = rtol;
}

realtype ODESolver::GetAbsoluteTolerance () const {
  return abstol;
}

void ODESolver::SetAbsoluteTolerance (realtype atol) {
  if (atol > 0)
    abstol = atol;
}

integration_mode ODESolver::GetIntegrationMode () const {
  return mode;
}

void ODESolver::SetIntegrationMode (integration_mode m) {
  mode = m;
}

void ODESolver::IsStiff(bool stiffness) {
  stiff = stiffness;
}

int ODESolver::GetMaxNumberOfIntersections() const {
  return max_intersections;
}

void ODESolver::SetMaxNumberOfIntersections(int intersections) {
  if(intersections > 0)
    max_intersections = intersections;
}

bool ODESolver::SaveOrbit(const char *filename) const {
  if(mode != BOTH)
    return false;

  std::ofstream ofs(filename);
  if(ofs.fail())
    return false;

  if(nturns == max_intersections) {
    DUMPBUFFER(ofs);
  }
  else {
    int i, cnt, start, stop;
    cnt = nturns;
    i = rows*cols-1;
    while(cnt >= 0 && i >= 0) {
      i--;
      if(buffer[i] == class_event) {
	if(cnt == nturns)
	  stop = i;
	cnt--;
      }
    }
    start = i+1;
    
    for(i=start; i<=stop; i++) {
      if(((i+1) % (neq+2)) == 0)
	ofs << "\n";
      else if((i % (neq+2)) == 0)
	ofs << std::scientific << (buffer[i] - buffer[start]) << " ";
      else
	ofs << std::scientific << buffer[i] << " ";
    }
  }

  ofs.close();
  return true;
}

void ODESolver::HaltAtEquilibrium(bool halt) {
  halt_at_equilibrium = halt;
}

bool ODESolver::HaltsAtEquilibrium() const {
  return halt_at_equilibrium;
}

void ODESolver::HaltAtCycle(bool halt) {
  halt_at_cycle = halt;
}

bool ODESolver::HaltsAtCycle() const {
  return halt_at_cycle;
}

void ODESolver::SetClassificationEvent(int ce) {
  if(ce > 0)
    class_event = ce;
}

int ODESolver::GetClassificationEvent() const {
  return class_event;
}

int ODESolver::GetNumberOfTurns() const {
  return nturns;
}

realtype ODESolver::GetEquilibriumTolerance() const {
  return equilibrium_tolerance;
}

void ODESolver::SetEquilibriumTolerance(realtype tol) {
  if(tol > 0)
    equilibrium_tolerance = tol;
}

realtype ODESolver::GetCycleTolerance() const {
  return cycle_tolerance;
}

void ODESolver::SetCycleTolerance(realtype tol) {
  if(tol > 0)
    cycle_tolerance = tol;
}

N_Vector ODESolver::GetX() const {
  return x;
}

N_Vector ODESolver::GetXdot() const {
  return xdot;
}

N_Vector ODESolver::GetX0() const {
  return x0;
}

realtype* ODESolver::GetXEnd() const {
  if(rows == 0)
    return NULL;
  return buffer.get() + (rows-1)*cols + 1;
}

void ODESolver::SetX0(N_Vector X0, int n) {
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

void ODESolver::SetX0(realtype * X0, int n) {
  if(X0 != NULL && dynsys != NULL) {
    int stop;
    if(n == -1)
      stop = dynsys->GetDimension();
    else
      stop = n;
    printf("stop = %d\n", stop);
    for(int i=0; i<stop; i++)
      Ith(x0,i) = X0[i];
  }
}

inline void ODESolver::SetOrthonormalBaseIC() {
  int length = dynsys->GetOriginalDimension();
  for(int i=0; i<length; i++) {
    for(int j=0; j<length; j++)
      Ith(x0,length+i*length+j) = (i==j ? 1.0 : 0.0) ;
  }
}

boost::shared_array<realtype> ODESolver::GetLyapunovExponents() const {
  return lyapunov_exponents;
}

bool ODESolver::Setup() {
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
  flag = CVodeMalloc (cvode_mem, DynamicalSystem::RHSWrapper, 0.0, x, CV_SS, reltol, &abstol);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeMalloc.\n");
    return false;
  }
#endif
#ifdef CVODE26
  flag = CVodeInit (cvode_mem, DynamicalSystem::RHSWrapper, 0.0, x);
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
  flag = CVodeSetUserData (cvode_mem, static_cast<void *>(dynsys.get()));
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
    flag = CVDenseSetJacFn (cvode_mem, DynamicalSystem::JacobianWrapper, dynsys);
    if (flag != CV_SUCCESS) {
      fprintf (stderr, "Error on CVDenseSetJacFn.\n");
      return false;
    }
#endif
#ifdef CVODE26
    flag = CVDlsSetDenseJacFn (cvode_mem, DynamicalSystem::JacobianWrapper);
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

bool ODESolver::AllocateSolutionBuffer() {
  int lrows;

  switch(mode) {
  case LYAP:
    lrows = (int) ceil ((tfinal - ttran) / lyap_tstep);
    if(lrows < 0) lrows = 0;
    lrows += 2;
    break;
  case TRAJ:
    lrows = (int) ceil ((tfinal - ttran) / tstep);
    if(lrows < 0) lrows = 0;
    lrows += 2;
    break;
  case EVENTS:
    lrows = max_intersections + 2;
    break;
  case BOTH:
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
    buffer = boost::shared_array<realtype>(new realtype[bufsize]);
  } catch (std::bad_alloc&) {
    fprintf (stderr, "Not enough memory to allocate for the solution buffer...\n");
    return false;
  }

  return true;
}

void ODESolver::ResetInitialCondition() {
  t = t0;
  for(int i=0; i<neq; i++) {
    Ith(x,i) = Ith(x0,i);
  }
  ResetPositionInBuffer();
  StoreRecordInBuffer(START);
}

void ODESolver::ResetPositionInBuffer() {
  rows = 0;
}

void ODESolver::StoreRecordInBuffer(int lbl) {
  int actual_pos = rows*cols;
  buffer[actual_pos] = t;
  for(int i=0; i<neq; i++)
    buffer[actual_pos+i+1] = Ith(x,i);
  buffer[actual_pos+neq+1] = (realtype) lbl;
  rows++;
}

void ODESolver::ChangeCurrentLabel(int lbl) {
  buffer[rows*cols-1] = (realtype) lbl;
}

void ODESolver::SkipTransient(bool *equilibrium, bool *error) {
  int flag;
  realtype tout;
  *equilibrium = false;
  *error = false;

  tout = t0 + ttran;
  while (t < tout) {
    flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
    if (flag < 0 && flag != CV_TOO_MUCH_WORK && (flag != CV_ILL_INPUT || mode == TRAJ)) {
      *error = true;
      break;
    }
    
    /* 
     * if events are enabled, we give the dynamical system the possibility to
     * change its internal structure also during the transient evolution.
     */
    if ((mode == EVENTS || mode == BOTH) && flag == CV_ROOT_RETURN) {
      CVodeGetRootInfo (cvode_mem, events.get());
      if(dynsys->HasEventsConstraints()) {
	dynsys->EventsConstraints(t,x,events_constraints.get(),static_cast<void *>(dynsys.get()));
	/* 
	 * call ManageEvents so that the dynamical system can (optionally)
	 * change something in its internal structure. See for example
	 * PLL, a switch system.
	 */
	dynsys->ManageEvents(t,x,events.get(),events_constraints.get());
      }
      else {
	dynsys->ManageEvents(t,x,events.get());
      }
    }

    flag = CheckEquilibrium();
    if(flag != EQUIL_FALSE)
      *equilibrium = true;
    if(flag == EQUIL_BREAK)
      break;
  }
  // if it's an equilibrium point, the label is changed at the end...
  StoreRecordInBuffer(TRAN_END);
}

int ODESolver::CheckEquilibrium() {
  dynsys->RHS (t, x, xdot, static_cast<void *>(dynsys.get()));
  if(EuclideanDistance(neq,xdot) < equilibrium_tolerance) {
    nturns = 0;
    if(halt_at_equilibrium) {
      ChangeCurrentLabel(EQUIL);
      return EQUIL_BREAK;
    }
    return EQUIL_TRUE;
  }
  return EQUIL_FALSE;
}

int ODESolver::CheckCycle(int guess) {
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

void ODESolver::SetSolutionLength(int length) {
  if (length > 0) {
    rows = length / cols;
  }
}

realtype ODESolver::EuclideanDistance(int length, N_Vector x, N_Vector y) const {
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

bool ODESolver::GramSchmidtOrthonorm(realtype * znorm) const {
  int i, n;
  realtype *xx, *xnorm;
  n = dynsys->GetOriginalDimension();
  xx = new realtype[n*n];
  xnorm = new realtype[n*n];
  for(i=0; i<n*n; i++)
    xx[i] = Ith(x,n+i);
  GramSchmidtOrthonorm(xx,xnorm,znorm);
  for(i=0; i<n*n; i++)
    Ith(x,i+n) = xnorm[i];
  delete [] xnorm;
  delete [] xx;
  return true;
}

bool ODESolver::GramSchmidtOrthonorm(realtype * x, realtype * xnorm, realtype * znorm) const {
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

  delete [] tmp;
  return true;
} 

inline realtype ODESolver::Norm(int length, realtype * x) const {
  realtype norm = 0.0;
  for(int i = 0; i<length; i++)
    norm += x[i]*x[i];
  return sqrt(norm);
}  

inline realtype ODESolver::DotProduct(int length, realtype * x, realtype * y) const {
  realtype res = 0.0;
  for(int i=0; i<length; i++)
    res += x[i]*y[i];
  return res;
}

bool ODESolver::ResetCVode() {
  int flag;

  switch(mode) {
  case TRAJ:
  case LYAP:
#ifdef CVODE25
    flag = CVodeRootInit (cvode_mem, 0, NULL, dynsys);
#endif
#ifdef CVODE26
    flag = CVodeRootInit (cvode_mem, 0, NULL);
#endif
    break;
  case EVENTS:
  case BOTH:
#ifdef CVODE25
    flag = CVodeRootInit (cvode_mem, nev, DynamicalSystem::EventsWrapper, dynsys);
#endif
#ifdef CVODE26
    flag = CVodeRootInit (cvode_mem, nev, DynamicalSystem::EventsWrapper);
#endif
  }
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeRootInit: flag = %d.\n", flag);
    return false;
  }

#ifdef CVODE25
  flag = CVodeReInit (cvode_mem, DynamicalSystem::RHSWrapper, t0, x, CV_SS, reltol, &abstol);
#endif
#ifdef CVODE26
  flag = CVodeReInit (cvode_mem, t0, x);
#endif
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return false;
  }
  return true;
}

bool ODESolver::Solve() {
  bool retval;
  
  if ((mode == LYAP && lyap_tstep == 0) || tstep == 0) {
    fprintf (stderr, "Can't continue: the integration time step is 0.\n");
    throw "Zero integration time step";
  }

  if(mode == LYAP) {
    if(!dynsys->IsExtended()) {
      dynsys->Extend(true);
      int i, n = dynsys->GetOriginalDimension();
      realtype *x0_tmp = new realtype[n];
      for(i=0; i<n; i++)
	x0_tmp[i] = Ith(x0,i);
      SetDynamicalSystem(dynsys.get());
      for(i=0; i<n; i++)
	Ith(x0,i) = x0_tmp[i];
      delete [] x0_tmp;
    }
  }

  if(! setup)
    Setup();
  dynsys->Reset();

  switch (mode) {
  case TRAJ:
    retval = SolveWithoutEvents();
    break;
  case EVENTS:
  case BOTH:
    retval = SolveWithEvents();
    break;
  case LYAP:
    retval = SolveWithoutEventsLyapunov();
    break;
  }
  return retval;
}

bool ODESolver::SolveWithoutEventsLyapunov() {
  int i, j, flag;
  int n, N;
  realtype tout;
  bool eq, err;
  realtype *znorm, *cum;

  // dimension of the system
  n = dynsys->GetOriginalDimension();
  // dimension of the extended system
  N = n*(n+1);

  znorm = new realtype[n];
  cum = new realtype[n];
  for (i=0; i<n; i++)
    cum[i] = 0.0;

  lyapunov_exponents = boost::shared_array<realtype>(new realtype[n]);

  // set the initial condition
  for(i=n; i<N; i++)
    Ith(x0,i) = 0.0;

  // check whether the buffer is allocated or it needs to be
  // reallocated
  AllocateSolutionBuffer();
  // set initial conditions (this automatically saves the i.c.
  // in the first row of the solution buffer)
  ResetInitialCondition();

  // reset CVode
  if(! ResetCVode()) {
    delete [] znorm;
    delete [] cum;
    return false;
  }
  
  /************* TRANSIENT *****************/
  SkipTransient(&eq,&err);

  for(i=0; i<n; i++) {
    for(j=0; j<n; j++)
      Ith(x,(i+1)*n+j) = (i==j ? 1. : 0.);
  }
  rows--;
  StoreRecordInBuffer(TRAN_END);
  
  /************* TRAJECTORY *****************/
  if(!err) {
    tout = t0 + ttran + lyap_tstep;
    while (tout < t0+tfinal+lyap_tstep) {
      // the integrator must be re-initialised every time the state is modified
#ifdef CVODE25
      flag = CVodeReInit (cvode_mem, DynamicalSystem::RHSWrapper, t, x, CV_SS, reltol, &abstol);
#endif
#ifdef CVODE26
      flag = CVodeReInit (cvode_mem, t, x);
#endif
      if (flag != CV_SUCCESS) {
	fprintf (stderr, "Error on CVodeReInit.\n");
	delete [] znorm;
	delete [] cum;
	return false;
      }

      flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
      StoreRecordInBuffer(REGUL);
      
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

      GramSchmidtOrthonorm(znorm);
      for(i=0; i<n; i++)
	cum[i] += log(std::max(znorm[i],1e-12))/log(2.0);

      tout += lyap_tstep;
    }
    for(i=0; i<n; i++){
      lyapunov_exponents[i] = cum[i]/t;
    }
  }

  delete [] znorm;
  delete [] cum;

  /* if the integrator stopped because of an error, we change the label of the last row in
   * the integration buffer and stop the integration procedure */
  if (err) {
    ChangeCurrentLabel(ERROR);
    return false;
  }

  return true;
}

/*** DEPRECATED ***/
bool ODESolver::SolveLyapunov() {
  fprintf(stderr, "ODESolver::SolveLyapunov>> This function is deprecated.\n");
  SetIntegrationMode(TRAJ);

  int i;
  realtype * temp_x0 = new realtype[dynsys->GetOriginalDimension()];
  for(i=0; i<dynsys->GetOriginalDimension(); i++)
    temp_x0[i] = Ith(x0,i);
  realtype temp_ttran = ttran;
  realtype tend = tfinal;
  realtype t_ = 0.0;
  tfinal = ttran;
  SolveWithoutEvents();
  
  realtype * new_x0 = new realtype [dynsys->GetDimension()];
  for (i=0; i<dynsys->GetDimension(); i++)
    new_x0[i] = Ith(x,i);
  
  dynsys->Extend(true);
  SetDynamicalSystem(dynsys.get());
  
  int N = dynsys->GetDimension();
  int n = dynsys->GetOriginalDimension();
  realtype * x_ = new realtype[N];
  realtype * xnorm = new realtype[n*n];
  realtype * znorm = new realtype[n];
  realtype * cum = new realtype[n];
  lyapunov_exponents = boost::shared_array<realtype>(new realtype[n]);
  for (i=0; i<n; i++)
    cum[i] = 0.0;
  
  SetX0(new_x0,n);
  delete [] new_x0;
  SetOrthonormalBaseIC();
  
  tfinal = lyap_tstep;
  ttran = 0.0;
  
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
      cum[i] += log(std::max(znorm[i],1e-12))/log(2.0);
    SetX0(x_);
  }
  
  for(i=0; i<n; i++)
    lyapunov_exponents[i] = cum[i]/t_;
  
  ttran = temp_ttran;
  tfinal = tend;
  dynsys->Extend(false);
  SetDynamicalSystem(dynsys.get());
  Setup();
  SetX0(temp_x0);
  
  delete [] temp_x0;
  delete [] x_;
  delete [] xnorm;
  delete [] znorm;
  delete [] cum;

  SetIntegrationMode(LYAP);
  return true;
}

bool ODESolver::SolveWithoutEvents() {
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
  
  static bool print = false;

  /************* TRAJECTORY *****************/
  if(!(eq && halt_at_equilibrium) && !err) {
    tout = t+tstep;
    while (tout < t0+tfinal+tstep) {
      flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
      StoreRecordInBuffer(REGUL);
      
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
  print = true;

  /* if the integrator stopped because of an error, we change the label of the last row in
   * the integration buffer and stop the integration procedure */
  if (err) {
    ChangeCurrentLabel(ERROR);
    return false;
  }
  /* if the point is an equilibrium, we change the label of the last row in
   * the integration buffer and stop the integration procedure */
  if (eq || nturns == 0)
    ChangeCurrentLabel(EQUIL);

  return true;
}

bool ODESolver::SolveWithEvents() {
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
    if(mode == EVENTS)
      tout = t0+tfinal;
    else if(mode == BOTH)
      tout = t+tstep;
    intersections = 0;
    class_inters = 0;
    nturns_guess = -1;
    cycle = false;
    while (intersections < max_intersections && t < t0+tfinal) {
      flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
      if (flag < 0 && flag != CV_TOO_MUCH_WORK && flag != CV_ILL_INPUT) {
	/*
	 * this is because the flag CV_ILL_INPUT is returned when two
	 * events are found at a very small time interval, which is the
	 * case if at least one of the Poincare' sections corresponds to an
	 * extremum of a state variables, which is a condition always
	 * verified when the trajectory is on an equilibrium
	 */
	StoreRecordInBuffer(ERROR);
	err = true;
	break;
      }

      eq_flag = CheckEquilibrium();
      if(eq_flag == EQUIL_BREAK)
	break;
      else if(eq_flag == EQUIL_TRUE)
	eq = true;

      if (flag == CV_SUCCESS && mode == BOTH) {
	StoreRecordInBuffer(REGUL);
	tout += tstep;
      }
      else if (flag == CV_ROOT_RETURN && !eq) {
	CVodeGetRootInfo (cvode_mem, events.get());
	if(dynsys->HasEventsConstraints()) {
	  dynsys->EventsConstraints(t,x,events_constraints.get(),static_cast<void *>(dynsys.get()));
	  /* 
	   * call ManageEvents so that the dynamical system can (optionally)
	   * change something in its internal structure. See for example
	   * PLL, a switch system.
	   */
	  dynsys->ManageEvents(t,x,events.get(),events_constraints.get());
	}
	else {
	  dynsys->ManageEvents(t,x,events.get());
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
    ChangeCurrentLabel(ERROR);
    return false;
  }
  
  /* if we weren't able to detect a periodicity, there's a good chance 
   * that the system be chaotic */
  if (nturns == -1)
    nturns = max_intersections;
  
  /* if the point is an equilibrium, we change the label of the last row in
   * the integration buffer */
  if (eq || nturns == 0)
    ChangeCurrentLabel(EQUIL);
    
  return true;
}

ODESolver* ODESolver::Clone() const {
  return new ODESolver(*this);
}

} // namespace bal

