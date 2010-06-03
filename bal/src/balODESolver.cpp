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

#include "balODESolver.h"

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
	t = 0.0;
	tstep = STEP;
	ttran = T_TRAN;
	tfinal = T_END;
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
}

balODESolver::balODESolver(const balODESolver& solver) {
	buffer = NULL;
	delete_buffer = false;
	rows = 0;
	
	SetDynamicalSystem((solver.GetDynamicalSystem())->Copy());
	
	for(int i=0; i<dynsys->GetDimension(); i++)
		Ith(x0,i) = Ith(solver.x0,i);
	
	reltol = solver.GetRelativeTolerance();
	abstol = solver.GetAbsoluteTolerance();
	mode = solver.GetIntegrationMode();
	t = 0.0;
	tstep = solver.GetTimeStep();
	ttran = solver.GetTransientDuration();
	tfinal = solver.GetFinalTime();
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
	if(x != NULL) {
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
		if(x != NULL) {
			N_VDestroy_Serial(x);
			N_VDestroy_Serial(xdot);
			N_VDestroy_Serial(x0);
			N_VDestroy_Serial(x_inters);
		}
		x = N_VNew_Serial(neq);
		xdot = N_VNew_Serial(neq);
		x0 = N_VNew_Serial(neq);
		x_inters = N_VNew_Serial(neq);
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

void balODESolver::SetX0(N_Vector X0) {
	if(X0 != NULL && dynsys != NULL) {
		for(int i=0; i<dynsys->GetDimension(); i++)
			Ith(x0,i) = Ith(X0,i);
	}
}

void balODESolver::SetX0(realtype * X0) {
	if(X0 != NULL && dynsys != NULL) {
		for(int i=0; i<dynsys->GetDimension(); i++)
			Ith(x0,i) = X0[i];
	}
}

bool balODESolver::Setup() {
	int flag;

  if (cvode_mem != NULL) {
		CVodeFree (&cvode_mem);
	}
	
  cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
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
		}
		delete_buffer = false;
		switch(mode) {
			case balTRAJ:
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

bool balODESolver::Solve() {
	if(! setup)
		Setup();
	dynsys->Reset();
  switch (mode) {
    case balTRAJ:
      return SolveTrajectory();
    case balEVENTS:
			return SolveIntersections();
    case balBOTH:
			return SolveTrajectoryAndIntersections();
  }
	return false;
}

bool balODESolver::SolveTrajectory() {
  int flag;
  realtype tout;
	bool eq, err;

	// check whether the buffer is allocated or it needs to be
	// reallocated
	AllocateSolutionBuffer();
	// set initial conditions (this automatically saves the i.c.
	// in the first row of the solution buffer)
	ResetInitialCondition();

  // CVode initialization
#ifdef CVODE25
  flag = CVodeRootInit (cvode_mem, 0, NULL, dynsys);
  if (flag != CV_SUCCESS) {
		fprintf (stderr, "Error on CVodeRootInit.\n");
    return 0;
  }
  flag = CVodeReInit (cvode_mem, balDynamicalSystem::RHSWrapper, 0.0, x, CV_SS, reltol, &abstol);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return 0;
  }
#endif
#ifdef CVODE26
  flag = CVodeRootInit (cvode_mem, 0, NULL);
  if (flag != CV_SUCCESS) {
		fprintf (stderr, "Error on CVodeRootInit.\n");
    return 0;
  }
  flag = CVodeReInit (cvode_mem, 0.0, x);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return 0;
  }
#endif

  /************* TRANSIENT *****************/
	nturns = -1;
  tout = ttran;
	eq = false;
	err = false;
  while (t < tout) {
    flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
    if (flag < 0 && flag != CV_TOO_MUCH_WORK) {
			err = true;
	  	break;
		}
		dynsys->RHS(t, x, xdot, params);
		if(EuclideanDistanceFromOrigin(xdot,neq) < equilibrium_tolerance) {
			nturns = 0;
			if(halt_at_equilibrium) {
				eq = true;
				break;
			}
		}
  }
	StoreRecordInBuffer(balTRAN_END);

  /************* TRAJECTORY *****************/
	if(!eq && !err) {
		tout += tstep;
		while (t < tfinal) {
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
				// NON DOVREBBE RILEVARE GLI EVENTI IN MODALITA' balTRAJ
				if(flag == CV_ILL_INPUT) {
					dynsys->RHS (t, x, xdot, params);
					if(EuclideanDistanceFromOrigin(xdot,neq) < equilibrium_tolerance) {
						nturns = 0;
						ChangeCurrentLabel(balEQUIL);
						if(halt_at_equilibrium)
							break;
					}
				}
				/* this is done if the error is not CV_ILL_INPUT */
				err = true;
		  	break;
			}

			dynsys->RHS (t, x, xdot, params);
			if(EuclideanDistanceFromOrigin(xdot,neq) < equilibrium_tolerance) {
				nturns = 0;
				ChangeCurrentLabel(balEQUIL);
				if(halt_at_equilibrium)
					break;
			}
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
  if (eq || nturns == 0) {
		ChangeCurrentLabel(balEQUIL);
    return true;
  }

	return true;
}

bool balODESolver::SolveIntersections() {
  int flag, i, j;
  int intersections; // the total number of intersections with Poincare' sections
	int inters1; // the number of intersections with the first Poincare' section (used to detect cycles)
	int nturns_guess;
  realtype tout;
	bool eq, err, cycle, restart;

	// check whether the buffer is allocated or it needs to be
	// reallocated
	AllocateSolutionBuffer();
	// set initial conditions
	ResetInitialCondition();

#ifdef CVODE25
	flag = CVodeRootInit (cvode_mem, nev, balDynamicalSystem::EventsWrapper, dynsys);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeRootInit.\n");
    return false;
  }
  flag = CVodeReInit (cvode_mem, balDynamicalSystem::RHSWrapper, 0.0, x, CV_SS, reltol, &abstol);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return 0;
  }
#endif
#ifdef CVODE26
  flag = CVodeRootInit (cvode_mem, nev, balDynamicalSystem::EventsWrapper);
  if (flag != CV_SUCCESS) {
		fprintf (stderr, "Error on CVodeRootInit.\n");
    return 0;
  }
  flag = CVodeReInit (cvode_mem, 0.0, x);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return 0;
  }
#endif

  /************* TRANSIENT *****************/
	nturns = -1;
  tout = ttran;
  eq = false;
	err = false;
  t = 0.0;
  while (t < tout) {
		flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
    if (flag < 0 && flag != CV_TOO_MUCH_WORK && flag != CV_ILL_INPUT) {
			err = true;
			break;
		}
		dynsys->RHS (t, x, xdot, params);
		if(EuclideanDistanceFromOrigin(xdot,neq) < equilibrium_tolerance) {
			nturns = 0;
			if(halt_at_equilibrium) {
				eq = true;
				break;
			}
		}
  }
	// if it's an equilibrium point, the label is changed at the end...
	StoreRecordInBuffer(balTRAN_END);
  
  /************* INTERSECTIONS ************/
	restart = false;
	if(!eq && !err) {
		tout = tfinal;
  	intersections = 0;
		inters1 = 0;
		nturns_guess = -1;
		cycle = false;
  	while (intersections < max_intersections && t < tout) {
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

			dynsys->RHS (t, x, xdot, params);
			if(EuclideanDistanceFromOrigin(xdot,neq) < equilibrium_tolerance) {
				nturns = 0;
				if(halt_at_equilibrium) {
					StoreRecordInBuffer(balEQUIL);
					break;
				}
			}

			if (flag == CV_ROOT_RETURN) {
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
							if(i == 0) {
								inters1++;
								if(inters1 == 1) {
									// save the first intersection
									//printf("first intersection saved @ t = %f: [ ", t);
									for(j=0; j<neq; j++) {
										Ith(x_inters,j) = Ith(x,j);
										//printf("%e ", Ith(x_inters,j));
									}
									//printf("]\n");
								}
								else {
									if(nturns_guess == -1) {
										if(EuclideanDistance(x,x_inters,neq) < cycle_tolerance) {
											// guess the number of turns
											nturns_guess = inters1-1;
											//printf("the guess @ t = %f is %d\n", t, nturns_guess);
										}
										else if(intersections >= round(max_intersections/2) && !restart) {
											// restart to keep track of intersections if we reach
											// half of the allowed intersections and we haven't
											// made a guess yet
											inters1 = 0;
											restart = true;
											//printf("restart @ t = %f\n", t);
										}
									}
									else if(inters1 == 2*nturns_guess + 1) {
										if(EuclideanDistance(x,x_inters,neq) < cycle_tolerance) {
											nturns = nturns_guess;
											//printf("The number of turns is %d\n", nturns);
											if(halt_at_cycle) {
												cycle = true;
												break;
											}
										}
										else {
											// restart to keep track of intersections and see if
											// we can detect a new cycle
											nturns_guess = -1;
											inters1 = 0;
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
  /* if the point is an equilibrium, we change the label of the last row in
	 * the integration buffer */
  if (eq || nturns == 0) {
		ChangeCurrentLabel(balEQUIL);
    return true;
  }

	if (nturns == -1)
		nturns = max_intersections;	// if we weren't able to detect a periodicity, there's a good chance that the system be chaotic

	return true;
}

bool balODESolver::SolveTrajectoryAndIntersections() {
  int flag, i, j;
  int intersections; // the total number of intersections with Poincare' sections
	int inters1; // the number of intersections with the first Poincare' section (used to detect cycles)
	int nturns_guess;
  realtype tout;
	bool eq, err, cycle, restart;

	// check whether the buffer is allocated or it needs to be
	// reallocated
	AllocateSolutionBuffer();
	// set initial conditions
	ResetInitialCondition();

#ifdef CVODE25
	flag = CVodeRootInit (cvode_mem, nev, balDynamicalSystem::EventsWrapper, dynsys);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeRootInit.\n");
    return 0;
  }
  flag = CVodeReInit (cvode_mem, balDynamicalSystem::RHSWrapper, 0.0, x, CV_SS, reltol, &abstol);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return 0;
  }
#endif
#ifdef CVODE26
  flag = CVodeRootInit (cvode_mem, nev, balDynamicalSystem::EventsWrapper);
  if (flag != CV_SUCCESS) {
		fprintf (stderr, "Error on CVodeRootInit.\n");
    return 0;
  }
  flag = CVodeReInit (cvode_mem, 0.0, x);
  if (flag != CV_SUCCESS) {
    fprintf (stderr, "Error on CVodeReInit.\n");
    return 0;
  }
#endif

  /************* TRANSIENT *****************/
	nturns = -1;
  tout = ttran;
	eq = false;
	err = false;
	t = 0.0;
  while (t < tout) {
    flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);
    if (flag < 0 && flag != CV_TOO_MUCH_WORK && flag != CV_ILL_INPUT) {
			err = true;
			break;
		}
		dynsys->RHS (t, x, xdot, params);
		if(EuclideanDistanceFromOrigin(xdot,neq) < equilibrium_tolerance) {
			nturns = 0;
			if(halt_at_equilibrium) {
				eq = true;
				break;
			}
		}
  }
	// if it's an equilibrium point, the label is changed at the end...
	StoreRecordInBuffer(balTRAN_END);

  /************* TRAJECTORY AND INTERSECTIONS ************/
	restart = false;
	if(!eq && !err) {
	  tout += tstep;
	  intersections = 0;
		inters1 = 0;
		nturns_guess = -1;
		cycle = false;
	  while (intersections < max_intersections && t < tfinal) {
	    flag = CVode (cvode_mem, tout, x, &t, CV_NORMAL);

			if (flag < 0 && flag != CV_TOO_MUCH_WORK && flag != CV_ILL_INPUT) {
				/*
				 * this is because the flag CV_ILL_INPUT is returned when two
				 * events are found at a very small time interval, which is the
				 * case if at least one of the Poincare' sections corresponds to
				 * extrema of a state variables, which is a condition always
				 * verified when the trajectory is on an equilibrium
				 */
				StoreRecordInBuffer(balERROR);
				err = true;
				break;
			}
		
			dynsys->RHS (t, x, xdot, params);
			if(EuclideanDistanceFromOrigin(xdot,neq) < equilibrium_tolerance) {
				nturns = 0;
				if(halt_at_equilibrium) {
					StoreRecordInBuffer(balEQUIL);
					break;
				}
			}

			if (flag == CV_SUCCESS) {
				StoreRecordInBuffer(balSTEP);
			  tout += tstep;
			}
	    else if (flag == CV_ROOT_RETURN) {
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
							if(i == 0) {
								inters1++;
								if(inters1 == 1) {
									// save the first intersection
									//printf("first intersection saved @ t = %f: [ ", t);
									for(j=0; j<neq; j++) {
										Ith(x_inters,j) = Ith(x,j);
										//printf("%e ", Ith(x_inters,j));
									}
									//printf("]\n");
								}
								else {
									if(nturns_guess == -1) {
										//printf("the distance @ t = %f is %e\n", t, EuclideanDistance(x,x_inters,neq));
										if(EuclideanDistance(x,x_inters,neq) < cycle_tolerance) {
											// guess the number of turns
											nturns_guess = inters1-1;
											//printf("the guess @ t = %f is %d\n", t, nturns_guess);
										}
										else if(intersections >= round(max_intersections/2) && !restart) {
											// restart to keep track of intersections if we reach
											// half of the allowed intersections and we haven't
											// made a guess yet
											inters1 = 0;
											restart = true;
											//printf("restart @ t = %f\n", t);
										}
									}
									else if(inters1 == 2*nturns_guess + 1) {
										if(EuclideanDistance(x,x_inters,neq) < cycle_tolerance) {
											nturns = nturns_guess;
											//printf("The number of turns @ t = %f is %d\n", t, nturns);
											if(halt_at_cycle) {
												cycle = true;
												break;
											}
										}
										else {
											// restart to keep track of intersections and see if
											// we can detect a new cycle
											nturns_guess = -1;
											inters1 = 0;
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
  /* if the point is an equilibrium, we change the label of the last row in
	 * the integration buffer and stop the integration procedure */
  if (eq || nturns == 0) {
		ChangeCurrentLabel(balEQUIL);
    return true;
  }

	if (nturns == -1)
		nturns = max_intersections;	// if we weren't able to detect a periodicity, there's a good chance that the system be chaotic

	return true;
}

void balODESolver::SetSolutionLength(int length) {
  if (length > 0) {
		rows = length / cols;
  }
}

realtype balODESolver::EuclideanDistanceFromOrigin(N_Vector x, int length) {
	realtype dst = 0.0;
	for(int i=0; i<length; i++) {
		dst += Ith(x,i) * Ith(x,i);
	}
	return sqrt(dst);
}

realtype balODESolver::EuclideanDistance(N_Vector x, N_Vector y, int length) {
	realtype dst = 0.0;
	for(int i=0; i<length; i++) {
		dst += (Ith(x,i) - Ith(y,i)) * (Ith(x,i) - Ith(y,i));
	}
	return sqrt(dst);
}

