/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balODESolver.h
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
 * \file balODESolver.h
 * \brief Definition of the class balODESolver
 */


#ifndef _BALODESOLVER_
#define _BALODESOLVER_

#include "balObject.h"
#include "balCommon.h"
#include "balParameters.h"
#include "balDynamicalSystem.h"
#include "balSolution.h"

#include <cstdio>
#include <cstdlib>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#ifdef CVODE25
#include <sundials/sundials_direct.h>
#endif

/** Default relative tolerance */
#define RTOL (1.0E-7)

/** Default absolute tolerance */
#define ATOL (1.0E-10)

/** Default integration step   */
#define STEP (1.0E-1)

/** Default final integration time */
#define T_END (1.0E3)

/** Default transient time */
#define T_TRAN (0.0)

/** Default number of intersection with the Poincare' section */
#define DEFAULT_INTERSECTIONS (100)

/** Maximum number of steps (for CVode) */
#define MAX_NUM_STEPS ((int)1E6)

/** Default equilibrium tolerance */
#define EQUIL_TOL (1.0E-5)

/** 
 * Default cycle detection tolerance: this value is two orders
 * of magnitude bigger than the equilibrium tolerance because the
 * check for the detection of a closed cycle is performed twice.
 */
#define CYCLE_TOL (1.0E-3)

/** Error file for CVode */
#define ERROR_FILE "bal.log"

#define DUMPBUFFER(fid)					\
  {							\
    int i, j;						\
    for(i=0; i<rows; i++) {				\
      for(j=0; j<cols; j++)				\
	fprintf(fid, "%e ", buffer[i*cols+j]);		\
      fprintf(fid, "\n");				\
    }							\
  }



/** The integration mode of the solver. */
typedef enum
  {
    /** Records the trajectory of the system */
    balTRAJ = 1,
    /** Records the events, i.e. the crossings of Poincare' sections */
    balEVENTS,
    /** Records both trajectory and events */
    balBOTH,
    /** Compute Lyapunov spectrum */
    balLYAP
  } integration_mode;


/**
 * The labels associated to every row in the solution matrix. Events have
 * labels that start from 1 and that correspond to the order in which the
 * corresponding Poincare' section has been loaded into the solver.
 */
typedef enum
  {
    /** Integration has terminated due to an error in CVode */
    balERROR = -10,
    /** Integration has terminated because the system has reached an
     * equilibrium point */
    balEQUIL = -3,
    /** Initial point */
    balSTART = -2,
    /** Point at the end of the transient evolution */
    balTRAN_END = -1,
    /** A normal trajectory point */
    balSTEP = 0
  } label;

enum { EQUIL_FALSE, EQUIL_TRUE, EQUIL_BREAK };
enum { CYCLE_FALSE, CYCLE_TRUE, CYCLE_BREAK };

/**
 * \class balODESolver
 * \brief Class for integration of dynamical systems
 * 
 * \sa balDynamicalSystem balBifurcationParameters balSolution
 */
class balODESolver : public balObject {

 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  /** Creates a new balODESolver. */
  static balODESolver * Create();
  /** Creates a copy of balODESolver. */
  static balODESolver * Copy (balODESolver * solver);
  /** Destroys a balODESolver. */
  virtual void Destroy();
  
  realtype * GetBuffer() const;
  int GetBufferSize() const;
  void GetBufferSize(int * r, int * c) const;
  
  void SetDynamicalSystem(balDynamicalSystem * ds);
  balDynamicalSystem * GetDynamicalSystem() const;
  
  void SetDynamicalSystemParameters(balParameters * par);
  
  balSolution * GetSolution() const;
  
  realtype GetInitialTime() const;
  void SetInitialTime(realtype T0);
  realtype GetTransientDuration () const;
  void SetTransientDuration (realtype tran);
  realtype GetFinalTime () const;
  void SetFinalTime (realtype final);
  realtype GetTimeStep () const;
  void SetTimeStep (realtype step);
  realtype GetLyapunovTimeStep () const;
  void SetLyapunovTimeStep (realtype tstep);
  realtype GetRelativeTolerance () const;
  void SetRelativeTolerance (realtype rtol);
  realtype GetAbsoluteTolerance () const;
  void SetAbsoluteTolerance (realtype atol);
  integration_mode GetIntegrationMode () const;
  void SetIntegrationMode (integration_mode m);
  void IsStiff(bool stiffness);
  int GetMaxNumberOfIntersections() const;
  void SetMaxNumberOfIntersections(int intersections);
  void HaltAtEquilibrium(bool halt);
  bool HaltsAtEquilibrium() const;
  void HaltAtCycle(bool halt);
  bool HaltsAtCycle() const;
  void SetClassificationEvent(int ce);
  int GetClassificationEvent() const;
  int GetNumberOfTurns() const;
  realtype GetEquilibriumTolerance() const;
  void SetEquilibriumTolerance(realtype tol);
  realtype GetCycleTolerance() const;
  void SetCycleTolerance(realtype tol);

  realtype * GetLyapunovExponents() const;
  
  N_Vector GetX() const;
  N_Vector GetXdot() const;
  N_Vector GetX0() const;
  realtype * GetXEnd() const;
  void SetX0(N_Vector X0, int n = -1);
  void SetX0(realtype * X0, int n = -1);
  
  bool SaveOrbit(const char *filename) const;

  virtual bool Setup();
  virtual bool Solve();
  
 protected:
  balODESolver();
  balODESolver (const balODESolver & solver);
  virtual ~balODESolver();
  bool AllocateSolutionBuffer();
  void StoreRecordInBuffer(int lbl);
  void ResetPositionInBuffer();
  void ResetInitialCondition();
  bool ResetCVode();
  void ChangeCurrentLabel(int lbl);
  void SkipTransient(bool *equilibrium, bool *error);
  bool SolveWithoutEvents();
  bool SolveWithoutEventsLyapunov();
  bool SolveWithEvents();
  bool SolveLyapunov();
  void SetSolutionLength(int length);
  int CheckEquilibrium();
  int CheckCycle(int guess);
  
  realtype EuclideanDistance(int length, N_Vector x, N_Vector y = NULL) const;
  inline realtype Norm(int length, realtype* x) const;
  inline realtype DotProduct(int length, realtype* x, realtype* y) const;
  bool GramSchmidtOrthonorm(realtype* x, realtype* xnorm, realtype* znorm) const;
  bool GramSchmidtOrthonorm(realtype* znorm) const;
  inline void SetOrthonormalBaseIC();
  
 private:
  /** Number of equations */
  int neq;
  /** Number of parameters */
  int npar;
  /** Number of events */
  int nev;
  
  /** Buffer for integration data */
  realtype * buffer;
  bool delete_buffer;
  /** The total size of the buffer */
  unsigned long bufsize;
  /** The (virtual) dimensions of the buffer */
  unsigned long rows;
  int cols;
  /** The actual position in the buffer */
  int position_in_buffer;
  
  /** Initial integration time */
  realtype t0;
  /** Current time */
  realtype t;
  /** Transient duration */
  realtype ttran; // SG
  /** Integration duration */
  realtype tfinal; // SG
  /** Integration time step */
  realtype tstep; // SG
  /** Lypunov exponents calculus tstep*/
  realtype lyap_tstep;
	
  realtype * lyapunov_exponents;
  bool delete_lyapunov_exponents;
	
  /** Relative tolerance */
  realtype reltol; // SG
  /** Absolute tolerance */
  realtype abstol; // SG
  
  /** Integration mode (step by step or event-driven */
  integration_mode mode; // SG
  /** Indicates whether we are solving a stiff or non-stiff system (the former one is the default */
  bool stiff; // S

  /** The memory for CVode */
  void *cvode_mem;
  
  /** The file for error logging */
  FILE *errfp;
  
  /** The maximum number of intersections */
  int max_intersections; // SG
  /** Event location array (each element has value 1 when the corresponding
   * event has been located */
  int *events;
  bool delete_events;
  /** Constraints on the location of events */
  int *events_constraints;
  bool delete_events_constraints;

  /** Initial state */
  N_Vector x0; // SG
  /** Current state */
  N_Vector x; // G
  /** Current derivative */
  N_Vector xdot; // G
  /** Coordinates of the first intersection with the Poincare' section
   * (used to detect cycles) */
  N_Vector x_inters;
  bool nvectors_allocated;

  /** Indicates whether setup has been performed or not */
  bool setup;
  
  /**
   * Indicates whether integration should stop if an equilibrium point is
   * reached.
   */
  bool halt_at_equilibrium; // SG
  
  /** Indicates whether integration should stop if a cycle is detected. */
  bool halt_at_cycle;
  /** The number of turns of the detected cycle. */
  int nturns;

  /** The event that should be checked to detect a cycle */
  int class_event; // SG

  /** The tolerance to say that a point is an equilibrium */
  realtype equilibrium_tolerance; // SG
  /** The tolerance to say that a cycle has been detected */
  realtype cycle_tolerance; // SG
  
  /** Parameters of the system */
  balParameters * params;
  /** The dynamical system that has to be integrated */
  balDynamicalSystem * dynsys;
};


#endif
