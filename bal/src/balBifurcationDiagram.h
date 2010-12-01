/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balBifurcationDiagram.h
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
 * \file balBifurcationDiagram.h
 * \brief Definition of the class balBifurcationDiagram.
 */

#ifndef _BALBIFURCATIONDIAGRAM_
#define _BALBIFURCATIONDIAGRAM_

#include <csignal>

#include "balODESolver.h"
#include "balDynamicalSystem.h"
#include "balParameters.h"
#include "balBifurcationParameters.h"
#include "balLogger.h"
#include "balObject.h"
#include "balCommon.h"

#include <cstdlib>
#include <cstdio>
#include <list>
#include <boost/thread.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>
#include <boost/ref.hpp>

using std::list;

#define DEBUG

void ResetColours(int d);

enum { balPARAMS, balIC };

/**
 * \class balSummaryEntry
 * \brief Object used to store an entry of the summary list.
 */
class balSummaryEntry : public balObject {
 public:
  balSummaryEntry(balSolution *sol, int mode = balPARAMS);
  virtual ~balSummaryEntry();
  int GetN() const;
  int GetID() const;
  double* GetData() const;

 private:
  double* data;
  int n;
  int id;
};

bool CompareBalSummaryEntry(balSummaryEntry *entry1, balSummaryEntry *entry2);

/**
 * \class balBifurcationDiagram
 * \brief Class to calculate brute-force bifurcation diagrams.
 *
 * This class is used to compute brute-force bifurcation diagrams: it
 * manages 4 different kind of objects: a balDynamicalSystem that defines
 * the RHS function of the system to integrate, a balODESolver that
 * performs the actual integration, a balLogger to store data to file and a
 * balBifurcationParameters that is mainly used to iterate over all the
 * combinations of parameters that the user is interested to simulate.
 *
 * \sa balDynamicalSystem balBifurcationParameters balODESolver
 */
class balBifurcationDiagram : public balObject {
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const;
  /** Creates a new balBifurcationDiagram. */
  static balBifurcationDiagram * Create();
  /** Destroys a balBifurcationDiagram. */
  virtual void Destroy();

  /**
   * Sets the dynamical system to integrate. balBifurcationDiagram
   * assumes that the dynamical contain an instance of balBifurcationParameters,
   * instead of the simpler balParameters.
   * @param sys A dynamical system: any instance of a class inherited
   * from balDynamicalSystem.
   */
  void SetDynamicalSystem(balDynamicalSystem * sys);

  /**
   * Gets the dynamical system to integrate.
   * @return The dynamical system used in the computation of the
   * bifurcation diagram.
   */
  balDynamicalSystem * GetDynamicalSystem() const;

  /**
   * Sets the logger used for saving integration data to file. A logger
   * is automatically instantiated when a balBifurcationDiagram is
   * created. This method should be used only if the user wants to use a
   * different kind of logger (such as one, for example, that saves data
   * in a particular format). The default logger uses H5 files.
   * @param log An instance of one of the classes inherited by balLogger.
   */
  void SetLogger(balLogger * log);

  /**
   * @return The logger used for saving data to file.
   */
  balLogger * GetLogger() const;

  /**
   * Sets the ODESolver used to integrate the dynamical system. An ODE
   * solver is automatically instantiated when a balBifurcationDiagram is
   * created. This method should be called only if the user has developed
   * their own ODE solver.
   * \param sol An instance of an ODE solver.
   */
  void SetODESolver(balODESolver * sol);

  /**
   * This method returns a pointer to the ODE solver used to integrate
   * the system. It is useful to set parameters of the ODE solver.
   * @return The ODE solver used to integrate the system.
   */
  balODESolver * GetODESolver() const;

  /**
   * Sets the name of the file where data will be saved.
   * @param filename File where the bifurcation diagram will be saved.
   * @param compress Flag that tells wheter the file should be compressed.
   */
  void SetFilename(const char * filename, bool compress = false);

  /**
   * @return The name of the file where data is saved.
   */
  const char * GetFilename();

  /**
   * Performs the actual computation of the brute-force bifurcation
   * diagram.
   */
  void ComputeDiagram();

  /**
   * Asks whether each new integration is restarted from the original
   * initial conditions.
   */
  bool RestartsFromX0() const;

  /**
   * Sets whether each new integration should restart or not from the
   * original initial conditions.
   */
  void RestartFromX0(bool restart);

  /**
   * Sets the number of parallel threads launched to compute the
   * bifurcation diagrams. The ideal value for _nthreads is equal to the
   * number of cores of the processor. The default value is 2 (who
   * doesn't own a dual-core these days?!?).
   */
  void SetNumberOfThreads(int _nthreads);

  /**
   * Gets the number of parallel threads launched to compute the
   * bifurcation diagrams.
   */
  int GetNumberOfThreads() const;

  bool SaveSummaryData(const char *filename) const;
  double** GetSummaryData(int *size = NULL) const;

  bool SetMode(int _mode);
  int GetMode() const;
  void SetInitialConditions(int nx0, double **x0);

 protected:
  balBifurcationDiagram();
  virtual ~balBifurcationDiagram();

 private:

  //void ComputeDiagramSingleThread();
  void ComputeDiagramMultiThread();
  void IntegrateAndEnqueue(balODESolver *sol, int solutionId);
  double* BuildSummaryEntry(balSolution *sol);

  /** The ODE solver used to integrate the system */
  balODESolver * solver;
  /** The dynamical system to integrate */
  balDynamicalSystem * system;
  /** The parameters of the dynamical system */
  balParameters * parameters;
  /**
   * The object used to save data to a file: by default H5 file logging
   * is used, i.e., logger is an instance of the balH5Logger class:
   * if the user wants to use another logger, they should provide it
   * by using the SetLogger member.
   */
  balLogger * logger;
  /** Tells whether the logger should be deleted when a new one is set. */
  bool destroy_logger;
  /** Tells whether the solver should be deleted when a new one is set. */
  bool destroy_solver;

  /** The number of dimensions of the system. */
  int ndim;
  /** The total number of parameters of the system. */
  int npar;
  /** The number of parameters to be varied for the bifurcation diagram. */
  int n_var_par;
  /** The mode of computation of the diagram. */
  int mode;

  int nX0;
  double **X0;

  /** Tells wheter each integration should be restarted from the initial
   * condition set by the user (true, default value) or from the final value of the
   * previous integration (false) */
  bool restart_from_x0;

  /*** multithreading stuff ***/
  /**
   * A list containing the results of the numerical integrations: when the list
   * is full, the threads that integrate stop and another thread saves data to file
   */
  list<balSolution *> *solutions;
  /**
   * A list containing the summary of the bifurcation diagram in
   * terms of number of turns of the solution. Each entry contains
   * (n+1) values, where n is the number of parameters of the system.
   * The first n values are the values of the parameters, the last value
   * is the number of turns of the limit cycle. The value '0' corresponds
   * to an equilibrium solution.
   */
  list<balSummaryEntry *> *summary;
  /** tells whether solutions and summary have been allocated */
  bool destroy_lists;

  boost::thread * logger_thread;
  boost::mutex list_mutex;
  boost::condition_variable q_empty;
  boost::condition_variable q_full;

  /** The number of threads that will be created to perform the integrations. */
  int nthreads;
};

#endif
