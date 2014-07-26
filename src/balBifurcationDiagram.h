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
 * \brief Definition of the class BifurcationDiagram.
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

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <list>
#include <boost/thread.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>
#include <boost/ref.hpp>

//#define DEBUG

namespace bal {

void ResetColours(int d);

/** The type of dynamical system analysis. */
typedef enum {
	/** Bifurcation analysis in parameter domain specified by BifurcationParameters object. */		
	PARAMS = 0,
	/** Analysis of dynamical system behaviour (parameters fixed) by varying initial conditions. */
	IC,
	/** Computation of the maximum Lyapunov exponent. */
	MLE
} diagram_mode;

/**
 * \class SummaryEntry
 * \brief Class used to store only the essential results of an integration (depending on choosen \ref diagram_mode and \ref integration_mode).
 * \ref data vector is organized as follows:
 * \f[
 * \begin{tabular}{|c|cccccccccc|}
 * \hline
 * diagram\_mode & & &	&	&	&	& & & &\\
 * \hline
 *	$IC$ & $p_0$ & $\dots$ & $p_m$ & $x_0^0$ & $\dots$ & $x_0^n$ & $x_{final}^0$ & $\dots$ & $x_{final}^n$ & $\#turns$ \\
 * \hline
 *  $PARAM$ & $p_0$	&	$\dots$	&	$p_m$	& $\#turns$ & & & & & & \\
 * \hline
 * $PARAM_{Lyap}$ & $p_0$	&	$\dots$	&	$p_m$	& $MLE$ & & & & & & \\
 * \hline
 * \end{tabular}
 * \f]
 *
 * \sa bal::BifurcationDiagram::SaveSummaryData bal::BifurcationDiagram::GetSummaryData
 */
class SummaryEntry : public Object {
 public:
  SummaryEntry(Solution *sol, diagram_mode mode = PARAMS);
  virtual ~SummaryEntry();
  int GetN() const;
  int GetID() const;
  double* GetData() const;

 private:
  double* data;
  int n;
  int id;
};

bool CompareSummaryEntry(SummaryEntry *entry1, SummaryEntry *entry2);

/**
 * \class BifurcationDiagram
 * \brief Class to calculate brute-force bifurcation diagrams.
 *
 * This class is used to compute brute-force bifurcation diagrams: it
 * manages 4 different kind of objects: a DynamicalSystem that defines
 * the RHS function of the system to integrate, a ODESolver that
 * performs the actual integration, a Logger to store data to file and a
 * BifurcationParameters that is mainly used to iterate over all the
 * combinations of parameters that the user is interested to simulate
 * (\link bifdiag.cpp \endlink).
 *
 * Otherwise the user can perform a basin of attraction analysis
 * utilising \ref SetInitialConditions method to define a set
 * of initial conditions over which the integrations will be performed (\link basin.cpp \endlink).
 *
 * \example basin.cpp
 * \example bifdiag.cpp 
 *
 * \sa DynamicalSystem BifurcationParameters ODESolver
 */

class BifurcationDiagram : public Object {
 public:
  BifurcationDiagram();
  BifurcationDiagram(const BifurcationDiagram& bifd);
  virtual ~BifurcationDiagram();

  std::string ToString() const;
  Object* Clone() const;

  /**
   * Sets the dynamical system to integrate. BifurcationDiagram
   * assumes that the dynamical system contains an instance of BifurcationParameters,
   * instead of the simpler Parameters.
   * \param sys A dynamical system: any instance of a class inherited
   * from DynamicalSystem.
   * \sa bal::DynamicalSystem
   */
  void SetDynamicalSystem(DynamicalSystem *sys);

  /**
   * Gets the dynamical system to integrate.
   * \return The dynamical system used in the computation of the
   * bifurcation diagram.
   */
  DynamicalSystem* GetDynamicalSystem() const;

  /**
   * Sets the logger used for saving integration data to file. A logger
   * is automatically instantiated when a BifurcationDiagram is
   * created. This method should be used only if the user wants to use a
   * different kind of logger (such as one, for example, that saves data
   * in a particular format). The default logger uses H5 files.
   * \param log An instance of one of the classes inherited by Logger.
   * \sa bal::Logger
   */
  void SetLogger(Logger *log);

  /**
   * \return The logger used for saving data to file.
   */
  Logger* GetLogger() const;

  /**
   * Sets the ODESolver used to integrate the dynamical system. An ODE
   * solver is automatically instantiated when a BifurcationDiagram is
   * created. This method should be called only if the user has developed
   * their own ODE solver.
   * \param sol An instance of an ODE solver.
   * \sa bal::ODESolver
   */
  void SetODESolver(ODESolver *sol);

  /**
   * This method returns a pointer to the ODE solver used to integrate
   * the system. It is useful to set parameters of the ODE solver.
   * \return The ODE solver used to integrate the system.
   */
  ODESolver* GetODESolver() const;

  /**
   * Sets the name of the file where data will be saved.
   * \param filename File where the bifurcation diagram will be saved.
   * \param compress Flag that tells wheter the file should be compressed.
   */
  void SetFilename(const char *filename, bool compress = false);

  /**
   * \return The name of the file where data is saved.
   */
  const char* GetFilename();

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
	 * \param restart If it is set to false each integration starts from the final state of
	 * the previous one. This may give an additional speed-up because the system
	 * could start already at steady state.
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

	/**
   * Saves the list of SummaryEntry resulting from ComputeDiagram in an ASCII file, sorted by parameters.
	 * \param filename Name of the file where results are saved.
   */
  bool SaveSummaryData(const char *filename) const;

	/**
   * Returns the results of ComputeDiagram as a matrix of SummaryEntry, sorted by parameters.
   */
  double** GetSummaryData(int *size = NULL) const;
	
  /** Defines what kind of analysis will be performed. */
  bool SetMode(diagram_mode _mode);
  int GetMode() const;
	
	/**
   * Initializes the set of initial conditions used for basins of attraction analysis.
	 * \param nx0 Number of initial conditions.
	 * \param x0  Matrix of initial conditions where the generic entry \f$x[i][j]\f$ refers to
	 *						the \f$j_{th}\f$ component of the \f$i_{th}\f$ initial condition.
   */
  void SetInitialConditions(int nx0, double **x0);

 private:

  void ComputeDiagramMultiThread();
  void IntegrateAndEnqueue(ODESolver *sol, int solutionId);
  double* BuildSummaryEntry(Solution *sol);

 private:

  /** The ODE solver used to integrate the system */
  boost::shared_ptr<ODESolver> solver;
  /** The dynamical system to integrate */
  boost::shared_ptr<DynamicalSystem> system;
  /** The parameters of the dynamical system */
  Parameters * parameters;

  /**
   * The object used to save data to a file: by default H5 file logging
   * is used, i.e.\ logger is an instance of the H5Logger class:
   * if the user wants to use another logger, they should provide it
   * by using the SetLogger member.
   */
  boost::shared_ptr<Logger> logger;

  /** The number of dimensions of the system. */
  int ndim;
  /** The total number of parameters of the system. */
  int npar;
  /** The number of parameters to be varied for the bifurcation diagram. */
  int n_var_par;
  /** The mode of computation of the diagram. */
  diagram_mode mode;

  int nX0;
  double **X0;

  /** Tells wheter each integration should be restarted from the initial
   * condition set by the user (true, default value) or from the final value of the
   * previous integration (false) */
  bool restart_from_x0;

  //// multithreading stuff ////
  /**
   * A list containing the results of the numerical integrations: when the list
   * is full, the threads that integrate stop and another thread saves data to file
   */
  std::list<Solution *> solutions;
  /**
   * A list containing the summary of the bifurcation diagram in
   * terms of number of turns of the solution. Each entry contains
   * (n+1) values, where n is the number of parameters of the system.
   * The first n values are the values of the parameters, the last value
   * is the number of turns of the limit cycle. The value '0' corresponds
   * to an equilibrium solution.
   */
  std::list<SummaryEntry *> summary;

  boost::thread logger_thread;
  boost::mutex list_mutex;
  boost::condition_variable q_empty;
  boost::condition_variable q_full;

  /** The number of threads that will be created to perform the integrations. */
  int nthreads;
};

} // namespace bal

#endif

