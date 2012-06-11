/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balLogger.h
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
 * \file balLogger.h
 * \brief Definition of classes Logger and H5Logger.\
 * Logger and its inherited classes are used for saving data to files.
 */

#ifndef _BALLOGGER_
#define _BALLOGGER_

#include <cstring>
#include <cstdio>
#include <sundials/sundials_types.h>

#include <list>
#include <algorithm>
#include <boost/thread.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>
#include <boost/ref.hpp>

#include "hdf5.h"
#include "H5LTpublic.h"
#include "balObject.h"
#include "balParameters.h"
#include "balSolution.h"
#include "balCommon.h"

#define FILENAME_LENGTH (200)
#define DATASETNAME_LENGTH (10)

using std::list;

namespace bal {

/**
 * \class Logger 
 * \brief Base class for saving integration data to file.
 * More specifically, each data buffer is assumed to be a matrix of realtype values,
 * where every row has the following structure:
 \f[
 \begin{tabular}{|c|ccccc|c|}
 \hline
 time	&	$x^1$	&	$\dots$	&	$x^i$	&	$\dots$	&	$x^n$	&	label	\\
 \hline
 $t_{start}$	&	$x_0^1$	&	$\dots$	&	$x_0^i$	&	$\dots$	&	$x_0^n$	&	$-2$	\\
 $t_{tran}$	&	$x_1^1$	&	$\dots$	&	$x_1^i$	&	$\dots$	&	$x_1^n$	&	$-1$	\\
 $t_j$	&	$x_j^1$	&	$\dots$	&	$x_j^i$	&	$\dots$	&	$x_j^n$	&	$0$		\\
 $\vdots$	&	$\vdots$	&	$\vdots$	&	$\vdots$	&	$\vdots$	&	$\vdots$	&	$\vdots$	\\
 $t_k$	&	$x_k^1$	&	$\dots$	&	$x_k^i$	&	$\dots$	&	$x_k^n$	&	$l$		\\
 \hline
 \end{tabular}
 \f]
 * where \a time represents the instant of time at which data is recorded, \f$x_i\f$
 * represent the value of the i-th state variable and \a label contains information about each entry: \a -2  for the initial condition, \a -1  at the
 * end of the transient period, \a 0  for a standard entry, every \a l  integer grater than zero indicates
 * a zero value of the \a l-th  event function at that point. Other possible values are: \a -3  at an equilibrium point 
 * and \a -10  if an error occurred during integration.
 * 
 * Each derived class can define how data are organized into the file and which is the data storage format used,
 * according to user and software requirements. 
 * 
 * \sa bal::state_label ODESolver Solution
 *
 * \example loggers.cpp
 * 
 */
	
class Logger : public Object {
 public:

  virtual const char * GetClassName () const;
  virtual void SetFilename(const char * fname, bool compress = false);
  const char * GetFilename() const;
  void SetParameters(Parameters * p);
  Parameters * GetParameters() const;
  void SetNumberOfColumns(int c);
  int GetNumberOfColumns() const;
  bool IsFileOpen() const;
  
	/** Base method to save an entry to file.\ Must be overriden by each derived class. */ 
  virtual bool SaveBuffer(realtype * buffer, int rows, int id = 1);
	/** It saves a Solution entry using \ref SaveBuffer method. */
  bool SaveSolution(Solution * solution);
	/** In multithread mode, try to get the mutex for a safe reading of the shared
	 *	solutions list (filled by integration threads), then calls \ref SortAndWriteSolutionList
	 *	method and release mutex.
	 */
  bool SaveSolutionThreaded(list <Solution *> * sol_list,  /// pubblica?
			    boost::mutex * list_mutex,
			    boost::condition_variable * q_empty,
			    boost::condition_variable * q_full);
  
 protected:
  Logger();
  virtual ~Logger();
  virtual bool OpenFile() = 0;
  virtual bool CloseFile() = 0;
  void SetFileIsOpen(bool open);
	/** Sorts and writes to file the first element popped from the shared solutions list (actually a queue). */
  virtual bool SortAndWriteSolutionList(list <Solution *> * sol_list);
  
 private:
  /** Tells whether the logging file is open or not */
  bool opened;
  char filename[FILENAME_LENGTH];
  int cols;
  Parameters * params;
};

/**
 * \class H5Logger 
 * \brief Class for saving integration data to H5 (compressed) file.
 *
 *	Please note that are available two MATLAB funtions to read and classify H5 datasets at [BAL_folder]/matlab/.
 *
 *	More informations on H5 format and HDF5 library can be found <a href="http://www.hdfgroup.org/HDF5/">here</a>.
 * \sa Logger ODESolver
 */
class H5Logger : public Logger {
 public:
  virtual const char * GetClassName () const;
  static H5Logger * Create();
  virtual void Destroy();
	/** Sets H5 dataset filename and allows to enable data compression. */
  virtual void SetFilename(const char * fname, bool compress = false);
	/** \brief Performs data storage using HDF5 library functions. */
  virtual bool SaveBuffer(realtype * buffer, int rows, int id = 1);
  
 protected:
  H5Logger();
  virtual ~H5Logger();
  virtual bool OpenFile();
  virtual bool CloseFile();
  
 private:
  // the handle of the file
  hid_t h5_fid;
  // dataset creation property list
  hid_t dcpl;
  // chunk size
  hsize_t chunk[2];
  /** Indicates if data compression is enabled (default: no). */
  bool compressed;
  /** The name of the dataset in the H5 file. */
  char datasetname[DATASETNAME_LENGTH];
};

} // namespace bal

#endif

