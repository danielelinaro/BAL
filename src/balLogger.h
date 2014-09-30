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


#include <sundials/sundials_types.h>

#include <string>
#include <list>
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
  Logger();
  Logger(const Logger& logger);
  virtual ~Logger();
  virtual void Open(const std::string& fname, bool compress = false) = 0;
  virtual void Close() = 0;

  std::string GetFilename() const;
  bool IsOpen() const;
  bool SaveSolution(Solution *solution);

  virtual Logger* Clone() const = 0;
  
 protected:
  virtual bool SaveBuffer(const Parameters* params,
			  const realtype *buffer, int rows, int columns,
			  int id) = 0;
  
 protected:
  /** The name of the file */
  std::string filename;
  /** Tells whether the logging file is open or not */
  bool file_is_open;
  /** Whether data compression should be enabled (default: no) */
  bool compressed;
};

void LoggerThread(Logger *logger, std::list<Solution*>& solutions,
		  boost::mutex& list_mutex,
		  boost::condition_variable& q_empty, boost::condition_variable& q_full);

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
  H5Logger();
  H5Logger(const std::string& fname, bool compress);
  H5Logger(const H5Logger& logger);
  virtual ~H5Logger();
  
  // inherited from Logger
  virtual void Open(const std::string& fname, bool compress = false);
  virtual void Close();
  
  virtual Logger* Clone() const;

 protected:
  // inherited from Logger
  virtual bool SaveBuffer(const Parameters* params,
			  const realtype *buffer, int rows, int columns,
			  int id);
 
 private:
  void EnableCompression();
  void DisableCompression();

 private:
  // the handle of the file
  hid_t h5_fid;
  // dataset creation property list
  hid_t dcpl;
  // chunk size
  hsize_t chunk[2];
};

} // namespace bal

#endif
