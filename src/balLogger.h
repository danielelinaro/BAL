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
 * \brief Definition of classes Logger and H5Logger
 *
 * Logger and its inherited classes are used for saving data to files.
 * More specifically, data is assumed to be a matrix of realtype values,
 * where every row has the following structure:
 *
 * TIME X1 X2 ... Xi ... Xn LABEL
 *
 * Here TIME represents the instant of time at which data is recorded, Xi
 * represent the value of the i-th state variable and LABEL is an
 * additional field that gives information on the type of record, i.e.
 * whether it is a normal integration step or it is an intersection with a
 * Poincare' section. For more details on this aspect and for the possible
 * values of the label column, see the documentation of the class
 * balSolver.
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
 * \brief Base class for saving integration data to file 
 * \sa ODESolver
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
 * \brief Class for saving integration data to H5 (compressed) files
 * \sa Logger ODESolver
 */

class H5Logger : public Logger {
 public:
  H5Logger();
  H5Logger(const std::string& fname, bool compress);
  H5Logger(const H5Logger& logger);
  virtual ~H5Logger();
  
  // inherited from Object
  Object* Clone() const;
  std::string ToString() const;

  // inherited from Logger
  virtual void Open(const std::string& fname, bool compress = false);
  virtual void Close();
  
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
