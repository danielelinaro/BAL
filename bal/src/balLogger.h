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
 * \brief Definition of classes balLogger and balH5Logger
 *
 * balLogger and its inherited classes are used for saving data to files.
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

#define DONTCOMPRESS

using std::list;

/**
 * \class balLogger 
 * \brief Base class for saving integration data to file 
 * \sa balODESolver
 */
class balLogger : public balObject {
 public:
  virtual const char * GetClassName () const;
  static balLogger * Create();
  virtual void Destroy();
  
  virtual bool SetFilename(const char * fname, bool open = false);
  const char * GetFilename() const;
  void SetParameters(balParameters * p);
  balParameters * GetParameters() const;
  void SetNumberOfColumns(int c);
  int GetNumberOfColumns() const;
  bool IsFileOpen() const;
  
  virtual bool SaveBuffer(realtype * buffer, int rows, int id = 1);
  bool SaveSolution(balSolution * solution);
  bool SaveSolutionThreaded(list <balSolution *> * sol_list,
			    boost::mutex * list_mutex,
			    boost::condition_variable * q_empty,
			    boost::condition_variable * q_full);
  
 protected:
  balLogger();
  virtual ~balLogger();
  virtual bool OpenFile();
  virtual bool CloseFile();
  void SetFileIsOpen(bool open);

  virtual bool SortAndWriteSolutionList(list <balSolution *> * sol_list);
  
 private:
  /** Tells whether the logging file is open or not */
  bool opened;
  char filename[FILENAME_LENGTH];
  int cols;
  balParameters * params;
};

/**
 * \class balH5Logger 
 * \brief Class for saving integration data to H5 (compressed) files
 * \sa balLogger balODESolver
 */
class balH5Logger : public balLogger {
 public:
  virtual const char * GetClassName () const;
  static balH5Logger * Create();
  virtual void Destroy();
  
  virtual bool SaveBuffer(realtype * buffer, int rows, int id = 1);
  
 protected:
  balH5Logger();
  virtual ~balH5Logger();
  virtual bool OpenFile();
  virtual bool CloseFile();
  
 private:
  // the handle of the file
  hid_t h5_fid;
  // dataset creation property list
  hid_t dcpl;
  // chunk size
  hsize_t chunk[2];
  // whether data compression should be enabled (default: yes)
  bool compressed;
  // the name of the dataset in the H5 file
  char datasetname[DATASETNAME_LENGTH];
};


#endif

