/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balLogger.h
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

// .NAME balLogger - base class (balLogger) and inherited classes
// (balH5Logger and balASCIILogger) for saving integration data to files.
// 
// .SECTION Description
// balLogger and its inherited classes are used for saving data to files.
// More specifically, data is assumed to be a matrix of realtype values,
// where every row has the following structure:
//
// TIME X1 X2 ... Xi ... Xn LABEL
//
// Here TIME represents the instant of time at which data is recorded, Xi
// represent the value of the i-th state variable and LABEL is an
// additional field that gives information on the type of record, i.e.
// whether it is a normal integration step or it is an intersection with a
// Poincare' section. For more details on this aspect and for the possible
// values of the label column, see the documentation of the class
// balSolver.
//
// .SECTION See also
// balSolver
//

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
#define FORMAT "%15.6e "
#define NEWLINE "\n"


using std::list;

class balLogger : public balObject {
 public:
  virtual const char * GetClassName () const { return "balLogger"; }
  static balLogger * Create() { return new balLogger; }
  virtual void Destroy() { this->~balLogger(); }
  
  virtual bool SetFilename(const char * fname, bool open = false);
  const char * GetFilename() const { return filename; }
  void SetParameters(balParameters * p) { params = p; }
  balParameters * GetParameters() const { return params; }
  void SetNumberOfColumns(int c) { cols = c; }
  int GetNumberOfColumns() const { return cols; }
  bool IsFileOpen() const { return opened; }
  
  virtual bool SaveBuffer(realtype * buffer, int rows) { 
    if(!IsFileOpen()) 
      OpenFile(); 
    return false; 
  }
  
  virtual bool SaveBufferThreaded(list <balSolution *> * sol_list,
				  boost::mutex * list_mutex,
				  boost::condition_variable * q_empty,
				  boost::condition_variable * q_full) {
    if(!IsFileOpen()) 
      OpenFile();
    return false; 
  }
  
  virtual bool SaveSolution(balSolution * solution) { 
    SetNumberOfColumns(solution->GetColumns());
    SetParameters(solution->GetParameters());
    return SaveBuffer(solution->GetData(), solution->GetRows()); 
  }
  
  virtual bool SaveSolutionThreaded(list <balSolution *> * sol_list,
				    boost::mutex * list_mutex,
				    boost::condition_variable * q_empty,
				    boost::condition_variable * q_full) {
    return SaveBufferThreaded(sol_list,list_mutex,q_empty,q_full); 
  }
  
 protected:
 balLogger() : opened(false), cols(-1), params(NULL) {}
  virtual ~balLogger() {}
  virtual bool OpenFile() { return false; }
  virtual bool CloseFile() { return false; }
  void IsFileOpen(bool open) { opened = open; }
  
 private:
  /** Tells whether the logging file is open or not */
  bool opened;
  char filename[FILENAME_LENGTH];
  int cols;
  balParameters * params;
};

class balH5Logger : public balLogger {
 public:
  virtual const char * GetClassName () const { return "balH5Logger"; }
  static balH5Logger * Create() { return new balH5Logger; }
  virtual void Destroy() { this->~balH5Logger(); }
  
  virtual bool SaveBuffer(realtype * buffer, int rows);
  virtual bool SaveBufferThreaded(list <balSolution *> * sol_list,
				  boost::mutex * list_mutex,
				  boost::condition_variable * q_empty,
				  boost::condition_variable * q_full);
  
 protected:
 balH5Logger() : h5_fid(-1), counter(-1) {}
  virtual ~balH5Logger();
  virtual bool OpenFile();
  virtual bool CloseFile();
  
 private:
  
  bool SortAndWriteSolutionList(list <balSolution *> * sol_list);
  
  hid_t h5_fid;
  char datasetname[DATASETNAME_LENGTH];
  int counter;
};

class balASCIILogger : public balLogger {
 public:
  virtual const char * GetClassName () const { return "balASCIILogger"; }
  static balASCIILogger * Create() { return new balASCIILogger; }
  virtual void Destroy() { this->~balASCIILogger(); }
  
  virtual bool SaveBuffer(realtype * buffer, int rows);
  /*
    virtual bool SaveBufferThreaded(queue<balSolution *> *q,
    boost::mutex * list_mutex,
    boost::condition_variable *q_empty,
    boost::condition_variable *q_full,
    volatile bool *finished);
  */
  
 protected:
 balASCIILogger() : fp(NULL) {}
  virtual ~balASCIILogger();
  virtual bool OpenFile();
  virtual bool CloseFile();
  
 private:
  FILE * fp;
};


#endif

