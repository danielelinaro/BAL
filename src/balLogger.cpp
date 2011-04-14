/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balLogger.cpp
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
 * \file balLogger.cpp
 * \brief Implementation of classes Logger and H5Logger
 */

#include "balLogger.h"

namespace bal {

///// BALLOGGER /////

Logger::Logger() : opened(false), cols(-1), params(NULL) {}

Logger::~Logger() {}

void Logger::SetFileIsOpen(bool open) {
  opened = open;
}

const char * Logger::GetClassName () const {
  return "Logger";
}

const char * Logger::GetFilename() const {
  return filename;
}

void Logger::SetParameters(Parameters * p) {
  params = p;
}

Parameters * Logger::GetParameters() const {
  return params;
}

void Logger::SetNumberOfColumns(int c) {
  cols = c;
}

int Logger::GetNumberOfColumns() const {
  return cols;
}

bool Logger::IsFileOpen() const {
  return opened;
}

bool Logger::SaveBuffer(realtype * buffer, int rows, int id) { 
  if(!IsFileOpen()) 
    OpenFile(); 
  return false; 
}

bool Logger::SaveSolution(Solution * solution) { 
  SetNumberOfColumns(solution->GetColumns());
  SetParameters(solution->GetParameters());
  return SaveBuffer(solution->GetData(), solution->GetRows(), solution->GetID()); 
}

bool Logger::SaveSolutionThreaded(list <Solution *> * sol_list,
				     boost::mutex * list_mutex,
				     boost::condition_variable * q_empty,
				     boost::condition_variable * q_full) {
  if(!IsFileOpen())
    OpenFile();
  
  // interrupt enabled by default
  
  /* data writing routine. */
  while (true) {
    try {
      {
	boost::mutex::scoped_lock lock(*list_mutex);
	/* 
	 * if the data queue is empty the thread is set on wait 
	 * on condition variable q_full waiting queue and released
	 * when ComputeDiagramMultiThread() calls notify_one()
	 * on q_full when new solution is available 
	 */
	while (sol_list->size() < LIST_MAX_SIZE) {
	  
	  q_full->wait(lock);  // INTERRUPTION POINT
	  /* 
	   * Atomically call lock.unlock() and blocks the current thread.
	   * The thread will unblock when notified by a call
	   * to this->notify_one() or this->notify_all(). 
	   * When the thread is unblocked (for whatever reason),
	   * the lock is reacquired by invoking lock.lock() before the call to 
	   * wait returns. The lock is also reacquired by
	   * invoking lock.lock() if the function exits with an exception.
	   */
	}
	
	SortAndWriteSolutionList(sol_list);
	
      }
      /* 
       * this scope has been introduced due to scoped_lock class implementation: 
       * mutex list_mutex is locked in lock initialization and unlocked while exiting the scope 
       */
      /* notifies all threaded solvers the queue is now empty, giving them the control */
      q_empty->notify_all();
      
    }
    catch (boost::thread_interrupted&) {
      SortAndWriteSolutionList(sol_list);
      break;
    }
  }
  return true;
}

bool Logger::SortAndWriteSolutionList(list <Solution *> * sol_list) {
  Solution * solution;
  
  /* SolutionComparer is a struct defined in balSolution.h defining a method on operator() *
   * to compare two Solution pointers */
  //sol_list->sort(SolutionComparer());
  sol_list->sort(CompareBalSolutions);
  
  while (!sol_list->empty()) {
    solution = sol_list->front();
    sol_list->pop_front();
    SaveSolution(solution);
    solution->Destroy();
  }
  return true;
}

void Logger::SetFilename(const char *fname, bool compress) {
  strncpy(filename, fname, FILENAME_LENGTH);
}

///// BALH5LOGGER /////

const char * H5Logger::GetClassName () const {
  return "H5Logger";
}

H5Logger * H5Logger::Create() {
  return new H5Logger;
}

void H5Logger::Destroy() {
  delete this;
}

H5Logger::H5Logger() {
  h5_fid = -1;
  chunk[0] = chunk[1] = -1;
}

H5Logger::~H5Logger() {
  if(IsFileOpen())
    CloseFile();
}

void H5Logger::SetFilename(const char * fname, bool compress) {
  Logger::SetFilename(fname,compress);

  compressed = false;
  if(compress) {
    htri_t avail;
    herr_t status;
    unsigned int filter_info;

    printf("Checking whether GZIP compression is available...");
    // check if gzip compression is available
    avail = H5Zfilter_avail (H5Z_FILTER_DEFLATE);
    if (!avail) {
      printf("\nGZIP compression is not available on this system.\n");
      return;
    }
    printf(" ok.\nGetting filter info...");
    status = H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info);
    if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ) {
      printf("\nUnable to get filter info: disabling compression.\n");
      return;
    }
    
    printf(" ok.\nChecking whether the shuffle filter is available...");
    // check for availability of the shuffle filter.
    avail = H5Zfilter_avail(H5Z_FILTER_SHUFFLE);
    if (!avail) {
      printf("\nThe shuffle filter is not available on this system.\n");
      return;
    }
    printf(" ok.\nGetting filter info...");
    status = H5Zget_filter_info (H5Z_FILTER_SHUFFLE, &filter_info);
    if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ) {
      printf("Unable to get filter info: disabling compression.\n");
      return;
    }
    printf(" ok.\nCompression is enabled.\n");
    // enable gzip compression
    compressed = true;
  }
  else {
    printf("Compression is disabled.\n");
  }
}

bool H5Logger::OpenFile() {
  if(IsFileOpen()) CloseFile();

  h5_fid = H5Fcreate(GetFilename(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(h5_fid == -1) {
    SetFileIsOpen(false);
    return false;
  }
  SetFileIsOpen(true);
  if(compressed) {
    /*
     * Create the dataset creation property list and add the shuffle
     * filter and the gzip compression filter.
     * The order in which the filters are added here is significant -
     * we will see much greater results when the shuffle is applied
     * first.  The order in which the filters are added to the property
     * list is the order in which they will be invoked when writing
     * data.
     */
    herr_t status;
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_shuffle (dcpl);
    status = H5Pset_deflate (dcpl, 9);
  }
  return true;
}

bool H5Logger::CloseFile() {
  if(!IsFileOpen())
    return false;
  
  if(compressed) {
    H5Pclose (dcpl);
  }
  int flag = H5Fclose(h5_fid);
  if(flag == 0)
    SetFileIsOpen(false);
  return flag == 0;
}

bool H5Logger::SaveBuffer(realtype * buffer, int rows, int id) {
  if(buffer == NULL || rows <= 0 || GetNumberOfColumns() <= 0 || GetParameters() == NULL) {
    return false;
  }
  
  if(!IsFileOpen())
    OpenFile();
  
  hsize_t dims[2];
  herr_t status;

  // describe the size of the dataset
  dims[0] = rows;
  dims[1] = GetNumberOfColumns();

  // create the name of the dataset: the letter C in the name means that the H5 file has been saved in C++
  sprintf(datasetname, "C%06d", id);
  
  if(compressed) { // do this if compression is enabled
    // set chunk size
    if(dims[0] > 10)
      chunk[0] = floor(dims[0]/10);
    else
      chunk[0] = dims[0];
    chunk[1] = dims[1];

    status = H5Pset_chunk (dcpl, 2, chunk);
    if(status < 0)
      return false;

    // create dataspace: setting maximum size to NULL sets the maximum
    // size to be the current size.
    hid_t space = H5Screate_simple (2, dims, NULL);
    if(space < 0)
      return false;

    // create the dataset.
    hid_t dset = H5Dcreate (h5_fid, datasetname, H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if(dset < 0) {
      H5Sclose(space);
      return false;
    }

    // write the data to the dataset.
    status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    
    // close and release resources.
    H5Dclose (dset);
    H5Sclose (space);
  }
  else { // do this if compression is disabled
    status = H5LTmake_dataset_double(h5_fid, datasetname, 2, dims, buffer);
  }
  if(status < 0)
    return false;

  // If there are no parameters the function is not called
  if (GetParameters()->GetNumber() > 0)
    status = H5LTset_attribute_double(h5_fid, datasetname, "parameters", GetParameters()->GetParameters(), GetParameters()->GetNumber());

  return status >= 0;
}

} // namespace bal

