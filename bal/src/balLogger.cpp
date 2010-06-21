/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balLogger.cpp
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

#include "balLogger.h"

///// BALLOGGER /////

bool balLogger::SetFilename(const char *fname, bool open) {
  if(IsFileOpen())
    CloseFile();
  strncpy(filename, fname, FILENAME_LENGTH);
  if(open) {
    opened = OpenFile();
    return opened;
  }
  return true;
}

///// BALH5LOGGER /////

balH5Logger::balH5Logger() {
  h5_fid = -1;
  counter = -1;
  chunk[0] = chunk[1] = -1;

  htri_t avail;
  herr_t status;
  unsigned int filter_info;

  // check if gzip compression is available
  avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE);
  if (!avail) {
    printf("Disabling compression...\n");
    compressed = false;
    return;
  }
  status = H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info);
  if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ) {
    printf("Disabling compression...\n");
    compressed = false;
    return;
  }

  // check for availability of the shuffle filter.
  avail = H5Zfilter_avail(H5Z_FILTER_SHUFFLE);
  if (!avail) {
    printf("Disabling compression...\n");
    compressed = false;
    return;
  }
  status = H5Zget_filter_info (H5Z_FILTER_SHUFFLE, &filter_info);
  if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ) {
    printf("Disabling compression...\n");
    compressed = false;
    return;
  }

  // enable gzip compression
  compressed = true;
}

balH5Logger::~balH5Logger() {
  if(IsFileOpen())
    CloseFile();
}

bool balH5Logger::OpenFile() {
  h5_fid = H5Fcreate(GetFilename(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(h5_fid == -1) {
    IsFileOpen(false);
    return false;
  }
  IsFileOpen(true);
  counter = 0;
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

bool balH5Logger::CloseFile() {
  if(!IsFileOpen())
    return false;
  
  if(compressed) {
    H5Pclose (dcpl);
  }
  int flag = H5Fclose(h5_fid);
  if(flag == 0)
    IsFileOpen(false);
  return flag == 0;
}

bool balH5Logger::SaveBuffer(realtype * buffer, int rows) {
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

  counter++;
  // create the name of the dataset: the letter C in the name means that the H5 file has been saved in C++
  sprintf(datasetname, "C%06d", counter);
  
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

  status = H5LTset_attribute_double(h5_fid, datasetname, "parameters", GetParameters()->GetParameters(), GetParameters()->GetNumber());

  return status >= 0;
}

bool balH5Logger::SaveBufferThreaded (list <balSolution*> * sol_list,
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
	 * on condition variable q_full waiting queue and released when ComputeDiagramMultiThread() calls notify_one()
	 * on q_full when new solution is available 
	 */
	while (sol_list->size() < LIST_MAX_SIZE) {
	  
	  q_full->wait(lock);  // INTERRUPTION POINT
	  /* 
	   * Atomically call lock.unlock() and blocks the current thread. The thread will unblock when notified by a call
	   * to this->notify_one() or this->notify_all(). When the thread is unblocked (for whatever reason),
	   * the lock is reacquired by invoking lock.lock() before the call to wait returns. The lock is also reacquired by
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

bool balH5Logger::SortAndWriteSolutionList(list <balSolution *> * sol_list) {
  balSolution * solution;
  
  /* balSolutionComparer is a struct defined in balSolution.h definig a method on operator() *
   * to compare two balSolution pointers */
  sol_list->sort(balSolutionComparer());
  
  while (!sol_list->empty()) {
    solution = sol_list->front();
    sol_list->pop_front();
    SaveSolution(solution);
    solution->Destroy();
  }
  return true;
}

/*
///// BALASCIILOGGER /////

balASCIILogger::~balASCIILogger() {
  if(IsFileOpen())
    CloseFile();
}

bool balASCIILogger::OpenFile() {
  fp = fopen(GetFilename(),"w");
  if(fp == NULL) {
    IsFileOpen(false);
    return false;
  }
  IsFileOpen(true);
  return true;
}

bool balASCIILogger::CloseFile() {
  if(!IsFileOpen())
    return false;
  if(fclose(fp) == 0) {
    IsFileOpen(false);
    return true;
  }
  return false;
}

bool balASCIILogger::SaveBuffer(realtype * buffer, int rows) {
  int i, j;
  
  if(!IsFileOpen())
    OpenFile();
  
  if(fp == NULL || buffer == NULL || GetParameters() == NULL || rows <= 0) {
    return false;
  }
  
  // prints the parameters
  for(i=0; i<GetParameters()->GetNumber(); i++) {
    fprintf(fp, FORMAT, GetParameters()->At(i));
    if(i % GetNumberOfColumns() == 0 && i != 0) {
      fprintf(fp, NEWLINE);
    }
  }
  
  // ends the line
  while(i % GetNumberOfColumns() != 0) {
    fprintf(fp, FORMAT, 0.0);
    i++;
  }
  fprintf(fp, NEWLINE);
  
  // prints the buffer containing the data
  for(i=0; i<rows; i++) {
    for(j=0; j<GetNumberOfColumns(); j++) {
      fprintf(fp, FORMAT, buffer[i*GetNumberOfColumns() + j]);
    }
    fprintf(fp, "\n");
  }
  return true;
}
*/
