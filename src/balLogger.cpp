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

#include <iostream>
#include <sstream>
#include "balLogger.h"

namespace bal {

///// THREAD /////
void LoggerThread(Logger *logger, std::list<Solution*>& solutions,
		  boost::mutex& list_mutex,
		  boost::condition_variable& q_empty, boost::condition_variable& q_full) {
  if(!logger->IsOpen())
    return;
  
  Solution *s;
  while (true) {
    try {
      {
	boost::mutex::scoped_lock lock(list_mutex);
	while (solutions.size() < LIST_MAX_SIZE) {
	  q_full.wait(lock);  // INTERRUPTION POINT
	}
	solutions.sort(CompareSolutions);
	while (! solutions.empty()) {
	  s = solutions.front();
	  solutions.pop_front();
	  logger->SaveSolution(s);
	  delete s;
	}
      }
      /* notifies all threaded solvers the queue is now empty, giving them the control */
      q_empty.notify_all();
    }
    catch (boost::thread_interrupted&) {
      solutions.sort(CompareSolutions);
      while (! solutions.empty()) {
	s = solutions.front();
	solutions.pop_front();
	logger->SaveSolution(s);
	delete s;
      }
      break;
    }
  }
}

///// BALLOGGER /////

Logger::Logger() : file_is_open(false), compressed(false) {
#ifdef DEBUG
  std::cout << "Logger constructor.\n";
#endif
}

Logger::Logger(const Logger& logger)
  : file_is_open(logger.file_is_open), compressed(logger.compressed), filename(logger.filename) {
#ifdef DEBUG
  std::cout << "Logger copy constructor.\n";
#endif
}

Logger::~Logger() {
#ifdef DEBUG
  std::cout << "Logger destructor.\n";
#endif
}

std::string Logger::GetFilename() const {
  return filename;
}

bool Logger::IsOpen() const {
  return file_is_open;
}

bool Logger::SaveSolution(Solution *solution) {
  return SaveBuffer(solution->GetParameters(),
		    solution->GetData(), solution->GetRows(), solution->GetColumns(),
		    solution->GetID()); 
}

///// BALH5LOGGER /////

H5Logger::H5Logger() {
#ifdef DEBUG
  std::cout << "H5Logger constructor.\n";
#endif
  h5_fid = dcpl = -1;
  chunk[0] = chunk[1] = -1;
}

H5Logger::H5Logger(const std::string& fname, bool compress) : Logger(){
#ifdef DEBUG
  std::cout << "H5Logger copy constructor.\n";
#endif
  Open(fname,compress);
}

H5Logger::H5Logger(const H5Logger& logger) : Logger(logger) {
  h5_fid = logger.h5_fid;
  dcpl = logger.dcpl;
  chunk[0] = logger.chunk[0];
  chunk[1] = logger.chunk[1];
}

H5Logger::~H5Logger() {
#ifdef DEBUG
  std::cout << "H5Logger destructor.\n";
#endif
}

std::string H5Logger::ToString() const {
  return "H5Logger";
}

void H5Logger::DisableCompression() {
  compressed = false;
}

void H5Logger::EnableCompression() {
  // compression has already been enabled
  if(compressed)
    return;

  htri_t avail;
  herr_t status;
  unsigned int filter_info;

  std::cerr << "Checking whether GZIP compression is available...";
  // check if gzip compression is available
  avail = H5Zfilter_avail (H5Z_FILTER_DEFLATE);
  if (!avail) {
    std::cerr << "\nGZIP compression is not available on this system.\n";
    return;
  }
  std::cerr << " ok.\nGetting filter info...";
  status = H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info);
  if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ) {
    std::cerr << "\nUnable to get filter info: disabling compression.\n";
    return;
  }
    
  std::cerr << " ok.\nChecking whether the shuffle filter is available...";
  // check for availability of the shuffle filter.
  avail = H5Zfilter_avail(H5Z_FILTER_SHUFFLE);
  if (!avail) {
    std::cerr << "\nThe shuffle filter is not available on this system.\n";
    return;
  }
  std::cerr << " ok.\nGetting filter info...";
  status = H5Zget_filter_info (H5Z_FILTER_SHUFFLE, &filter_info);
  if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ) {
    std::cerr << "Unable to get filter info: disabling compression.\n";
    return;
  }
  std::cerr << " ok.\nCompression is enabled.\n";
  // enable gzip compression
  compressed = true;
}

void H5Logger::Open(const std::string& fname, bool compress) {
  if(file_is_open)
    Close();

  if(compress)
    EnableCompression();
  else
    DisableCompression();

  h5_fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(h5_fid == -1)
    return;
  if(compressed) {
    // Create the dataset creation property list and add the shuffle
    // filter and the gzip compression filter.
    // The order in which the filters are added here is significant -
    // we will see much greater results when the shuffle is applied
    // first.  The order in which the filters are added to the property
    // list is the order in which they will be invoked when writing
    // data.
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if(dcpl == -1)
      return;
    if(H5Pset_shuffle(dcpl) < 0)
      return;
    if(H5Pset_deflate(dcpl, 9) < 0)
      return;
  }
  file_is_open = true;
  filename = fname;
}

void H5Logger::Close() {
  if(file_is_open) {
    if(compressed)
      H5Pclose(dcpl);
    H5Fclose(h5_fid);
    file_is_open = false;
    filename = "";
  }
}

bool H5Logger::SaveBuffer(const Parameters *params,
			  const realtype *buffer, int rows, int columns,
			  int id) {
  if(!file_is_open || params == NULL || buffer == NULL || rows <= 0 || columns <= 0)
    return false;
  
  hsize_t dims[2];
  herr_t status;

  // describe the size of the dataset
  dims[0] = rows;
  dims[1] = columns;

  // create the name of the dataset: the letter C in the name means that the H5 file has been saved in C++
  char datasetname[14];
  sprintf(datasetname, "C%012d", id);

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
  else // do this if compression is disabled
    status = H5LTmake_dataset_double(h5_fid, datasetname, 2, dims, buffer);

  if(status < 0)
    return false;

  // If there are no parameters the function is not called
  int np = params->GetNumber();
  if(np > 0) {
    double *p = params->GetParameters();
    status = H5LTset_attribute_double(h5_fid, datasetname, "parameters", p, np);
  }

  return status >= 0;
}

} // namespace bal

