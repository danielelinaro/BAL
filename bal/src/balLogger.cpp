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

bool balLogger::SetFilename(const char *fname, bool open) {
	//printf("balLogger::SetFilename\n");
	if(IsFileOpen())
		CloseFile();
	strncpy(filename, fname, FILENAME_LENGTH);
	if(open) {
		opened = OpenFile();
		//printf("balLogger::SetFilename>> opened = %d\n", opened);
		return opened;
	}
	return true;
}

balH5Logger::~balH5Logger() {
	if(IsFileOpen())
		CloseFile();
}

bool balH5Logger::OpenFile() {
	//printf("balH5Logger::OpenFile %s\n", GetFilename());
	h5_fid = H5Fcreate(GetFilename(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(h5_fid == -1) {
		IsFileOpen(false);
		return false;
	}
	IsFileOpen(true);
	counter = 0;
	return true;
}

bool balH5Logger::CloseFile() {
	//printf("balH5Logger::CloseFile %s\n", GetFilename());
	if(!IsFileOpen())
		return false;
	int flag = H5Fclose(h5_fid);
	if(flag == 0)
		IsFileOpen(false);
	//printf("balH5Logger::CloseFile>> flag = %d\n", flag);
	return flag == 0;
}

bool balH5Logger::SaveBuffer(realtype * buffer, int rows) {
  hsize_t dims[2];
  herr_t status;

	if(buffer == NULL || rows <= 0 || GetNumberOfColumns() <= 0 || GetParameters() == NULL) {
		return false;
	}

	if(!IsFileOpen())
		OpenFile();

	// the name of the dataset
	counter++;
	sprintf(datasetname, "%06d", counter);
	
  // describe the size of the dataset
	dims[0] = rows;
	dims[1] = GetNumberOfColumns();

	status = H5LTmake_dataset_double(h5_fid, datasetname, 2, dims, buffer);
	status = H5LTset_attribute_double(h5_fid, datasetname, "parameters", GetParameters()->GetParameters(), GetParameters()->GetNumber());

  return true;
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
	
	hsize_t dims[2];
	herr_t status;
	balSolution * solution;
	
	/* balSolutionComparer is a struct defined in balSolution.h definig a method on operator() *
	 * to compare two balSolution pointers */
	sol_list->sort(balSolutionComparer());
	
	while (!sol_list->empty()) {
		
		solution = sol_list->front();
		sol_list->pop_front();
		
		if(solution->GetData() == NULL || solution->GetRows() <= 0 ||
		   solution->GetColumns() <= 0 || solution->GetParameters() == NULL) {
			return false;
		}
		
		// the name of the dataset
		counter++;
		sprintf(datasetname, "%06d", counter);
		
		// describe the size of the dataset
		dims[0] = (hsize_t) solution->GetRows();
		dims[1] = (hsize_t) solution->GetColumns();
		
		status = H5LTmake_dataset_double(h5_fid, datasetname, 2, dims, solution->GetData());
		status = H5LTset_attribute_double(h5_fid, datasetname,
										  "parameters", solution->GetParameters()->GetParameters(),
										  (size_t) solution->GetParameters()->GetNumber());
		
		solution->Destroy();
	}
	
	return true;
	
}


balASCIILogger::~balASCIILogger() {
	if(IsFileOpen())
		CloseFile();
}

bool balASCIILogger::OpenFile() {
	//printf("balASCIILogger::OpenFile\n");
	fp = fopen(GetFilename(),"w");
	if(fp == NULL) {
		IsFileOpen(false);
		return false;
	}
	IsFileOpen(true);
	return true;
}

bool balASCIILogger::CloseFile() {
	//printf("balASCIILogger::CloseFile\n");
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

