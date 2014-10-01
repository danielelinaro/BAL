/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    loggers.cpp
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

#include <iostream>
#include "balObject.h"
#include "balParameters.h"
#include "balSolution.h"
#include "balLogger.h"
using namespace bal;

// TEST H5Logger and ASCIILogger
int main(int argc, char *argv[]) {

  // data
  const int rows = 2;
  const int columns = 5;
  double buffer[rows*columns] = {0, 0, 0, 0, -2, 1, 1, 1, 1, -1}; 
  
  // Parameters
  Parameters pars(4);
  pars[0] = 3.0;
  pars[1] = 5.0;
  pars[2] = 0.01;
  pars[3] = 4.0;
  
  // Solution
  Solution sol(rows, columns, buffer);
  sol.SetParameters(&pars);

  int copies = 250;
  // H5Logger
  std::string filename;
  H5Logger logger;
  filename = "test.1.h5";
  logger.Open(filename,true);
  if(logger.IsOpen())
    std::cout << "Opened " << filename << " with compression enabled.\n";
  else {
    std::cout << "Error opening " << filename << ". Aborting.\n";
    return 1;
  }
  for(int i=1; i<=copies; i++) {
    sol.SetID(i);
    logger.SaveSolution(&sol);
  }
  std::cout << "Saved first file.\n";
  logger.Close();
  filename = "test.2.h5";
  logger.Open(filename,false);
  if(logger.IsOpen())
    std::cout << "Opened " << filename << " with compression disabled.\n";
  else {
    std::cout << "Error opening " << filename << ". Aborting.\n";
    return 1;
  }
  for(int i=1; i<=copies; i++) {
    sol.SetID(i);
    logger.SaveSolution(&sol);
  }
  std::cout << "Saved second file.\n";
  logger.Close();

  return 0;
}

