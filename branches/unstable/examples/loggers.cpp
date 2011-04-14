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

#include "balObject.h"
#include "balParameters.h"
#include "balLogger.h"
using namespace bal;

// TEST H5Logger and ASCIILogger
int main(int argc, char *argv[]) {

  // data
  double buffer[10] = {0, 0, 0, 0, -2, 1, 1, 1, 1, -1}; 
  
  // parameters
  Parameters pars(4);
  pars[0] = 3.0;
  pars[1] = 5.0;
  pars[2] = 0.01;
  pars[3] = 4.0;
  
  // H5Logger
  H5Logger logger;
  logger.SetFilename("test.1.h5");
  logger.SetNumberOfColumns(5);
  logger.SetParameters(pars);
  logger.SaveBuffer(buffer, 2, 5683);
  logger.SetFilename("test.2.h5");
  logger.SaveBuffer(buffer, 2, 7583);

  return 0;
}

