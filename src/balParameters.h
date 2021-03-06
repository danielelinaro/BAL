/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balParameters.h
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
 * \file balParameters.h
 * \brief Definition of the class balParameters
 */

#ifndef _BALPARAMETERS_
#define _BALPARAMETERS_

#include "balObject.h"
#include "balCommon.h"
#include <string>
#include <boost/shared_array.hpp>

namespace bal {

/**
 * \class Parameters 
 * \brief Base class for objects that contain parameters of
 * a dynamical system.
 * 
 * Parameters is the (general-purpose) base class to represent
 * parameters of a dynamical system. It simply contains an array that
 * represents the parameters of a system.  
 * This class can be inherited to contain more sophisticated definition of
 * parameters.
 *
 * \sa DynamicalSystem
 */
class Parameters : public Object {
public:
  Parameters (int numpars);
  Parameters (const Parameters& param);
  virtual ~Parameters ();
  virtual Object* Clone() const;
  virtual std::string ToString() const;

  int GetNumber () const;
  double* GetParameters () const;
  double At (int k) const;
  double& operator[] (int k);
  void operator=(const Parameters& param);

protected:
  int p;

private:
  boost::shared_array<double> pars;  
};


} // namespace bal

#endif

