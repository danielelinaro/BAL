/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balObject.h
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
 * \file balObject.h 
 * \brief Definition of the class balObject
 */

#ifndef _BALOBJECT_
#define _BALOBJECT_

#include <string>

namespace bal {

/**
 *  \class balObject
 *  \brief Base class for all BAL objects.
 *  
 *  Object is the base class for all objects in the Bifurcation Analysis
 *  Library. Every object in the library should be a subclass of balObject.
 */
class Object {

 public:
  /** Constructor of the class. */
  Object();
  /** Destructor of the class. */
  virtual ~Object();
  
  virtual std::string ToString() const = 0;

};

} // namespace bal

#endif

