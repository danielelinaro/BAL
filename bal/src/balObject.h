/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balObject.h
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

/** 
 * \file balObject.h 
 * \brief Definition of the class balObject
 */

#ifndef _BALOBJECT_
#define _BALOBJECT_

#include <cstdio>
#include <cstring>

/**
 *  \class balObject
 *  \brief Base class for all BAL objects.
 *  
 *  balObject is the base class for all objects in the Bifurcation Analysis
 *  Library. Every object in the library should be a subclass of balObject.
 *  Constructor and destructor of the subclasses of balObject
 *  should be protected, so that only Create() and Destroy() actually
 *  call them.
 *  Note: Objects of subclasses of balObjects should always be
 *  created with the Create() method and deleted with the Destroy()
 *  method. They cannot be allocated off the stack (i.e., automatic
 *  objects) because the constructor is a protected method.
 */
class balObject {

 public:

  /** Returns the name of the class. */
  virtual const char * GetClassName() const;

  /** Creates a new balObject. */
  static balObject * Create();

  /** Destroys a balObject. */
  virtual void Destroy();

  /** Checks whether this object is of a particular type. */
  virtual bool IsA(const char * name) const;

 protected:

  /** Protected constructor of the class. */
  balObject();

  /** Protected destructor of the class. */
  virtual ~balObject();
};

#endif

