/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balParameters.h
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

// .NAME balParameters - base class for objects that contain parameters of
// a dynamical system.
// 
// .SECTION Description
// balParameters is the (general-purpose) base class to represent
// parameters of a dynamical system. It simply contains an array that
// represents the parameters of a system.  
// This class can be inherited to contain more sophisticated definition of
// parameters.
//
// .SECTION See also
// balDynamicalSystem
//

#ifndef _BALPARAMETERS_
#define _BALPARAMETERS_

#include "balObject.h"
#include <fstream>
using namespace std;

class balParameters : public balObject {
 public:
  virtual const char * GetClassName () const { return "balParameters"; }
  static balParameters * Create () { return new balParameters; }
  static balParameters * Copy (balParameters * params) { return new balParameters(*params); }
  virtual void Destroy () { this->~balParameters(); }
  virtual void SetNumber (int);
  int GetNumber () const;
  
  double & At (int k);
  double * GetParameters() const;
  
  void CopyValues(balParameters* _par);
  
  friend ostream & operator<< (ostream & out, const balParameters & bp) {
    out << "(";
    for(int i=0; i<bp.p-1; i++)
      out << bp.pars[i] << ",";
    out << bp.pars[bp.p-1] << ")";
    return out;
  }
  
 protected:
  balParameters ();
  balParameters (const balParameters & param);
  virtual ~balParameters () { if (pars != NULL) delete pars; }
  
 private:
  int p;
  double * pars;
  bool _dealloc;
};

#endif

