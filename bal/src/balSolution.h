/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balSolution.h
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

#ifndef _BALSOLUTION_
#define _BALSOLUTION_

#include "balObject.h"
#include "balCommon.h"
#include "balParameters.h"

class balSolution : public balObject {
 public:
  /** Returns the name of the class. */
  virtual const char * GetClassName() const { return "balSolution" ; }
  /** Creates a new balSolution. */
  static balSolution * Create() { return new balSolution; }
  /** Copies a balSolution */
  static balSolution * Copy(balSolution * solution) { return new balSolution(*solution); }
  /** Destroys a balSolution. */
  virtual void Destroy() { delete this; }
  /** Checks whether this object is of a particular type. */
  virtual bool IsA(const char * name) const { return (strcmp(name, this->GetClassName()) == 0); }
  
  int GetRows() const { return rows; }
  int GetColumns() const { return columns; }
  void GetSize(int * r, int * c) const {
    *r = rows;
    *c= columns;
  }
  void SetSize(int r, int c);
  
  balParameters * GetParameters() const { return parameters; }
  void SetParameters(balParameters * p);
  
  realtype * GetData() { return buffer; }	
  void SetData(int r, int c, realtype * data);
  
  int GetNumberOfTurns() { return nturns; }
  void SetNumberOfTurns(int _nturns) { nturns = _nturns; }
  
 protected:
  /* Protected destructor of the class. */
  virtual ~balSolution();
  balSolution();
  balSolution(const balSolution & solution);
  
 private:
  balParameters * parameters;
  realtype * buffer;
  int rows, columns;
  int nturns;
};

/** Used to compare two balSolution objects when sorting a list */
struct balSolutionComparer {
  bool operator() (balSolution * sol1, balSolution * sol2) {
    
    for (int i = 0; i < sol1->GetParameters()->GetNumber(); i++) {
      if (sol1->GetParameters()->At(i) < sol2->GetParameters()->At(i))
	return true;
    }
    
    return false;
  }
};



#endif


