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

/** 
 * \file balSolution.cpp
 * \brief Implementation of the class balSolution
 */

#include "balSolution.h"

balSolution::balSolution() {
  rows = columns = 0;
  buffer = NULL;
  parameters = NULL;
  ID = 0;
}

balSolution::balSolution(const balSolution& solution) {
  balSolution();
  solution.GetSize(&rows,&columns);
  SetSize(rows,columns);
  memcpy(buffer,solution.buffer,rows*columns*sizeof(realtype));
  nturns = solution.nturns;
  ID = solution.ID;
}

balSolution::~balSolution() {
  if(buffer != NULL) delete buffer;
  if(parameters != NULL) parameters->Destroy();
}

const char * balSolution::GetClassName() const {
  return "balSolution";
}

balSolution * balSolution::Create() {
  return new balSolution;
}

balSolution * balSolution::Copy(balSolution * solution) {
  return new balSolution(*solution);
}

void balSolution::Destroy() {
  delete this;
}

int balSolution::GetRows() const {
  return rows;
}

int balSolution::GetColumns() const {
  return columns;
}

void balSolution::SetSize(int r, int c) {
  if(buffer != NULL) delete buffer;
  rows = r;
  columns = c;
  buffer = new realtype[rows*columns];
}

void balSolution::GetSize(int * r, int * c) const {
  *r = rows;
  *c= columns;
}

void balSolution::SetData(int r, int c, realtype * data) {
  SetSize(r,c);
  memcpy(buffer,data,rows*columns*sizeof(realtype));
}

realtype * balSolution::GetData() const {
  return buffer;
}	

void balSolution::SetParameters(balParameters * p) {
  if(parameters != NULL) parameters->Destroy();
  parameters = balParameters::Copy(p);
}

balParameters * balSolution::GetParameters() const {
  return parameters;
}

void balSolution::SetNumberOfTurns(int _nturns) {
  nturns = _nturns;
}
  
int balSolution::GetNumberOfTurns() const {
  return nturns;
}

void balSolution::SetID(int id) {
  ID = id;
}

int balSolution::GetID() const {
  return ID;
}

bool CompareBalSolutions(balSolution *sol1, balSolution *sol2) {
  return sol1->GetID() < sol2->GetID();
}
