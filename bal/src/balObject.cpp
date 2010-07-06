/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balObject.cpp
 *
 *   Copyright (C) 2010 Daniele Linaro
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

balObject::balObject() {
}

balObject::~balObject() {
}

balObject* balObject::Create() {
  return new balObject;
}

void balObject::Destroy() {
  delete this;
}

const char* balObject::GetClassName() const {
  return "balObject";
}

bool balObject::IsA(const char * name) const {
  return (strcmp(name, this->GetClassName()) == 0);
}
