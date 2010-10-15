/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    obj.cpp
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
using namespace std;

// TEST balObject
int main(int argc, char *argv[]) {
  bool flag;
  const char name[] = "balObject";
  balObject *obj = balObject::Create();

  if(argc > 1)
    flag = obj->IsA(argv[1]);
  else
    flag = obj->IsA(name);

  cout << obj->GetClassName() << endl;
  cout << "I am " << (flag ? "" : "not ") << "a " << (argc>1 ? argv[1] : name) << "." << endl;
  obj->Destroy();

  return 0;
}

