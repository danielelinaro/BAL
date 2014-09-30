/**
 *  \file balInterpolator.cpp
 *  \brief Implementation of abstract class Interpolator. 
 */

#include "balInterpolator.h"

namespace bal {

int Interpolator::Init() {
  return 0;
 }

Interpolator::Interpolator(const Interpolator & interp) {
  nnd = interp.nnd;
  nnf = interp.nnf;
}

Interpolator::Interpolator() {
  nnd = 0;
  nnf = 0;
}

Interpolator::~Interpolator() {}

int Interpolator::GetDomainDimensions(){
  return nnd;
}

int Interpolator::GetCodomainDimensions(){
  return nnf;
}

} // namespace bal

