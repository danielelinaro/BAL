#include "balInterpolator.h"

const char * balInterpolator::getClassName() const {
    return "balInterpolator";
 }

//virtual balInterpolator * Clone(balInterpolator *interp);

void balInterpolator::Destroy() {
   delete this;
 }

int balInterpolator::Init() {
  return 0;
 }

balInterpolator::balInterpolator(const balInterpolator & interp) {
  nnd = interp.nnd;
  nnf = interp.nnf;
}

balInterpolator::balInterpolator() {
  nnd = 0;
  nnf = 0;
}

balInterpolator::~balInterpolator() {}

int balInterpolator::GetDomainDimensions(){
  return nnd;
}

int balInterpolator::GetCodomainDimensions(){
  return nnf;
}
