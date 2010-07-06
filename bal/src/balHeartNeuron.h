#ifndef _BALHeartNeuron_
#define _BALHeartNeuron_

#include "balObject.h"
#include "balParameters.h"
#include "balDynamicalSystem.h"
#include <cvode/cvode.h>

class balHeartNeuron : public balDynamicalSystem {
 public:
  virtual const char * GetClassName () const { return "balHeartNeuron"; }
  static balHeartNeuron * Create () { return new balHeartNeuron; }
  virtual void Destroy () { delete this; }
  
  int RHS (realtype t, N_Vector x, N_Vector xdot, void * data);
#ifdef CVODE25
  int Jacobian (long int N, DenseMat J, realtype t, N_Vector x, N_Vector fy, 
		void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef CVODE26
  int Jacobian (int N, realtype t, N_Vector x, N_Vector fy, DlsMat J, 
		void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
  int Events (realtype t, N_Vector x, realtype * event, void * data);
  void EventsConstraints (realtype t, N_Vector x, int * constraints, void * data);
  
  bool HasJacobian() const { return false; }
  bool HasEvents() const { return true; }
  bool HasEventsConstraints() const { return true; }
  
  static realtype BoltzmannF(realtype a, realtype b, realtype V) {
    return 1./(1. + exp(a*(b + V)));
  }
  
  static realtype BoltzmannDFDV(realtype a, realtype b, realtype V) {
    realtype e = exp(a*(b+V));
    return - (a*e) / ((1.+e)*(1.+e));			
  }
  
 protected:
  balHeartNeuron();
  virtual ~balHeartNeuron();
  
 private:
  N_Vector xderiv;
  const realtype C,gK2,EK,ENa,gNa,E1,g1,tauNa;	// Constant parameters
  realtype A[3], B[3];
};


#ifdef __cplusplus
extern "C" {
#endif

balDynamicalSystem* balHeartNeuronFactory();
	
#ifdef __cplusplus
}
#endif

#endif
