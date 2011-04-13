
#include <cstdio>
#include <nvector/nvector_serial.h>
#include "balObject.h"
#include "balDynamicalSystem.h"
#include "balRossler.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"
#include "balBifurcationParameters.h"
using namespace std;

int main(int argc, char *argv[]) {
	
  int steps[3] = {1,4,5};
  realtype x0[3] = {0.5,0.5,0.5};
  BifurcationParameters * bp = BifurcationParameters::Create();
  bp->SetNumber(3);
  bp->SetIthParameter(0,0.25);
  bp->SetIthParameterLowerBound(1,0.2);
	bp->SetIthParameterUpperBound(1,0.3);
	bp->SetIthParameterLowerBound(2,9.3);
	bp->SetIthParameterUpperBound(2,9.6);
	//bp->SetIthParameter(2,4.3);
  bp->SetNumberOfSteps(steps);
  
	Rossler * ros = Rossler::Create();
  ros->SetParameters(bp);
  
	BifurcationDiagram * bifd = BifurcationDiagram::Create();
  bifd->SetDynamicalSystem(ros);
  bifd->GetODESolver()->SetIntegrationMode(balLYAP);
  bifd->GetODESolver()->SetTransientDuration(1e3);
	bifd->GetODESolver()->SetLyapunovTimeStep(1);
	bifd->GetODESolver()->SetTimeStep(1);
  bifd->GetODESolver()->SetFinalTime(1e4);
  bifd->GetODESolver()->SetX0(x0);
  
  bifd->SetNumberOfThreads(argc > 1 ? atoi(argv[1]) : 2);
	
  bifd->ComputeDiagram();
  bifd->SaveSummaryData("Rossler2.classified");
	
  bifd->Destroy();
  ros->Destroy();
  bp->Destroy();
  
  return 0;
}
