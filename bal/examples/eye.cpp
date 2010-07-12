#include "balParameters.h"
#include "balEye.h"
#include "balODESolver.h"
#include "balSolution.h"

int main(int argc, char *argv[]) {

  double x0[2] = {1.0,-1.0};

  balParameters * pars = balParameters::Create();
  pars->SetNumber(0);

  balEye * eye = balEye::Create();

  if(! eye->SpecialOptions((void *) argv[1])) {
    printf("Unable to read vector field file. Aborting.\n");
    pars->Destroy();
    eye->Destroy();
    exit(1);
  }
  eye->SetParameters(pars);

  balODESolver * solver = balODESolver::Create();
  solver->SetDynamicalSystem(eye);
  solver->SetTransientDuration(0e7);
  solver->SetFinalTime(1e7);
  solver->SetTimeStep(500);
  solver->HaltAtEquilibrium(true);
  solver->SetEquilibriumTolerance(1e-2 * solver->GetEquilibriumTolerance());
  solver->SetIntegrationMode(balBOTH);
  solver->SetX0(x0);
  solver->Solve();
  
  FILE *fid = fopen("eye.dat","w");
  int r, c, i, j;
  double * buffer;
  balSolution * sol = solver->GetSolution();
  buffer = sol->GetData();
  sol->GetSize(&r,&c);
  for(i=0; i<r; i++) {
    for(j=0; j<c-1; j++)
      fprintf(fid, "%e ", buffer[i*c+j]);
    fprintf(fid,"%d\n", (int) buffer[i*c+j]);
  }
  fclose(fid);

  sol->Destroy();
  solver->Destroy();
  eye->Destroy();
  pars->Destroy();

  return 0;
}
