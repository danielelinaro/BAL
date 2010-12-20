#include "balParameters.h"
#include "balInterpSystem.h"
#include "balODESolver.h"
#include "balInterp3D.h"

#define DIM 3


int main(int argc, char *argv[]) {

  int i,j;
  int dim;
  double *xi1, *xi2, *xi3, **yi; 
  int nx1, nx2, nx3, nx12, nx123;
  double x0[3] = {0.1,-0.1, 0.1}; 


  // Instantiate eye
  balInterpSystem * interpsystem = balInterpSystem::Create();
  
  // Read vector field
  FILE *fdat = fopen(argv[1],"r");
  if (fdat == NULL) {
    fprintf(stderr,"Unable to read vector field. Aborting.\n");
    exit(0);
  }
  fscanf(fdat,"%d",&dim);
  fscanf(fdat,"%d %d %d",&nx1, &nx2, &nx3);
  nx12 = nx1*nx2;
  nx123 = nx1*nx2*nx3;
  xi1 = new double[nx1];
  xi2 = new double[nx2];
  xi3 = new double[nx3];
  yi = new double * [DIM];
  for (i=0; i<DIM; i++)
    yi[i] = new double[nx123];
  for (i=0; i<nx123; i++) {
    fscanf(fdat,"%le %le %le %le %le %le",&xi1[i%nx1],&xi2[(i%nx12)/nx1],&xi3[i/nx12],&yi[0][i],&yi[1][i],&yi[2][i]);
  }
  fclose(fdat);
    
  balSplineInterp3D * interp = balSplineInterp3D::Create();
  interp->SetInterpolationPoints(xi1,xi2,xi3,yi,nx1,nx2,nx3,DIM);
  interp->SetWindow(10);
  interp->Init();

  interpsystem->SetInterpolator(interp);
  
  bool backward = true;
  interpsystem->SpecialOptions(&backward);

  balODESolver * solver = balODESolver::Create();
  solver->SetDynamicalSystem(interpsystem);
  solver->SetTransientDuration(0e6);
  solver->SetFinalTime(1e6);
  solver->SetTimeStep(10);
  solver->HaltAtEquilibrium(true);
  solver->SetEquilibriumTolerance(1e-2 * solver->GetEquilibriumTolerance());
  solver->SetIntegrationMode(balTRAJ);
  solver->SetX0(x0);
  solver->Solve();

  FILE *fid = fopen("interpsystem.dat","w");
  int r, c;
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

  
 // Free memory
  solver->Destroy();
  interp->Destroy();
  interpsystem->Destroy();
  for (i=0; i<DIM; i++)
    delete [] yi[i];
  delete [] yi;
  delete [] xi1;
  delete [] xi2;
  delete [] xi3;

  return 0;
}
