#include "balInterp2D.h"

double balBilinearInterp2D::interp(double x1p, double x2p) {
	int i, j;
	double yy, t, u;
	 
	i = x1terp->cor ? x1terp->hunt(x1p) : x1terp->locate(x1p);
	j = x2terp->cor ? x2terp->hunt(x2p) : x2terp->locate(x2p);

	t = (x1p-x1terp->xx[i])/(x1terp->xx[i+1]-x1terp->xx[i]); 
	u = (x2p-x2terp->xx[j])/(x2terp->xx[j+1]-x2terp->xx[j]); 
	yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j]
		+ (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
	return yy;
}

