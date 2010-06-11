#!/usr/bin/env python

from bal import *
from pylab import figure, plot, xlabel, ylabel, title, show, axis

### START FUNCTIONS ###

def gsr(x):
    from numpy import zeros, dot, sqrt
    from numpy.linalg import norm
    N = len(x)
    n = int(sqrt(N))
    znorm = zeros(n)
    znorm[0] = norm(x[0:n])
    xnorm = zeros(N)
    xnorm[0:n] = x[0:n] / znorm[0]

    for i in range(1,n):
        tmp = x[n*i:n*(i+1)].copy()
        for j in range(i):
            d = dot(x[n*i:n*(i+1)],xnorm[n*j:n*(j+1)])
            tmp = tmp - dot(x[n*i:n*(i+1)],xnorm[n*j:n*(j+1)])*xnorm[n*j:n*(j+1)]
        znorm[i] = norm(tmp)
        xnorm[n*i:n*(i+1)] = tmp / znorm[i]
    return xnorm,znorm

def lyapunov(dynsys,pars,tstart,tend,tstep,ic):
    from numpy import zeros, array, sqrt, log, reshape, eye
    from numpy.linalg import qr, norm
    n = len(ic)
    N = n*(n+1)
    x = zeros(N)
    cum = zeros(n)
    lp = zeros(n)
    x[0:n] = ic
    x[n:] = list(reshape(eye(n),[n**2,1]).flatten())
    dynsys.extended = 1
    solver = ODESolver(dynsys,pars)
    solver.x0 = list(x)
    solver.dt = tstep/5
    solver.ttran = 0.0
    solver.tstop = tstep
    solver.mode = 'trajectory'
    t = 0
    while t<tend:
        solver.run()
        s = solver.solution()
        t = t+tstep
        x = array(s.data['x'][-N:])

        # construct a new orthonormal basis by Gram-Schmidt #
        xnorm,znorm = gsr(x[n:])
        x[n:] = xnorm
        # update running vector magnitudes
        for i in range(n):
            cum[i] = cum[i] + log(znorm[i])/log(2.0)
        # normalize the exponents
        for i in range(n):
            lp[i] = cum[i]/(t-tstart)
        #print 'lp = ', lp
        solver.x0 = list(x)
    return lp

#### END FUNCTIONS ####

from numpy import array, reshape
from numpy.linalg import qr
from numpy.random import random

# the system
sys = DynamicalSystem()
#sys.create('balLorenz')
sys.create('balHindmarshRose')

# the parameters
par = Parameters(sys.npar)

### LORENZ ###
# chaos
#par.setpars([10.,28,8./3.])
# limit cycle
#par.setpars([10.,99.96,8./3.])
# stable equilibrium
#par.setpars([10.,10,8./3.])

### HR ###
# chaos
#par.setpars([2.96,3,0.01,4])
# limit cycle
#par.setpars([4,5,0.01,4])
# equilibrium
#par.setpars([4,1,0.01,4])

solver = ODESolver(sys,par)
solver.dt = 0.01
solver.ttran = 1000
solver.tstop = 2000
solver.mode = 'trajectory + events'
solver.x0 = [0,0,0]

par.setpars([2.96,3,0.01,4])
par.bifpar(1,[3.42,3.72,301])
x0 = [0,1,0]
for p in par:
    # lyapunov exponents
    lp = lyapunov(sys,par,0,1000,5,x0)
    # number of turns
    solver.run()
    s = solver.solution()
    print p, lp, solver.nturns
