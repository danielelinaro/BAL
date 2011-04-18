#!/usr/bin/env python

from pybal import bal
from AUTOclui import *

# create a dynamical system and load the equations of the Hindmarsh-Rose neuron model
hr = bal.DynamicalSystem()
hr.create('balHindmarshRose')

# create the parameters of the model
par = bal.Parameters(hr.npar)
# set the fixed parameters
par.setpars([3,2.8,0.01,4])

# create an ODE solver
solver = bal.ODESolver(hr,par)
solver.x0 = [0,0,0]
solver.dt = 0.01
solver.ttran = 1000.0
solver.tstop = 2000.0
solver.mode = 'trajectory + events'
# integrate
solver.run()
# save the solution to file
solver.write_orbit('HR.dat')

# continue the 3-turns cycle and find a homoclinic orbit
T = 1500
r1 = run(e='HR', dat='HR', c='lc', NMX=5)
r2 = run(r1, ICP=[2,11,1], NMX=200, UZR={-2:2.95}, DS=0.01)
lc3 = run(r2('UZ1'), ICP=[1,11,2], DS=-0.01, ISP=0, UZR={-11:T}, NMX=2000, NPR=1000)
save(lc3,'lc.3')

cl()
