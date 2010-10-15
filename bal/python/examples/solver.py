#!/usr/bin/env python

from pybal import bal
from pylab import figure, plot, xlabel, ylabel, title, show, axis

# create a dynamical system and load the equations of the Hindmarsh-Rose neuron model
hr = bal.DynamicalSystem()
hr.create('balHindmarshRose')

# create the parameters of the model
par = bal.Parameters(hr.npar)
# set the fixed parameters
par.setpars([2.96,0.01,4],(0,2,3))
# set the bifurcation parameter
par.bifpar(1,[2.5,4.5,11])

# create an ODE solver
solver = bal.ODESolver(hr,par)
solver.x0 = [0,0,0]
solver.dt = 0.01
solver.ttran = 1000.0
solver.tstop = 2000.0
solver.mode = 'trajectory + events'

# iterate over all possible tuple of parameters
for p in par:
    # integrate
    solver.run()
    # get the solution
    s = solver.solution()
    # plot the results of the integration
    figure()
    plot(s.data['x'][5::3],s.data['x'][3::3],'k')
    xlabel('t (a.u.)')
    ylabel('x (a.u.)')
    title('I = '+str(p[1])+' '+str(s.parameters[1])+' # turns = '+str(solver.nturns))
    axis('tight')
    show()
