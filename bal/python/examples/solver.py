#!/usr/bin/env python

from pybal import bal
from pylab import figure, plot, xlabel, ylabel, title, show, axis

hr = bal.DynamicalSystem()
hr.create('balHindmarshRose')
par = bal.Parameters(hr.npar)
par.setpars([2.96,0.01,4],(0,2,3))
par.bifpar(1,[2.5,4.5,11])
solver = bal.ODESolver(hr,par)
solver.x0 = [0,0,0]
solver.dt = 0.01
solver.ttran = 500.0
solver.tstop = 2000.0
solver.mode = 'trajectory + events'

for p in par:
    solver.run()
    s = solver.solution()
    figure()
    plot(s.data['x'][5::3],s.data['x'][3::3],'k')
    xlabel('t (a.u.)')
    ylabel('x (a.u.)')
    title('I = '+str(p[1])+' '+str(s.parameters[1])+' # turns = '+str(solver.nturns))
    axis('tight')
    show()
