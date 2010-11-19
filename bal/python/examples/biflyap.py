#!/usr/bin/env python

from pybal import bal

hr = bal.DynamicalSystem()
hr.create('balHindmarshRose')

par = bal.Parameters(hr.npar)
par.bifpar(0,[2.5,3.5,11])
par.bifpar(1,[2.5,4.5,21])
par.setpars([0.01,4],(2,3))

bifdiag = bal.BifurcationDiagram(hr,par)
bifdiag.mode = 'lyapunov'
print bifdiag.mode
bifdiag.ttran = 5e2
bifdiag.tstop = 2e3
bifdiag.dt = 1
bifdiag.x0 = [0.5,0.5,0.5]
bifdiag.nthreads = 2

bifdiag.run()
bifdiag.classification('hr.lyap')
