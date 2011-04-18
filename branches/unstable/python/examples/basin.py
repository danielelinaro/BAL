#!/usr/bin/env python

from pybal import bal

outfile = 'hr.h5'

hr = bal.DynamicalSystem()
hr.create('HindmarshRose')

par = bal.Parameters(hr.npar)
par.setpars([2.88,2.6,0.01,4])

bifdiag = bal.BifurcationDiagram(hr,par)
bifdiag.outfile = outfile
bifdiag.mode = 'events'
bifdiag.diagram_mode = 'ic'
bifdiag.equilbreak = True
bifdiag.cyclebreak = True
bifdiag.ttran = 1e3
bifdiag.tstop = 1e4
bifdiag.intersections = 300
bifdiag.x0 = [[2,-4,1.5],[0,-4,1.5]]
bifdiag.nthreads = 2

bifdiag.run()
bifdiag.summary('hr.classified')
