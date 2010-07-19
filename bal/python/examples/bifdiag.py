#!/usr/bin/env python

from bal import *

outfile = 'hr.h5'

hr = DynamicalSystem()
hr.create('balHindmarshRose')

par = Parameters(hr.npar)
par.bifpar(0,[2.5,3.5,21])
par.bifpar(1,[2.5,4.5,101])
par.setpars([0.01,4],(2,3))

bifdiag = BifurcationDiagram(hr,par)
bifdiag.outfile = outfile
bifdiag.mode = 'events'
bifdiag.equilbreak = True
bifdiag.cyclebreak = False
bifdiag.ttran = 1e3
bifdiag.tstop = 1e4
bifdiag.intersections = 200
bifdiag.x0 = [0.5,0.5,0.5]
bifdiag.nthreads = 2

bifdiag.run()
bifdiag.classification('hr.dat')
