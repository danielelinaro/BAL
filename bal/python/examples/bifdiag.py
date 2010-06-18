#!/usr/bin/env python

from bal import *

outfile = 'hr.h5'

hr = DynamicalSystem()
hr.create('balHindmarshRose')

par = Parameters(hr.npar)
par.setpars([2.96,0.01,4],(0,2,3))
par.bifpar(1,[2.5,4.5,201])

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
