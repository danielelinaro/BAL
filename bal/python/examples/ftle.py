#!/usr/bin/env python

from pybal import bal
import numpy as np

outfile = 'gyre.h5'

gyre = bal.DynamicalSystem()
gyre.create('balDoubleGyre')

A = 0.1
omega = 2*np.pi/10
eps = 0.25
par = bal.Parameters(gyre.npar)
par.setpars([A,omega,eps])

bifdiag = bal.BifurcationDiagram(gyre,par)
bifdiag.outfile = outfile
bifdiag.mode = 'trajectory'
bifdiag.diagram_mode = 'ic'
bifdiag.equilbreak = False
bifdiag.cyclebreak = False
bifdiag.ttran = 15
bifdiag.tstop = 15
bifdiag.x0 = [[0.5,0.2],[1.5,0.8]]
bifdiag.nthreads = 2

bifdiag.run()
bifdiag.classification('gyre.classified')
