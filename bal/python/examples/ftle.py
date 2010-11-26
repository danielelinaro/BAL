#!/usr/bin/env python

from pybal import bal
import numpy as np

def createGrid(xlim,ylim,n,delta):
    dx = np.diff(xlim)/(n[0]-1)
    dy = np.diff(ylim)/(n[1]-1)
    N = np.prod(n)
    X,Y = np.mgrid[xlim[0]:xlim[1]+dx/2:dx,ylim[0]:ylim[1]+dy/2:dy]
    X = X.flatten()
    Y = Y.flatten()
    XY = np.zeros((3*N,2))
    for k in range(N):
        XY[3*k,:] = [X[k],Y[k]]
        if X[k]+delta[0] > xlim[1]:
            XY[3*k+1,:] = [X[k]-delta[0],Y[k]]
        else:   
            XY[3*k+1,:] = [X[k]+delta[0],Y[k]]
        if Y[k]+delta[1] > ylim[1]:
            XY[3*k+2,:] = [X[k],Y[k]-delta[1]]
        else:
            XY[3*k+2,:] = [X[k],Y[k]+delta[1]]
    return XY

outfile = 'gyre.h5'

gyre = bal.DynamicalSystem()
gyre.create('balDoubleGyre')

A = 0.1
omega = 2*np.pi/10
eps = 0.25
par = bal.Parameters(gyre.npar)
par.setpars([A,omega,eps])

xlim = [0.,2.]
ylim = [0.,1.]
n = [11,21]
delta = [.01,.01]
X0 = createGrid(xlim,ylim,n,delta)
bifdiag = bal.BifurcationDiagram(gyre,par)
bifdiag.outfile = outfile
bifdiag.mode = 'trajectory'
bifdiag.diagram_mode = 'ic'
bifdiag.equilbreak = False
bifdiag.cyclebreak = False
bifdiag.ttran = 0
bifdiag.tstop = 20
bifdiag.dt = 0.05
bifdiag.x0 = X0.tolist()
bifdiag.nthreads = 2

bifdiag.run()
#data = bifdiag.summary()
#X1 = np.array(data)
#del data
