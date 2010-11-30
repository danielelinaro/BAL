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

def computeFTLE(X0,X1,T):
    dx = X0[1,0] - X0[0,0]
    dy = X0[2,1] - X0[0,1]
    J = np.array([[(X1[1,0]-X1[0,0])/dx,(X1[2,0]-X1[0,0])/dy],\
                  [(X1[1,1]-X1[0,1])/dx,(X1[2,1]-X1[0,1])/dy]])
    #print J
    delta = np.dot(J.transpose(),J)
    w,v = np.linalg.eig(delta)
    ftle = 1./T * np.log(np.sqrt(np.max(w)))
    return ftle

outfile = 'gyre.h5'
T = 15

gyre = bal.DynamicalSystem()
gyre.create('balDoubleGyre')

A = 0.25
omega = 2*np.pi/10
eps = 0.1
par = bal.Parameters(gyre.npar)
par.setpars([A,omega,eps])

xlim = [0.,2.]
ylim = [0.,1.]
#n = [201,101]
n = [5,5]
N = np.prod(n)
delta = [.001,.001]
X0 = createGrid(xlim,ylim,n,delta)
bifdiag = bal.BifurcationDiagram(gyre,par)
bifdiag.outfile = outfile
bifdiag.mode = 'trajectory'
bifdiag.diagram_mode = 'ic'
bifdiag.equilbreak = False
bifdiag.cyclebreak = False
bifdiag.ttran = 0
bifdiag.tstop = T
bifdiag.dt = 0.05
bifdiag.x0 = X0.tolist()
bifdiag.nthreads = 2

bifdiag.run()
bifdiag.summary('gyre.ic')
data = np.loadtxt('gyre.ic')
X1 = data[:,5:7]

#print('Computing finite time Lyapunov exponents...')
#ftle = np.zeros(N)
#for k in range(N):
#    ftle[k] = computeFTLE(X0[3*k:3*(k+1)],X1[3*k:3*(k+1)],T)
    
#buffer = np.zeros((N,3))
#buffer[:,0:2] = X0[::3]
#buffer[:,2] = ftle
#np.savetxt('ftle.dat',buffer,'%.6f')
