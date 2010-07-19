#!/usr/bin/env python

### FUNCTIONS - START ###

def printHelp(scriptname):
    print('\nUsage: ' + scriptname + ' [options]\n')
    print('\twhere options are:\n');
    print('\t\t-f --vector-file: read the specified file.\n')
    print('\t\t-h --help: print this help message.\n')
    print('\nAuthor: Daniele Linaro -- daniele.linaro@unige.it\n')

#### FUNCTIONS - END ####

from sys import argv
from os.path import isfile

if len(argv) == 2 and (argv[1] == '-h' or argv[1] == '--help'):
    printHelp(argv[0])
    exit(1)

if len(argv) > 2 and (argv[1] == '-f' or argv[1] == '--vector-file'):
    vectorfile = argv[2]
else:
    printHelp(argv[0])
    exit(1)

if not isfile(vectorfile):
    printHelp(argv[0])
    print(vectorfile + ' does not exist. Aborting...')
    exit(1)

from pybal import bal

# the dynamical system
eye = bal.DynamicalSystem()
eye.create('balEye')
eye.options = vectorfile

# the parameters
par = bal.Parameters(eye.npar)

# the solver
solver = bal.ODESolver(eye,par)
solver.x0 = [-1,1]
solver.intersections = 1e7
solver.dt = 1e-2
solver.ttran = 0
solver.tstop = 10
solver.mode = 'trajectory'

solver.run()
#s = solver.solution()
#saveH5file([s],'eye.h5')

#from pylab import figure, plot, xlabel, ylabel, title, show, axis

#figure()
#plot(s.data['t'],s.data['x'][0::2],'k')
#xlabel('t')
#ylabel('x')
#axis('tight')
#show()
