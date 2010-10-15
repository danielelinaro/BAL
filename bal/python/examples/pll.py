#!/usr/bin/env python

### FUNCTIONS - START ###

def printHelp(scriptname):
    print('\nUsage: ' + scriptname + ' [options]\n')
    print('\twhere options are:\n');
    print('\t\t-f --conf-file: parse the specified configuration file.\n')
    print('\t\t-h --help: print this help message.\n')
    print('\nAuthor: Daniele Linaro -- daniele.linaro@unige.it\n')

#### FUNCTIONS - END ####

from sys import argv
from os.path import isfile

if len(argv) == 2 and (argv[1] == '-h' or argv[1] == '--help'):
    printHelp(argv[0])
    exit(1)

conffile = 'pll.cfg'
if len(argv) > 2:
    if argv[1] == '-f' or argv[1] == '--config-file':
        conffile = argv[2]
    else:
        printHelp(argv[0])
        exit(1)

from ConfigParser import ConfigParser
if not isfile(conffile):
    printHelp(argv[0])
    print(conffile + ' is not a valid configuration file. Aborting...')
    exit(1)

fid = open(conffile,'r')
config = ConfigParser()
config.readfp(fid)
fid.close()

from pybal import bal
from pybal import util

# the dynamical system
pll = bal.DynamicalSystem()
pll.create('balPLL')

# the parameters
par = bal.Parameters(pll.npar)
parnames = ['fref','r1','fvco','vdd','rho0','rhoap','k0','krho','kap','alpha','kvcoa','kvcob','kvcoc','tuning']
for k,p in enumerate(parnames):
    steps = config.getint(p,'steps')
    if steps == 1:
        pmin = config.getfloat(p,'value')
        pmax = pmin
    elif steps > 1:
        pmin = config.getfloat(p,'min')
        pmax = config.getfloat(p,'max')
    else:
        print('The number of steps for ' + p + ' must be equal to or greater than 1.')
        exit(1)
    par.bifpar(k,[pmin,pmax,steps])
    if p == 'vdd':
        vdd = pmin

# the solver
solver = bal.ODESolver(pll,par)
solver.x0 = [0,vdd,0,0]
solver.intersections = 1e7
solver.dt = 5e-11
solver.ttran = config.getfloat('Simulation','ttran')
solver.tstop = config.getfloat('Simulation','tout')
if config.getint('Simulation','trajectory'):
    solver.mode = 'trajectory + events'
else:
    solver.mode = 'events'

print('Starting simulation...')
solver.run()
print('Simulation finished...')
s = solver.solution()
util.saveH5file([s],'pll.h5')

from pylab import figure, plot, xlabel, ylabel, title, show, axis

figure()
plot(s.data['t'],s.data['x'][3::4],'k')
xlabel('t (s)')
ylabel('w (V)')
axis('tight')
show()

