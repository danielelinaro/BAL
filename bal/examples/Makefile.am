
LDADD = ../src/libbal.la

AM_CPPFLAGS = "@CVODEVERSION@" -I$(srcdir)/../src -I.

if EXAMPLES
noinst_PROGRAMS = obj params loggers dynsys bifparams solver bifdiag interp pll lorenz eye basin
interp_SOURCES = interp.cpp 
obj_SOURCES = obj.cpp 
params_SOURCES = params.cpp 
loggers_SOURCES = loggers.cpp 
dynsys_SOURCES = dynsys.cpp 
bifparams_SOURCES = bifparams.cpp 
solver_SOURCES = solver.cpp 
lorenz_SOURCES = lorenz.cpp 
eye_SOURCES = eye.cpp 
bifdiag_SOURCES = bifdiag.cpp 
basin_SOURCES = basin.cpp 
pll_SOURCES = pll.cpp ConfigFile.cpp
endif

