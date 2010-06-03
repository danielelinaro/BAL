
LDADD = ../src/libbal.la
AM_CPPFLAGS = -I$(srcdir)/../src

if EXAMPLES
noinst_PROGRAMS = obj params loggers dynsys bifparams solver bifdiag interp pll
interp_SOURCES = interp.cpp 
obj_SOURCES = obj.cpp 
params_SOURCES = params.cpp 
loggers_SOURCES = loggers.cpp 
dynsys_SOURCES = dynsys.cpp 
bifparams_SOURCES = bifparams.cpp 
solver_SOURCES = solver.cpp 
bifdiag_SOURCES = bifdiag.cpp 
pll_SOURCES = pll.cpp ConfigFile.cpp
endif

