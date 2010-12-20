
LDADD = ../src/libbal.la

# debug off
AM_CPPFLAGS = "@CVODEVERSION@" -I$(srcdir)/../src -I.

if EXAMPLES
#noinst_PROGRAMS = obj params loggers dynsys bifparams solver bifdiag interp pll lorenz interpsystem basin dynasty lorenzLyap hrLyap bifdiagLyap
bin_PROGRAMS = obj params loggers dynsys bifparams solver bifdiag interp pll lorenz interpsystem basin dynasty lorenzLyap hrLyap bifdiagLyap
interp_SOURCES = interp.cpp 
obj_SOURCES = obj.cpp 
params_SOURCES = params.cpp 
loggers_SOURCES = loggers.cpp 
dynsys_SOURCES = dynsys.cpp 
bifparams_SOURCES = bifparams.cpp 
solver_SOURCES = solver.cpp 
lorenz_SOURCES = lorenz.cpp 
lorenzLyap_SOURCES = lorenzLyap.cpp
hrLyap_SOURCES = hrLyap.cpp
interpsystem_SOURCES = interpsystem.cpp 
bifdiag_SOURCES = bifdiag.cpp
bifdiagLyap_SOURCES = bifdiagLyap.cpp
basin_SOURCES = basin.cpp 
dynasty_SOURCES = dynasty.cpp 
pll_SOURCES = pll.cpp ConfigFile.cpp
endif

