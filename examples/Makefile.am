
LDADD = ../src/libbal.la

AM_CPPFLAGS = "@CVODEVERSION@" -I$(srcdir)/../src -I. -g -O0

if EXAMPLES
#noinst_PROGRAMS = params loggers dynsys bifparams solver bifdiag pll lorenz basin dynasty lorenzLyap hrLyap bifdiagLyap
#bin_PROGRAMS = params loggers dynsys bifparams solver bifdiag pll lorenz basin dynasty lorenzLyap hrLyap bifdiagLyap rosslerLyap
bin_PROGRAMS = params loggers dynsys bifparams solver
params_SOURCES = params.cpp 
loggers_SOURCES = loggers.cpp 
dynsys_SOURCES = dynsys.cpp 
bifparams_SOURCES = bifparams.cpp 
solver_SOURCES = solver.cpp 
#lorenz_SOURCES = lorenz.cpp 
#lorenzLyap_SOURCES = lorenzLyap.cpp
#hrLyap_SOURCES = hrLyap.cpp
#bifdiag_SOURCES = bifdiag.cpp
#bifdiagLyap_SOURCES = bifdiagLyap.cpp
#basin_SOURCES = basin.cpp 
#dynasty_SOURCES = dynasty.cpp 
#pll_SOURCES = pll.cpp
#rosslerLyap_SOURCES = rosslerLyap.cpp
endif

