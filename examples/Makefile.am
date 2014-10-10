
LDADD = ../src/libbal.la

AM_CPPFLAGS = "@CVODEVERSION@" -I$(srcdir)/../src -I. -g3 -O0 -DDEBUG
AM_CFLAGS = -g3 -O0
AM_CXXFLAGS = -g3 -O0

if EXAMPLES
#noinst_PROGRAMS = params loggers dynsys bifparams solver bifdiag pll lorenz basin dynasty lorenzLyap hrLyap bifdiagLyap
bin_PROGRAMS = params dynsys bifparams loggers solver lorenz lorenzLyap hrLyap
params_SOURCES = params.cpp 
loggers_SOURCES = loggers.cpp 
dynsys_SOURCES = dynsys.cpp 
bifparams_SOURCES = bifparams.cpp 
solver_SOURCES = solver.cpp 
lorenz_SOURCES = lorenz.cpp 
lorenzLyap_SOURCES = lorenzLyap.cpp
hrLyap_SOURCES = hrLyap.cpp
#bifdiag_SOURCES = bifdiag.cpp
#bifdiagLyap_SOURCES = bifdiagLyap.cpp
#basin_SOURCES = basin.cpp 
#dynasty_SOURCES = dynasty.cpp 
#pll_SOURCES = pll.cpp
endif

