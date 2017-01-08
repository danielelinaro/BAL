
#AM_CPPFLAGS = @CVODEVERSION@ -DDEBUG -g -O0
AM_CPPFLAGS = @CVODEVERSION@ -O3

lib_LTLIBRARIES = libbal.la
libbal_la_SOURCES = balBifurcationParameters.cpp balDynamicalSystem.cpp balHindmarshRose.cpp \
	balLogger.cpp balObject.cpp balParameters.cpp balSolution.cpp balDynasty.cpp \
	balODESolver.cpp balBifurcationDiagram.cpp

#libbal_la_SOURCES += balInterpSystem.cpp balHeartNeuron.cpp balPLL.cpp balLorenz.cpp \
#	balRossler.cpp balDoubleGyre.cpp balInterp1D.cpp balInterp2D.cpp \
#	balInterp3D.cpp balInterpolator.cpp 

libbal_la_LDFLAGS = -version-info ${LIB_VER}

include_HEADERS = balBifurcationParameters.h balCommon.h balDynamicalSystem.h balHindmarshRose.h \
	balLogger.h balObject.h balParameters.h balSolution.h balDynasty.h balRandomLibrary.h \
	balODESolver.h balBifurcationDiagram.h

#include_HEADERS += balInterpSystem.h balInterp1D.h balInterp2D.h balInterp3D.h balInterpolator.h \
#	balHeartNeuron.h balPLL.h balLorenz.h balRossler.h balDoubleGyre.h
