
AM_CPPFLAGS = "@CVODEVERSION@"

lib_LTLIBRARIES = libbal.la
libbal_la_SOURCES = balBifurcationDiagram.cpp balBifurcationDiagram.h balBifurcationParameters.cpp balBifurcationParameters.h \
balCommon.h balDynamicalSystem.cpp balDynamicalSystem.h balHindmarshRose.cpp balHindmarshRose.h balHeartNeuron.cpp balHeartNeuron.h \
balLogger.cpp balLogger.h balODESolver.cpp balODESolver.h balObject.cpp balObject.h balParameters.cpp balParameters.h balPLL.cpp balPLL.h \
balInterp1D.cpp balInterp1D.h balInterp2D.cpp balInterp2D.h balSolution.cpp balSolution.h \
balLorenz.cpp balLorenz.h balEye.cpp balEye.h
libbal_la_LDFLAGS = -version-info ${LIB_VER}

include_HEADERS = balBifurcationDiagram.h balBifurcationParameters.h balCommon.h balDynamicalSystem.h balHindmarshRose.h \
balHeartNeuron.h balLogger.h balODESolver.h balObject.h balParameters.h balPLL.h balInterp1D.h balInterp2D.h balSolution.h \
balLorenz.h

