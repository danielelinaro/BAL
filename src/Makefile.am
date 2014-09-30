
AM_CPPFLAGS = @CVODEVERSION@ -DDEBUG -g -O0

lib_LTLIBRARIES = libbal.la
libbal_la_SOURCES = balBifurcationParameters.cpp balBifurcationParameters.h \
balCommon.h balDynamicalSystem.cpp balDynamicalSystem.h balHindmarshRose.cpp balHindmarshRose.h balHeartNeuron.cpp balHeartNeuron.h \
balLogger.cpp balLogger.h balODESolver.cpp balODESolver.h balObject.cpp balObject.h balParameters.cpp balParameters.h balPLL.cpp balPLL.h \
balInterp1D.cpp balInterp1D.h balInterp2D.cpp balInterp2D.h balInterp3D.cpp balInterp3D.h balInterpolator.cpp balInterpolator.h balSolution.cpp balSolution.h \
balLorenz.cpp balLorenz.h balInterpSystem.cpp balInterpSystem.h balDynasty.cpp balDynasty.h balRossler.cpp balRossler.h \
balDoubleGyre.cpp balDoubleGyre.h
libbal_la_LDFLAGS = -version-info ${LIB_VER}

include_HEADERS = balBifurcationDiagram.h balBifurcationParameters.h balCommon.h balDynamicalSystem.h balInterpSystem.h balHindmarshRose.h \
balHeartNeuron.h balLogger.h balODESolver.h balObject.h balParameters.h balPLL.h balInterp1D.h balInterp2D.h balInterp3D.h balInterpolator.h balSolution.h \
balLorenz.h balDynasty.h balRossler.h balRandomLibrary.h balDoubleGyre.h

