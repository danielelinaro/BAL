/*=========================================================================
 *
 *   Program:   Bifurcation Analysis Library
 *   Module:    balCommon.h
 *
 *   Copyright (C) 2009,2010 Daniele Linaro
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *=========================================================================*/

/**
 * \file balCommon.h
 * \brief Header file with common definitions for other bal classes.
 */

/*!
 * \mainpage A Bifurcation Analysis Library (BAL)
 *
 * \author Daniele Linaro <daniele.linaro@unige.it>
 * \version 0.9.5
 * \date 2009 - 2011
 * 
 * \section con Contents:
 *
 * -# \ref intro "Introduction"
 * -# \ref install "Installation"
 * -# \ref ex "Examples of use"
 * -# \ref links "Useful Links"
 *
 * \page intro Introduction
 * 
 * The study of ordinary differential equations (ODE) either autonomous or non autonomous,
 * generally nonlinear and depending on parameters, \f$\dot{x} = f (x,p,t)\f$,
 * requires good and powerful mathematical software.
 * One of the main concerns in this kind of studies is the so called bifurcation analysis,
 * that is the analysis of the qualitatively different asymptotic behaviors the dynamical
 * system can give rise to by varying one or more parameters.
 * 
 * Despite the several advantages offered by continuation methods, brute-force simulation
 * methods provide quite simply a global picture of the bifurcation scenario (at the cost
 * of extensive simulations) with respect to one or two parameters and so it's a crucial step
 * to start analyzing dynamics of an unknown system. \n
 * This is the reason why a bifurcation
 * scenario can be advantageously analyzed by combining simulations with continuation methods.
 *
 * This library aims to provide an easy way to compute brute-force bifurcation diagrams 
 * (by varying bifurcation parameters) and basins of attractions (by varying initial conditions),
 * trying, at the same time, to optimize the most common problems of brute-force methods:
 * time to solution (requested CPU time), memory use (efficiency of data storage), accuracy
 * of the results and flexibility.\n
 *
 * Bifurcation diagrams can be obtained by evaluating asymptotic trajectories
 * with a user-defined Poincare's section or by calculating their Lyapunov exponents
 * [see examples]. \n
 *
 * BAL has a multilayer organization with three external open-source C/C++ libraries at its
 * core: \a CVode for numerical integration, \a BOOST for extending
 * the functionality of C++ and handling parallel threads, and \a HDF5 for
 * efficient data storage.
 *
 * 
 * \page install Installation
 *
 * \section sub_1 I. Installation of auxiliary libraries
 * To install BAL, you need the following libraries:
 * -# \b CVode, for the solution of ordinary differentail equations.
 *		It is available <a href="https://computation.llnl.gov/casc/sundials/description/description.html">here</a>.
 * -# \b HDF5, to save data in compressed H5 files.
 *		It is available <a href="http://hdfgroup.org/HDF5/">here</a>.
 * -# \b BOOST, for multithreading operations.
 *		It is available <a href="http://www.boost.org">here</a>.
 * 
 * \remarks Please note that, if \b CVode, \b HDF5 and \b BOOST are not installed where the build
 * system can automatically find it, the user has to set the correct path to their \a include and \a lib
 * directory typing:
\verbatim
$ export CPPFLAGS="-I/path/to/the/include/directory"
$ export LDFLAGS"-L/path/to/the/libraries/directory" \endverbatim
 * For example, assuming that the libraries are installed in the directory
 * \a "local" in John's home directory, John should type:
\verbatim
$ export CPPFLAGS="-I/home/john/local/include"
$ export LDFLAGS="-L/home/john/local/lib" \endverbatim
 *	
 * \section sub_2 II. Installation Procedure 
 * The installation procedure is the usual one. Type:
\verbatim 
 $ ./configure [--with-examples] to compile also the examples
 $ make
 # make install\endverbatim
 *
 *
 * \page ex Examples of use
 * 
 * Here below are shown three typical applications of BAL library (for compilable example code
 * refer to examples section):
 * \n\n
 * - \ref bif_a "Bifurcation analysis"
 *		-# \ref sub1 "Classification of trajectories"
 *		-# \ref sub2 "Lyapunov exponents"
 * - \ref bas_a "Basins of attraction"
 * \n\n
 *
 * \section bif_a Bifurcation analysis
 *
 * Bifurcation analysis is performed on Rossler attractor (as implemeted in \ref bal::Rossler)
 * described by the following equations:
 * \f{eqnarray*}{
 \dot{x}_1&=&-x_2-x_3\\
 \dot{x}_2&=&x_1+ax_2\\
 \dot{x}_3&=&b+x_3(x_1-c)
 * \f}
 * with parameters set as follows: \f$a=0.25\f$, \f$b\in[0.2, 2.3]\f$ and \f$c\in [4.3, 9.6]\f$.\n
 * Parameters plane \f$(b,c)\f$ will be meshed with a \f$800\times800\f$ points regular grid.\n
 *
 * In this snippet of code are created and initialialized (as described above) the parameters for
 * bifurcation analysis:
 \code
 BifurcationParameters * bp = BifurcationParameters::Create();
 bp->SetNumber(3);
 bp->SetIthParameter(0,0.25);
 bp->SetIthParameterLowerBound(1,0.2);
 bp->SetIthParameterUpperBound(1,2.3);
 bp->SetIthParameterLowerBound(2,4.3);
 bp->SetIthParameterUpperBound(2,9.6);
 int steps[3] = {1,800,800};
 bp->SetNumberOfSteps(steps);
 \endcode
 *
 * then \ref bal::BifurcationParameters previously created is assigned to an instance 
 * of \ref bal::Rossler dynamical system: 
 \code
 Rossler * ros = Rossler::Create();
 ros->SetParameters(bp);
 \endcode
 * 
 * Last preliminar step is to create an instance of \ref bal::BifurcationDiagram which handles and
 * creates all the objects necessary to calculate a bifurcation diagram (\ref bal::ODESolver and
 * \ref bal::Logger)
 \code
 BifurcationDiagram * bifd = BifurcationDiagram::Create();
 bifd->SetDynamicalSystem(ros); 
 \endcode
 *
 * and the initial condition for the system:
 \code
 realtype x0[3] = {0.5,0.5,0.5};
 \endcode
 *
 * Now, by setting appropriately the fields of the solver, it is possible to obtain different types
 * of analysis, shown at \ref sub1 "I" and \ref sub2 "II".
 * \n\n
 * \subsection sub1 I. Bifurcation diagram from classification of trajectories
 * 
 * To be able to classify trajectories by number of turns, a valid Poincare' section must have been
 * implemented for the dynamical system; this can be achieved through:
 * - \ref bal::Rossler::Events describes the implicit equations of the surface.
 * - \ref bal::Rossler::EventsConstraints puts additional constraints on the "surface crossing" event.
 *		In this case disambiguates the intersection points by trajectory direction.
 * 
 * This section of code specifies the essential settings for the solver and data saving:
 \code
 bifd->GetODESolver()->SetIntegrationMode(BOTH);				// calculation of trajectories and events
 bifd->GetODESolver()->HaltAtEquilibrium(true);					// integration stops if trajectory converges to an equilibrium
 bifd->GetODESolver()->HaltAtCycle(true);						// integration stops if trajectory converges to a limit cycle
 bifd->GetODESolver()->SetTransientDuration(10e3);				// transient duration
 bifd->GetODESolver()->SetFinalTime(1e4);						// total integration time
 bifd->GetODESolver()->SetMaxNumberOfIntersections(500);		// maximum number of intersection with Poincare' section
 bifd->GetODESolver()->SetX0(x0);								// set initial condition 
 bifd->RestartFromX0(true);										// the solver use the same initial condition for each integration
 bifd->SetFilename("rossler.h5",true);							// trajectories will be saved in h5 compressed file format
 \endcode
 *
 * It must be set the number of integration threads used and start diagram computation:
 *
 \code
 bifd->SetNumberOfThreads(8);									// optimal value for a QuadCore CPU
 bifd->ComputeDiagram();										// computation starts
 \endcode
 *
 * Finally, classification file is saved: 
 *
 \code
 bifd->SaveSummaryData("rossler.classified");					// saves a file with data for bifurcation diagram (only parameters and #turns)
 \endcode
 *
 * \n The image shown below was obtained by imaging with MATLAB® the bifurcation diagram resulting
 * from data saved in file \a "rossler.classified".  
 *
 * \image html rossler_class.gif
 * \n\n
 *
 * \subsection sub2 II. Bifurcation diagram from Lyapunov exponents
 *
 * To perform a bifurcation analysis based on the value of the Maximum Lyapunov Exponent
 * is necessary to set up the solver as follows.\n
 * In this example, FinalTime and LyapunovTimeStep were choosen to have an accurate estimation
 * of the exponents at the expense of an increasing execution time.
 *
 \code
 bifd->GetODESolver()->SetIntegrationMode(LYAP);				// calculation of Lyapunov exponents
 bifd->GetODESolver()->SetTransientDuration(1e3);				// transient duration
 bifd->GetODESolver()->SetFinalTime(51e3);						// total integration time
 bifd->GetODESolver()->SetLyapunovTimeStep(50);					// time interval between each between each update procedure of exponents values
 bifd->GetODESolver()->SetTimeStep(0.5);						// time step of integration
 bifd->GetODESolver()->SetX0(x0);								// set initial condition
 \endcode 
 *
 * It's therefore possible to start computation as done previously:
 \code
 bifd->SetNumberOfThreads(8);									// optimal value for a QuadCore CPU
 bifd->ComputeDiagram();										// computation starts
 bifd->SaveSummaryData("rossler.lyap");							// saves a file with data for bifurcation diagram (only parameters and MLE)
 \endcode
 *
 * \n The image shown below was obtained by imaging with MATLAB® the bifurcation diagram resulting
 * from data saved in file \a "rossler.lyap".  
 *
 * \image html rossler_Lyap.gif
 *
 *
 * \remarks It's important to point out that, combining this two different types of bifurcation
 *					analysis supplied by BAL (as shown in figure below), it's possible to extract the maximum
 *					information obtainable from brute-force methods. The final bifurcation diagram can merge results
 *					that are inherently complementary: classification gives information on periodic solutions,
 *					while MLE allows to evaluate the "level" of chaos in chaotic regions.
 *					
 *					<br><br> \image html rossler_class+Lyap.gif <br><br>
 * 
 * \n\n\n
 * \section bas_a Basins of attraction
 *
 * In this example, it is shown how to calculate basins of attraction on a specific domain of the state space
 * of Hindmarsh-Rose dynamical system (as implemeted in \ref bal::HindmarshRose):
 * \f{eqnarray*}{
 \dot{x}_1&=&x_2-x_1^3+bx_1^2+I-x_3\\
 \dot{x}_2&=&1-5x_1^2-x_2\\
 \dot{x}_3&=&\mu(s(x_1-x_{rest})-x_3)
 * \f}
 * parameters will be set as follows: \f$b=2.88\f$, \f$I=2.6\f$, \f$\mu=0.01\f$ and \f$s=4\f$ (\f$x_{rest}=-1.6\f$ by default).\n
 * The choosen domain it's on a \f$x_2x_3\f$ plane with \f$x_1=1\f$ and \f$(x_2,x_3)\in\{[-10,2]\times[-1,2]\}\f$.
 * Each intervals of interest is divided in 6 subintervals, breaking the calculation in 36 blocks; each block is meshed with a
 * \f$100\times100\f$ points regular grid and results are saved in separated files, making data lighter and easier to manipulate.\n\n
 * Below are reported only the significative parts of code necessary to understand library functionalities.\n
 *
 * Definition of costants for domain's boundary and its subdivision:
 *
 \code
 #define	YBLOCKS	6				// number of blocks along y direction
 #define	ZBLOCKS	6				// number of blocks along z direction
 
 #define	YSTEPS	100				// number of points along y direction in each block
 #define	ZSTEPS	100				// number of points along z direction in each block
 
 #define	XMIN	1.
 #define	YMIN	-10.
 #define	YMAX	2.
 #define	ZMIN	-1.
 #define	ZMAX	2.
 \endcode
 *
 \code
 float y_blockstep = (YMAX-YMIN)/float(YBLOCKS);	// block spacing in y direction
 float z_blockstep = (ZMAX-ZMIN)/float(ZBLOCKS);	// block spacing in z direction
 float y_step = y_blockstep/float(YSTEPS);			// point spacing in y direction
 float z_step = z_blockstep/float(ZSTEPS);			// point spacing in z direction
 
 int nX0 = (YSTEPS+1)*(ZSTEPS+1);					// number of initial conditions in each block
 double **X0 = new double*[nX0];					// matrix for each block initial conditions
 for(i=0; i<nX0; i++)
	X0[i] = new double[3];
 \endcode
 *
 * Creation of \ref bal::HindmarshRose dynamical system and assignment of particular \ref bal::Parameters
 * value for which two stable solutions coexist:
 *
 \code
 HindmarshRose * hr = HindmarshRose::Create();
 
 Parameters * par = Parameters::Create();
 par->SetNumber(4);
 par->At(0) = 2.88;
 par->At(1) = 2.6;
 par->At(2) = 0.01;
 par->At(3) = 4.0;
 
 hr->SetParameters(par);
 \endcode
 *
 * Creation of a \ref bal::BifurcationDiagram instance to compute basins of attraction by setting operating mode
 * (\ref bal::diagram_mode ) to \a IC .
 *
 \code
 BifurcationDiagram * bifd = BifurcationDiagram::Create();
 bifd->SetDynamicalSystem(hr);
 
 bifd->SetMode(IC);
 \endcode
 *
 * Setting of solver options; in particular \ref bal::integration_mode is set to \a EVENTS to focus on final condition
 * and trajectory classification via Poincare's section:
 *
 \code
 bifd->GetODESolver()->SetIntegrationMode(EVENTS);
 
 bifd->GetODESolver()->HaltAtEquilibrium(true);
 bifd->GetODESolver()->HaltAtCycle(true);
 bifd->GetODESolver()->SetTransientDuration(5e3);
 bifd->GetODESolver()->SetFinalTime(1e5);
 bifd->GetODESolver()->SetMaxNumberOfIntersections(300);
 bifd->SetNumberOfThreads(16);
 \endcode
 *
 * Main code, iterating for each block:
 *
 \code
 for (by=0; by<YBLOCKS; by++) {
	y_start = YMIN + y_blockstep*by;
 
	for (bz=0; bz<ZBLOCKS; bz++) {
		z_start = ZMIN + z_blockstep*bz;
		
		//---- basin of attraction analysis for i-th block ---
 
		// creation and assignment of initial conditions matrix for i-th block
		for (i=0; i<ZSTEPS+1; i++) {
			for (j=0; j<YSTEPS+1; j++){ 
				index = j+(YSTEPS+1)*i;
				X0[index][0] = XMIN;
				X0[index][1] = y_start + y_step*j;
				X0[index][2] = z_start + z_step*i;
			}
		}
		bifd->SetInitialConditions(nX0,X0);
 
		// creating .h5 file for i-th block 	
		bifd->SetLogger(H5Logger::Create());
		sprintf(file_name,"[%d,%d]block.h5",by,bz);
		bifd->SetFilename(file_name,true); //compressed format
		
		bifd->ComputeDiagram();
		
		// saving summary data file for i-th block containing parameters
		// value, initial and final condition, number of turns
		sprintf(file_name,"[%d,%d]block.classified",by,bz);
		bifd->SaveSummaryData(file_name);
		
		//---------------------------------------------------
	}
 }
 \endcode
 *\n\n
 * The obtained results are reported below together with another basins plot calculated (applying small changes to the code)
 * in a complementary domain on \f$x_1x_2\f$ plane with \f$x_3=1.5\f$ and \f$(x_1,x_2)\in\{[0,2.5]\times[-10,2]\}\f$.\n
 * In third figure are plotted the two stable solution coexisting for the choosen parameters' values. 
 * <br><br>
 * \image html hr-basin.png
 * <br><br>
 * Solution, with 3 turns in each period, attracts the initial conditions displayed in dark grey tone. The other solution, with 4 turns in each period,
 * attracts the initial conditions displayed in light grey tone.\n
 * In both basin panels, the transient duration is ttran = 5000 and the integration is stopped as soon as a
 * periodic solution is detected. Since in this setup the classiﬁcation procedure never
 * fails, the integration is stopped after at most 8 intersections with the Poincare' section.
 *
 *
 *
 * \page links Links
 * - <a href="http://ncas.dibe.unige.it/research/lines/dynsys/">NCAS</a>
 * - ecc...
 *
 */

#ifndef _BALCOMMON_
#define _BALCOMMON_

#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#ifdef CVODE26
#include <sundials/sundials_direct.h>
#endif
#include <nvector/nvector_serial.h>

/** Ith numbers components 0..NEQ-1 */
#define Ith(v,i)    NV_Ith_S(v,i)
/** IJth numbers rows,cols 0..NEQ-1 */
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

/* colors */
#define ESC ''
#define RED "[31m"
#define GREEN "[32m"
#define YELLOW "[33m"
#define BLUE "[34m"
#define MAGENTA "[35m"
#define CYAN "[36m"
#define NORMAL "[00m"

#define PI				(3.141592653589793)

#define	LIST_MAX_SIZE	100

#include <exception>

namespace bal {

/**
 * \class Exception
 * \brief Class that implements the exceptions thrown by the BAL library
 */
class Exception : public std::exception {
 public:
 Exception(const char* description = NULL) : errorDescription(description) {}
  virtual const char* what() const throw() {
    if(errorDescription == NULL)
      return "bal::Exception";
    return errorDescription;
  }
 private:
  const char* errorDescription;
};

} // namespace bal

#endif

