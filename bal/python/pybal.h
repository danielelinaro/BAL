#ifndef PYBAL
#define PYBAL

#include <Python.h>
#include <structmember.h>

#include <sundials/sundials_types.h>
//#include <nvector/nvector_serial.h>
//#include <sundials/sundials_dense.h>
//#include <cvode/cvode.h>

#include "balObject.h"
#include "balCommon.h"
#include "balDynamicalSystem.h"
#include "balParameters.h"
#include "balBifurcationParameters.h"
#include "balODESolver.h"
#include "balBifurcationDiagram.h"

#include <dlfcn.h>

#ifdef __cplusplus
extern "C" {
#endif

/******************** pyBalDynamicalSystem ********************/

typedef balDynamicalSystem* (*Factory)();
	
typedef struct {
	PyObject_HEAD
	balDynamicalSystem * dynsys;
	void * lib;
} pyBalDynamicalSystem;

static PyObject * pyBalDynamicalSystem_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int pyBalDynamicalSystem_init(pyBalDynamicalSystem *self, PyObject *args, PyObject *kwds);
static void pyBalDynamicalSystem_dealloc(pyBalDynamicalSystem *self);
static PyObject * pyBalDynamicalSystem_getattro(pyBalDynamicalSystem *self, PyObject *name);
static int pyBalDynamicalSystem_setattro(pyBalDynamicalSystem *self, PyObject *name, PyObject *value);
static PyObject * pyBalDynamicalSystem_name(pyBalDynamicalSystem *self);
static PyObject * pyBalDynamicalSystem_create(PyObject *self, PyObject *args);

/******************** pyBalParameters ********************/

typedef struct {
	PyObject_HEAD
	balBifurcationParameters * bifparams;
} pyBalParameters;
	
static PyObject * pyBalParameters_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int pyBalParameters_init(pyBalParameters *self, PyObject *args, PyObject *kwds);
static void pyBalParameters_dealloc(pyBalParameters *self);
static PyObject * pyBalParameters_getattro(pyBalParameters *self, PyObject *name);
static int pyBalParameters_setattro(pyBalParameters *self, PyObject *name, PyObject *value);
static PyObject * pyBalParameters_par(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject * pyBalParameters_bifpar(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject * pyBalParameters_setpars(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject * pyBalParameters_hastuples(pyBalParameters *self);
static PyObject * pyBalParameters_hasnext(pyBalParameters *self);
static PyObject * pyBalParameters_next(pyBalParameters *self);
static PyObject * pyBalParameters_reset(pyBalParameters *self);
static PyObject * pyBalParameters_iter(pyBalParameters *self);
static PyObject * pyBalParameters_iternext(pyBalParameters *self);

/******************** pyBalODESolver ********************/
	
typedef struct {
	PyObject_HEAD
	balODESolver * solver;
	pyBalDynamicalSystem * dynsys;
	pyBalParameters *params;
} pyBalODESolver;
	
static PyObject * pyBalODESolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int pyBalODESolver_init(pyBalODESolver *self, PyObject *args, PyObject *kwds);
static void pyBalODESolver_dealloc(pyBalODESolver *self);
static PyObject * pyBalODESolver_getattro(pyBalODESolver *self, PyObject *name);
static int pyBalODESolver_setattro(pyBalODESolver *self, PyObject *name, PyObject *value);
static PyObject * pyBalODESolver_solve(pyBalODESolver *self);

/******************** pyBalSolution ********************/
	
typedef struct {
	PyObject_HEAD
	balParameters * params;
	PyObject * data;
} pyBalSolution;

static pyBalSolution* createPyBalSolution(balSolution *s);
static void pyBalSolution_dealloc(pyBalSolution *self);
static PyObject * pyBalSolution_getattro(pyBalSolution *self, PyObject *name);

/******************** pyBalBifurcationDiagram ********************/

typedef struct {
	PyObject_HEAD
	balBifurcationDiagram * diagram;
	pyBalDynamicalSystem * dynsys;
	pyBalParameters *params;
} pyBalBifurcationDiagram;

static PyObject * pyBalBifurcationDiagram_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int pyBalBifurcationDiagram_init(pyBalBifurcationDiagram *self, PyObject *args, PyObject *kwds);
static void pyBalBifurcationDiagram_dealloc(pyBalBifurcationDiagram *self);
static PyObject * pyBalBifurcationDiagram_getattro(pyBalBifurcationDiagram *self, PyObject *name);
static int pyBalBifurcationDiagram_setattro(pyBalBifurcationDiagram *self, PyObject *name, PyObject *value);
static PyObject * pyBalBifurcationDiagram_compute(pyBalBifurcationDiagram *self);
static PyObject * pyBalBifurcationDiagram_classification(pyBalBifurcationDiagram *self, PyObject *args, PyObject *kwds);
	
/******************** Module stuff ********************/

#define BALLIB "libbal.0.0.0.dylib"
#define BALLIBEXT "libbalext.so"
extern void * ballib;
extern void * ballibext;

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initbal(void);
	
	
#ifdef __cplusplus
}
#endif

#endif
