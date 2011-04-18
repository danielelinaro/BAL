#include "pybal.h"
//#include "balEye.h"

// libraries loaded at runtime
void * ballib;
void * ballibext;

/******************** pyBalDynamicalSystem ********************/

static PyObject * pyBalDynamicalSystem_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	pyBalDynamicalSystem *self;
	self = (pyBalDynamicalSystem *) type->tp_alloc(type,0);
	self->dynsys = NULL;
	return (PyObject *) self;
}

static int pyBalDynamicalSystem_init(pyBalDynamicalSystem *self, PyObject *args, PyObject *kwds) {
	self->dynsys = bal::DynamicalSystem::Create();
	self->lib = NULL;
	return 0;
}

static void pyBalDynamicalSystem_dealloc(pyBalDynamicalSystem *self) {
	if(self->dynsys != NULL)
		self->dynsys->Destroy();
	if(self->lib != NULL)
		dlclose(self->lib);
	self->ob_type->tp_free((PyObject *) self);
}

static PyObject * pyBalDynamicalSystem_name(pyBalDynamicalSystem *self) {
	return Py_BuildValue("s",self->dynsys->GetClassName());
}

static PyObject * pyBalDynamicalSystem_create(PyObject *self, PyObject *args) {
	const char *name;
	if (!PyArg_ParseTuple(args, "s", &name))
		return Py_BuildValue("i",1);
	
	void * lib = NULL;
	char libname[50], symbolname[50];
	sprintf(symbolname, "%sFactory", name);
	void * addr = dlsym(ballib, symbolname);
	if(addr == NULL) {
		if(ballibext != NULL) {
			addr = dlsym(ballibext, symbolname);
		}
		if(ballibext == NULL || addr == NULL) {
			sprintf(libname, "%s.so", name);
			lib = dlopen(libname, RTLD_LAZY);
			if(lib == NULL) {
				printf("Unable to locate the symbol %s.\n", symbolname);
				return Py_BuildValue("i",1);
			}
			addr = dlsym(lib, symbolname);
			if(addr == NULL) {
				printf("Unable to locate the symbol %s.\n", symbolname);
				dlclose(lib);
				return Py_BuildValue("i",1);
			}
		}
	}

	// +++
	Py_INCREF(self);
	
	// cast self to the appropriate type
	pyBalDynamicalSystem *ds = (pyBalDynamicalSystem *) self;
	// destroy the previously allocated dynamical system
	ds->dynsys->Destroy();
	// allocate a new dynamical system
	Factory f = (Factory) addr;
	ds->dynsys = f();
	if(lib != NULL) {
		//close the previous library (if there is one...)
		if(ds->lib != NULL)
			dlclose(ds->lib);
		//copy the pointer to the new library
		ds->lib = lib;
	}
	
	// ---
	Py_DECREF(self);
	return Py_BuildValue("i",0);
}

static PyObject * pyBalDynamicalSystem_special(PyObject *self, PyObject *args) {
	const char *opt;
	if (!PyArg_ParseTuple(args, "s", &opt))
		return Py_BuildValue("i",1);

	// +++
	Py_INCREF(self);
	// cast self to the appropriate type
	pyBalDynamicalSystem *ds = (pyBalDynamicalSystem *) self;
	bool flag = ds->dynsys->SpecialOptions((void *) opt);
	int retval = (flag ? 0 : 1);
	// ---
	Py_DECREF(self);

	return Py_BuildValue("i",0);
}

static PyObject * pyBalDynamicalSystem_getattro(pyBalDynamicalSystem *self, PyObject *name) {
	Py_INCREF(name);
	char *n = PyString_AsString(name);
	PyObject *result = NULL;
	if (strcmp(n, "ndim") == 0)
		result = Py_BuildValue("i",self->dynsys->GetDimension());
	else if(strcmp(n, "nev") == 0)
		result = Py_BuildValue("i",self->dynsys->GetNumberOfEvents());
	else if(strcmp(n, "npar") == 0)
		result = Py_BuildValue("i",self->dynsys->GetNumberOfParameters());
	else if (strcmp(n, "extended") == 0)
		result = Py_BuildValue("i",self->dynsys->IsExtended() ? 1 : 0);
	else
		result = PyObject_GenericGetAttr((PyObject*)self, name);
	Py_DECREF(name);
	return result;
}

static int pyBalDynamicalSystem_setattro(pyBalDynamicalSystem *self, PyObject *name, PyObject *value) {
	int err = 0;
	Py_INCREF(name);
	char *n = PyString_AsString(name);
	
	if (strcmp(n, "extended") == 0) {
		int extend;
		if(PyArg_Parse(value,"i",&extend)) {
			self->dynsys->Extend(extend);
		}
		else {
			PyErr_SetString(PyExc_ValueError,"Wrong value.");
			err = -1;
		}
	}
	else if (strcmp(n, "options") == 0) {
		if(PyString_Check(value)) {
			char *opt = PyString_AsString(value);
			printf("before calling SpecialOptions with opt = '%s'\n", opt);
			bool retval = false;
			retval = self->dynsys->SpecialOptions((void *) opt);
			/*
			Py_INCREF(self->dynsys);
			bal::Eye *eye = (bal::Eye *) self->dynsys;
			retval = eye->ReadVectorField(opt);
			Py_DECREF(self->dynsys);
			 */
			printf("after calling SpecialOptions, which returned %s\n", (retval ? "true" : "false"));
		}
		else {
			printf("Unknown option...\n");
			PyErr_SetString(PyExc_ValueError,"Wrong option.");
			err = -1;
		}
	}
	else {
		err = PyObject_GenericSetAttr((PyObject*)self, name, value);
	}
	Py_DECREF(name);
	return err;
}

static PyMethodDef pyBalDynamicalSystem_methods[] = {
	{"name", (PyCFunction) pyBalDynamicalSystem_name, METH_NOARGS, "Return the name of the class"},
	{"create", (PyCFunction) pyBalDynamicalSystem_create, METH_VARARGS, "Create a particular dynamical system"},
	{"special", (PyCFunction) pyBalDynamicalSystem_special, METH_VARARGS, "Sets special options"},
	{NULL}	/* Sentinel */
};

static PyTypeObject pyBalDynamicalSystemType = {
	PyObject_HEAD_INIT(NULL)
	0,									/* ob_size */
	"bal.DynamicalSystem",			/* ob_name */
	sizeof(pyBalDynamicalSystem),		/* tp_basicsize */
	0,									/* tp_itemsize */
	(destructor)pyBalDynamicalSystem_dealloc,	/* tp_dealloc */
	0,									/* tp_print */
	0,									/* tp_getattr */
	0,									/* tp_setattr */
	0,									/* tp_compare */
	0,									/* tp_repr */
	0,									/* tp_as_number */
	0,									/* tp_as_sequence */
	0,									/* tp_as_mapping */
	0,									/* tp_hash */
	0,									/* tp_call */
	0,									/* tp_str */
	(getattrofunc)pyBalDynamicalSystem_getattro,	/* tp_getattro  */
	(setattrofunc)pyBalDynamicalSystem_setattro,	/* tp_setattro  */
	0,									/* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/* tp_flags */
	"Dynamical System object",			/* tp_doc */
	0,									/* tp_traverse */
    0,									/* tp_clear */
    0,									/* tp_richcompare */
    0,									/* tp_weaklistoffset */
    0,									/* tp_iter */
    0,									/* tp_iternext */
    pyBalDynamicalSystem_methods,		/* tp_methods */
    0,									/* tp_members */
    0,									/* tp_getset */
    0,									/* tp_base */
    0,									/* tp_dict */
    0,									/* tp_descr_get */
    0,									/* tp_descr_set */
    0,									/* tp_dictoffset */
    (initproc)pyBalDynamicalSystem_init,/* tp_init */
    0,									/* tp_alloc */
    pyBalDynamicalSystem_new			/* tp_new */
};

/******************** pyBalParameters ********************/

static PyObject * pyBalParameters_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	pyBalParameters *self;
	self = (pyBalParameters *) type->tp_alloc(type,0);
	self->bifparams = NULL;
	return (PyObject *) self;
}

static int pyBalParameters_init(pyBalParameters *self, PyObject *args, PyObject *kwds) {
	int npar = 0;
	static char *kwlist[] = {"number", NULL};
	
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "|i", kwlist, &npar))
		return NULL;
	if (npar < 0) {
		printf("The number of parameters must be greater than or equal to zero.\n");
		return NULL;
	}
	
	self->bifparams = bal::BifurcationParameters::Create();
	self->bifparams->SetNumber(npar);
	return 0;
}

static void pyBalParameters_dealloc(pyBalParameters *self) {
	if(self->bifparams != NULL)
		self->bifparams->Destroy();
	self->ob_type->tp_free((PyObject *) self);
}

static PyObject * pyBalParameters_getattro(pyBalParameters *self, PyObject *name) {
	Py_INCREF(name);
	char *n = PyString_AsString(name);
	PyObject *result = NULL;

	int i, npar = self->bifparams->GetNumber();

	if (strcmp(n, "npars") == 0) {
		result = Py_BuildValue("i",npar);
	}
	else if(strcmp(n, "steps") == 0) {
		result = PyTuple_New((Py_ssize_t) npar);
		for(i=0; i<npar; i++)
			PyTuple_SET_ITEM(result,i,Py_BuildValue("i",self->bifparams->GetNumberOfSteps(i)));
	}
	else if(strcmp(n, "pars") == 0) {
		result = PyList_New((Py_ssize_t) npar);
		for(i=0; i<npar; i++)
			PyList_SET_ITEM(result,i,Py_BuildValue("d",self->bifparams->At(i)));
	}
	else if(strcmp(n, "bifpars") == 0) {
		PyObject *lower = PyList_New((Py_ssize_t) npar);
		PyObject *upper = PyList_New((Py_ssize_t) npar);
		PyObject *steps = PyTuple_New((Py_ssize_t) npar);
		for(i=0; i<npar; i++) {
			PyList_SET_ITEM(lower,i,Py_BuildValue("d",self->bifparams->GetIthParameterLowerBound(i)));
			PyList_SET_ITEM(upper,i,Py_BuildValue("d",self->bifparams->GetIthParameterUpperBound(i)));
			PyTuple_SET_ITEM(steps,i,Py_BuildValue("i",self->bifparams->GetNumberOfSteps(i)));
		}
		result = PyDict_New();
		PyDict_SetItemString(result,"lower",lower);
		PyDict_SetItemString(result,"upper",upper);
		PyDict_SetItemString(result,"steps",steps);
	}
	else {
		result = PyObject_GenericGetAttr((PyObject*)self, name);
	}
	Py_DECREF(name);
	return result;
}

static int pyBalParameters_setattro(pyBalParameters *self, PyObject *name, PyObject *value) {
	int err = 0;
	Py_INCREF(name);
	char *n = PyString_AsString(name);

	//printf("pyBalParameters_setattro name = %s\n", n);
	if (strcmp(n, "npars") == 0) {
		int npar;
		if(PyArg_Parse(value,"i",&npar) && npar >= 0) {
			self->bifparams->SetNumber(npar);
		}
		else {
			PyErr_SetString(PyExc_ValueError,"The number of parameters must be non-negative.");
			err = -1;
		}
	}
	else if (strcmp(n, "steps") == 0) {
		PyObject* steps;
		if(!PyArg_Parse(value,"O",&steps)) {
			PyErr_SetString(PyExc_ValueError,"ERROR.");
			err = -1;
		}
		else {
			Py_ssize_t (*objsize)(PyObject *);
			PyObject* (*getitem)(PyObject *, Py_ssize_t);
			if(PyTuple_Check(steps)) {
				objsize = PyTuple_Size;
				getitem = PyTuple_GetItem;
			}
			else if(PyList_Check(steps)) {
				objsize = PyList_Size;
				getitem = PyList_GetItem;
			}
			else {
				PyErr_SetString(PyExc_ValueError,"The argument must be either a tuple or a list.");
				err = -1;
			}
			if(! err) {
				if(objsize(steps) == self->bifparams->GetNumber()) {
					int *s = new int[self->bifparams->GetNumber()];
					for(int i=0; i<self->bifparams->GetNumber(); i++) {
						s[i] = PyInt_AsLong(getitem(steps,i));
						if(s[i] <= 0) {
							PyErr_SetString(PyExc_ValueError,"The number of steps must be positive.");
							err = -1;
							break;
						}
					}
					if(err != -1) {
						for(int i=0; i<self->bifparams->GetNumber(); i++)
							self->bifparams->SetNumberOfSteps(i,s[i]);
					}
					delete s;
				}
				else {
					PyErr_SetString(PyExc_ValueError,"The number of parameters provided does not match the number of parameters in the class.");
					err = -1;
				}
			}
		}
	}
	else {
		err = PyObject_GenericSetAttr((PyObject*)self, name, value);
	}
	Py_DECREF(name);
	return err;
}

static PyObject * pyBalParameters_iternext(pyBalParameters *self) {
	static bool first = true;
	if(self->bifparams->HasNext()) {
		if(first)
			first = false;
		else
			self->bifparams->Next();
		int npar = self->bifparams->GetNumber();
		PyObject * result = PyList_New((Py_ssize_t) npar);
		for(int i=0; i<npar; i++)
			PyList_SET_ITEM(result,i,Py_BuildValue("d",self->bifparams->At(i)));
		//self->bifparams->Next();
		return result;
	}
	self->bifparams->Reset();
	return NULL;
}

static PyObject * pyBalParameters_iter(pyBalParameters *self) {
	Py_INCREF(self);
	return (PyObject *) self;
}

static PyObject * pyBalParameters_par(PyObject *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"index","value",NULL};
	int index, retval = -1;
	double value;
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "id", kwlist, &index, &value))
		return Py_BuildValue("i",retval);
	Py_INCREF(self);
	pyBalParameters * bp = (pyBalParameters *) self;	
	if (bp->bifparams->SetIthParameter(index,value))
		retval = 0;
	Py_DECREF(self);
	return Py_BuildValue("i",retval);
}

static PyObject * pyBalParameters_bifpar(PyObject *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"index","values",NULL};
	int index, retval = 0;
	PyObject *value;
	
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "iO", kwlist, &index, &value))
		return Py_BuildValue("i",-1);
	
	Py_ssize_t (*objsize)(PyObject *) = NULL;
	PyObject* (*getitem)(PyObject *, Py_ssize_t) = NULL;
	if(PyTuple_Check(value)) {
		objsize = PyTuple_Size;
		getitem = PyTuple_GetItem;
	}
	else if(PyList_Check(value)) {
		objsize = PyList_Size;
		getitem = PyList_GetItem;
	}
	else {
		PyErr_SetString(PyExc_ValueError,"The argument must be either a tuple or a list.");
		retval = -1;
	}
	
	int sz = objsize(value);
	if(retval != -1 && (sz == 2 || sz == 3)) {
		Py_INCREF(self);
		pyBalParameters * bp = (pyBalParameters *) self;
		if(bp->bifparams->SetIthParameterLowerBound(index,PyFloat_AsDouble(getitem(value,0))) &&
		   bp->bifparams->SetIthParameterUpperBound(index,PyFloat_AsDouble(getitem(value,1))))
			retval = 0;
		if(sz == 3)
			bp->bifparams->SetNumberOfSteps(index,PyInt_AsLong(getitem(value,2)));
		Py_DECREF(self);
	}
	return Py_BuildValue("i",retval);
}

static PyObject * pyBalParameters_setpars(PyObject *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"value","index","type", NULL};
	PyObject *index = NULL, *value = NULL;
	char *type = NULL;
	
	int i, npar, retval = 0;
	int *idx = NULL;
	double *p = NULL;
	
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "|OOs", kwlist, &value, &index, &type))
		return Py_BuildValue("i",-1);

	if (value == NULL)
		return Py_BuildValue("i",-1);
	
	// the first argument must be a list containing the parameters
	if(!(PyList_Check(value) || PyFloat_Check(value) || PyInt_Check(value))) {
		PyErr_SetString(PyExc_TypeError,"The parameters must be contained in a list.");
		return Py_BuildValue("i",-2);
	}
	
	Py_INCREF(self);
	pyBalParameters * parameters = (pyBalParameters *) self;
	
	// the second argument may be either a list or a tuple containing the indices of the parameters
	if(index != NULL) {
		if(PyTuple_Check(index)) {
			npar = PyTuple_Size(index);
			idx = new int[npar];
			for(i=0; i<npar; i++)
				idx[i] = PyInt_AsLong(PyTuple_GetItem(index,i));
		}
		else if(PyList_Check(index)) {
			npar = PyList_Size(index);
			idx = new int[npar];
			for(i=0; i<npar; i++)
				idx[i] = PyInt_AsLong(PyList_GetItem(index,i));
		}
		else if(PyInt_Check(index)) {
			npar = 1;
			idx = new int[npar];
			idx[0] = PyInt_AsLong(index);
		}
		else {
			PyErr_SetString(PyExc_TypeError,"The indices must be contained in a tuple or in a list.");
			Py_DECREF(self);
			return Py_BuildValue("i",-3);
		}
	}
	else {
		if(PyList_Check(value))
			npar = PyList_Size(value);
		else
			npar = 1;
		if(npar != parameters->bifparams->GetNumber()) {
			PyErr_SetString(PyExc_ValueError,"Wrong number of parameters.");
			Py_DECREF(self);
			return Py_BuildValue("i",-4);			
		}
		idx = new int[npar];
		for(i=0; i<npar; i++)
			idx[i] = i;
	}
	
	if((PyList_Check(value) && npar != PyList_Size(value)) || npar > parameters->bifparams->GetNumber()) {
		PyErr_SetString(PyExc_ValueError,"Wrong number of parameters.");
		Py_DECREF(self);
		delete idx;
		return Py_BuildValue("i",-5);
	}

	p = new double[npar];
	if(PyList_Check(value)) {
		for(i=0; i<npar; i++)
			p[i] = PyFloat_AsDouble(PyList_GetItem(value,i));
	}
	else {
		if(PyFloat_Check(value))
			p[0] = PyFloat_AsDouble(value);
		else
			p[0] = (double) PyInt_AsLong(value);
	}
	
	if(type == NULL) {
		for(i=0; i<npar; i++)
			parameters->bifparams->SetIthParameter(idx[i],p[i]);
	}
	else {
		if(strcmp(type,"lower") == 0) {
			for(i=0; i<npar; i++)
				parameters->bifparams->SetIthParameterLowerBound(idx[i],p[i]);
		}
		else if(strcmp(type,"upper") == 0) {
			for(i=0; i<npar; i++)
				parameters->bifparams->SetIthParameterUpperBound(idx[i],p[i]);
		}
		else {
			PyErr_SetString(PyExc_ValueError,"Wrong type of parameters (either 'upper' or 'lower' is accepted).");
			retval = -6;
		}
	}
	
	if(retval == 0)
		parameters->bifparams->Reset();
	
	delete p;
	delete idx;
	Py_DECREF(self);
	return Py_BuildValue("i",retval);
}

static PyObject * pyBalParameters_hastuples(pyBalParameters *self) {
	return Py_BuildValue("b",self->bifparams->HasTuples());
}

static PyObject * pyBalParameters_hasnext(pyBalParameters *self) {
	return Py_BuildValue("b",self->bifparams->HasNext());
}

static PyObject * pyBalParameters_next(pyBalParameters *self) {
	self->bifparams->Next();
	return Py_None;
}

static PyObject * pyBalParameters_reset(pyBalParameters *self) {
	self->bifparams->Reset();
	return Py_None;
}

static PyMethodDef pyBalParameters_methods[] = {
	{"par", (PyCFunction) pyBalParameters_par, METH_VARARGS | METH_KEYWORDS, "Set the value of the i-th parameter"},
	{"bifpar", (PyCFunction) pyBalParameters_bifpar, METH_VARARGS | METH_KEYWORDS,
		"Set the upper and lower bounds of the i-th parameter (and optionally the number of steps)"},
	{"setpars", (PyCFunction) pyBalParameters_setpars, METH_VARARGS | METH_KEYWORDS, "Set the values of the parameters"},
	{"hastuples", (PyCFunction) pyBalParameters_hastuples, METH_NOARGS, "Check whether the parameters' tuples are not finished"},
	{"hasnext", (PyCFunction) pyBalParameters_hasnext, METH_NOARGS, "Check whether there's another tuple of parameters"},
	{"next", (PyCFunction) pyBalParameters_next, METH_NOARGS, "Increment the current tuple of parameters"},
	{"reset", (PyCFunction) pyBalParameters_reset, METH_NOARGS, "Reset the current tuple of parameters"},
	{NULL}	/* Sentinel */
};

static PyTypeObject pyBalParametersType = {
	PyObject_HEAD_INIT(NULL)
	0,									/* ob_size */
	"bal.Parameters",					/* ob_name */
	sizeof(pyBalParameters),			/* tp_basicsize */
	0,									/* tp_itemsize */
	(destructor)pyBalParameters_dealloc,/* tp_dealloc */
	0,									/* tp_print */
	0,									/* tp_getattr */
	0,									/* tp_setattr */
	0,									/* tp_compare */
	0,									/* tp_repr */
	0,									/* tp_as_number */
	0,									/* tp_as_sequence */
	0,									/* tp_as_mapping */
	0,									/* tp_hash */
	0,									/* tp_call */
	0,									/* tp_str */
	(getattrofunc)pyBalParameters_getattro,	/* tp_getattro  */
	(setattrofunc)pyBalParameters_setattro,	/* tp_setattro  */
	0,									/* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/* tp_flags */
	"Parameters object",				/* tp_doc */
	0,									/* tp_traverse */
    0,									/* tp_clear */
    0,									/* tp_richcompare */
    0,									/* tp_weaklistoffset */
    (getiterfunc)pyBalParameters_iter,	/* tp_iter */
    (iternextfunc)pyBalParameters_iternext,	/* tp_iternext */
    pyBalParameters_methods,			/* tp_methods */
    0,									/* tp_members */
    0,									/* tp_getset */
    0,									/* tp_base */
    0,									/* tp_dict */
    0,									/* tp_descr_get */
    0,									/* tp_descr_set */
    0,									/* tp_dictoffset */
    (initproc)pyBalParameters_init,		/* tp_init */
    0,									/* tp_alloc */
    pyBalParameters_new					/* tp_new */
};

/******************** pyBalODESolver ********************/

static PyObject * pyBalODESolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	pyBalODESolver *self;
	self = (pyBalODESolver *) type->tp_alloc(type,0);
	self->solver = NULL;
	self->dynsys = NULL;
	self->params = NULL;
	return (PyObject *) self;
}

static int pyBalODESolver_init(pyBalODESolver *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"system","parameters",NULL};
	PyObject *ds = NULL;
	PyObject *p = NULL;
	
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist, &ds, &p))
		return NULL;
	
	self->params = (pyBalParameters *) p;
	self->dynsys = (pyBalDynamicalSystem *) ds;
	if(self->params->bifparams->GetNumber() != self->dynsys->dynsys->GetNumberOfParameters()) {
		self->params = NULL;
		self->dynsys = NULL;
		return NULL;
	}
	self->dynsys->dynsys->SetParameters(self->params->bifparams);
	self->solver = bal::ODESolver::Create();
	self->solver->SetDynamicalSystem(self->dynsys->dynsys);
	self->solver->SetTransientDuration(0.0);
	self->solver->SetFinalTime(10.0);
	self->solver->SetIntegrationMode(bal::TRAJ);
	
	return 0;
}

static void pyBalODESolver_dealloc(pyBalODESolver *self) {
	if(self->solver != NULL)
		self->solver->Destroy();
	self->ob_type->tp_free((PyObject *) self);
}

static PyObject * pyBalODESolver_getattro(pyBalODESolver *self, PyObject *name) {
	Py_INCREF(name);
	char *n = PyString_AsString(name);
	PyObject *result = NULL;
	
	if (strcmp(n, "system") == 0) {
		Py_INCREF(self->dynsys);
		result = (PyObject *) self->dynsys;
	}
	else if (strcmp(n, "parameters") == 0) {
		Py_INCREF(self->params);
		result = (PyObject *) self->params;
	}
	else if(strcmp(n, "dt") == 0) {
		result = Py_BuildValue("d",self->solver->GetTimeStep());
	}
	else if(strcmp(n, "lyap_dt") == 0) {
		result = Py_BuildValue("d",self->solver->GetLyapunovTimeStep());
	}
	else if(strcmp(n, "t0") == 0) {
		result = Py_BuildValue("d",self->solver->GetInitialTime());
	}
	else if(strcmp(n, "tstop") == 0) {
		result = Py_BuildValue("d",self->solver->GetFinalTime());
	}
	else if(strcmp(n, "ttran") == 0) {
		result = Py_BuildValue("d",self->solver->GetTransientDuration());
	}
	else if(strcmp(n, "x0") == 0) {
		int i, ndim = self->dynsys->dynsys->GetDimension();
		N_Vector x0 = self->solver->GetX0();
		result = PyList_New(ndim);
		for(i=0; i<ndim; i++)
			PyList_SET_ITEM(result,i,Py_BuildValue("d",NV_Ith_S(x0,i)));
	}
	else if(strcmp(n, "mode") == 0) {
		switch(self->solver->GetIntegrationMode()) {
			case bal::TRAJ:
				result = Py_BuildValue("s","trajectory");
				break;
			case bal::EVENTS:
				result = Py_BuildValue("s","events");
				break;
			case bal::BOTH:
				result = Py_BuildValue("s","trajectory + events");
				break;
		}
	}
	else if(strcmp(n, "intersections") == 0) {
		result = Py_BuildValue("i",self->solver->GetMaxNumberOfIntersections());
	}
	else if(strcmp(n, "nturns") == 0) {
		result = Py_BuildValue("i",self->solver->GetNumberOfTurns());
	}
	else {
		result = PyObject_GenericGetAttr((PyObject*)self, name);
	}
	Py_DECREF(name);
	return result;
}

static int pyBalODESolver_setattro(pyBalODESolver *self, PyObject *name, PyObject *value) {
	int err = 0;
	Py_INCREF(name);
	char *n = PyString_AsString(name);

	//printf("pyBalODESolver_settatro %s\n", n);
	if (strcmp(n, "system") == 0) {
		Py_DECREF(self->dynsys);
		Py_INCREF(value);
		self->dynsys = (pyBalDynamicalSystem *) value;
		self->solver->SetDynamicalSystem(self->dynsys->dynsys);
	}
	else if (strcmp(n, "parameters") == 0) {
		Py_DECREF(self->params);
		Py_INCREF(value);
		self->params = (pyBalParameters *) value;
		self->dynsys->dynsys->SetParameters(self->params->bifparams);
	}
	else if(strcmp(n, "dt") == 0) {
		double dt = PyFloat_AsDouble(value);
		if(dt > 0)
			self->solver->SetTimeStep(dt);
		else
			err = -1;
	}
	else if(strcmp(n, "lyap_dt") == 0) {
		double dt = PyFloat_AsDouble(value);
		if(dt > 0)
			self->solver->SetLyapunovTimeStep(dt);
		else
			err = -1;
	}
	else if(strcmp(n, "t0") == 0) {
		double t0 = PyFloat_AsDouble(value);
		if(t0 >= 0)
			self->solver->SetInitialTime(t0);
		else
			err = -1;
	}
	else if(strcmp(n, "tstop") == 0) {
		double tstop = PyFloat_AsDouble(value);
		if(tstop >= 0)
			self->solver->SetFinalTime(tstop);
		else
			err = -1;
	}
	else if(strcmp(n, "ttran") == 0) {
		double ttran = PyFloat_AsDouble(value);
		if(ttran >= 0)
			self->solver->SetTransientDuration(ttran);
		else
			err = -1;
	}
	else if(strcmp(n, "x0") == 0) {
		int i, ndim = self->dynsys->dynsys->GetDimension();
		double *x0 = new double[ndim];
		for(i=0; i<ndim; i++)
			x0[i] = PyFloat_AsDouble(PyList_GET_ITEM(value,i));
		self->solver->SetX0(x0);
		delete x0;
	}
	else if(strcmp(n, "mode") == 0) {
		Py_INCREF(value);
		char *v = PyString_AsString(value);
		if(strcmp(v, "trajectory") == 0)
			self->solver->SetIntegrationMode(bal::TRAJ);
		else if(strcmp(v, "events") == 0)
			self->solver->SetIntegrationMode(bal::EVENTS);
		else if(strcmp(v, "trajectory + events") == 0 || strcmp(v, "both") == 0)
			self->solver->SetIntegrationMode(bal::BOTH);
		else if(strcmp(v, "lyap") == 0 || strcmp(v, "lyapunov"))
			self->solver->SetIntegrationMode(bal::LYAP);
		else
			err = -1;
		Py_DECREF(value);
	}	
	else if(strcmp(n, "intersections") == 0) {
		self->solver->SetMaxNumberOfIntersections(PyInt_AsLong(value));
	}
	else {
		err = -1;
	}
	Py_DECREF(name);
	return 0;
}

static PyObject * pyBalODESolver_solve(pyBalODESolver *self) {
	bool flag = self->solver->Solve();
	return Py_BuildValue("i",flag ? 0 : -1);
}

static PyObject * pyBalODESolver_getsolution(pyBalODESolver * self) {
	bal::Solution *sol = self->solver->GetSolution();
	if(sol == NULL)
		return Py_BuildValue("i",-1);
	pyBalSolution *pbs = createPyBalSolution(sol);
	sol->Destroy();
	return (PyObject *) pbs;
}

static PyObject * pyBalODESolver_writeorbit(pyBalODESolver *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"filename",NULL};
	char *filename = NULL;
	
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "|s", kwlist, &filename))
		return NULL;
	
	return Py_BuildValue("i", self->solver->SaveOrbit(filename) ? 0 : -1);
}

static PyMethodDef pyBalODESolver_methods[] = {
	{"run", (PyCFunction) pyBalODESolver_solve, METH_NOARGS, "Integrate the dynamical system"},
	{"solution", (PyCFunction) pyBalODESolver_getsolution, METH_NOARGS, "Get the solution of the integration"},
	{"write_orbit", (PyCFunction) pyBalODESolver_writeorbit, METH_VARARGS | METH_KEYWORDS, "Write a closed orbit to file"},
	{NULL}	/* Sentinel */
};

static PyTypeObject pyBalODESolverType = {
	PyObject_HEAD_INIT(NULL)
	0,									/* ob_size */
	"bal.ODESolver",					/* ob_name */
	sizeof(pyBalODESolver),				/* tp_basicsize */
	0,									/* tp_itemsize */
	(destructor)pyBalODESolver_dealloc,	/* tp_dealloc */
	0,									/* tp_print */
	0,									/* tp_getattr */
	0,									/* tp_setattr */
	0,									/* tp_compare */
	0,									/* tp_repr */
	0,									/* tp_as_number */
	0,									/* tp_as_sequence */
	0,									/* tp_as_mapping */
	0,									/* tp_hash */
	0,									/* tp_call */
	0,									/* tp_str */
	(getattrofunc)pyBalODESolver_getattro,/* tp_getattro  */
	(setattrofunc)pyBalODESolver_setattro,/* tp_setattro  */
	0,									/* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/* tp_flags */
	"ODE solver object",				/* tp_doc */
	0,									/* tp_traverse */
    0,									/* tp_clear */
    0,									/* tp_richcompare */
    0,									/* tp_weaklistoffset */
    0,									/* tp_iter */
    0,									/* tp_iternext */
    pyBalODESolver_methods,				/* tp_methods */
    0,									/* tp_members */
    0,									/* tp_getset */
    0,									/* tp_base */
    0,									/* tp_dict */
    0,									/* tp_descr_get */
    0,									/* tp_descr_set */
    0,									/* tp_dictoffset */
    (initproc)pyBalODESolver_init,		/* tp_init */
    0,									/* tp_alloc */
    pyBalODESolver_new					/* tp_new */
};


/******************** pyBalSolution ********************/


static PyMethodDef pyBalSolution_methods[] = {
	{NULL}	/* Sentinel */
};

static PyTypeObject pyBalSolutionType = {
	PyObject_HEAD_INIT(NULL)
	0,									/* ob_size */
	"bal.Solution",						/* ob_name */
	sizeof(pyBalSolution),				/* tp_basicsize */
	0,									/* tp_itemsize */
	(destructor)pyBalSolution_dealloc,	/* tp_dealloc */
	0,									/* tp_print */
	0,									/* tp_getattr */
	0,									/* tp_setattr */
	0,									/* tp_compare */
	0,									/* tp_repr */
	0,									/* tp_as_number */
	0,									/* tp_as_sequence */
	0,									/* tp_as_mapping */
	0,									/* tp_hash */
	0,									/* tp_call */
	0,									/* tp_str */
	(getattrofunc)pyBalSolution_getattro,/* tp_getattro  */
	0,									/* tp_setattro  */
	0,									/* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/* tp_flags */
	"Solution object",					/* tp_doc */
	0,									/* tp_traverse */
    0,									/* tp_clear */
    0,									/* tp_richcompare */
    0,									/* tp_weaklistoffset */
    0,									/* tp_iter */
    0,									/* tp_iternext */
    pyBalSolution_methods,				/* tp_methods */
    0,									/* tp_members */
    0,									/* tp_getset */
    0,									/* tp_base */
    0,									/* tp_dict */
    0,									/* tp_descr_get */
    0,									/* tp_descr_set */
    0,									/* tp_dictoffset */
    0,									/* tp_init */
    0,									/* tp_alloc */
    0									/* tp_new */
};

static pyBalSolution* createPyBalSolution(bal::Solution *s) {
	pyBalSolution *solution = (pyBalSolution *) pyBalSolutionType.tp_alloc(&pyBalSolutionType,0);
	
	// Parameters
	solution->params = bal::Parameters::Copy(s->GetParameters());
	
	// Data
	int r,c,i,j;
	s->GetSize(&r,&c);
	double *buffer = s->GetData();
	PyObject *t = PyList_New((Py_ssize_t) r);
	PyObject *x = PyList_New((Py_ssize_t) r*(c-2));
	PyObject *labels = PyList_New((Py_ssize_t) r);
	
	for(i=0; i<r; i++) {
		PyList_SET_ITEM(t,i,Py_BuildValue("d",buffer[i*c]));
		for(j=1; j<c-1; j++)
			PyList_SET_ITEM(x,i*(c-2)+(j-1),Py_BuildValue("d",buffer[i*c+j]));
		PyList_SET_ITEM(labels,i,Py_BuildValue("i",buffer[(i+1)*c-1]));
	}
	solution->data = PyDict_New();
	//Py_INCREF(t);
	PyDict_SetItemString(solution->data,"t",t);
	//Py_INCREF(x);
	PyDict_SetItemString(solution->data,"x",x);
	//Py_INCREF(labels);
	PyDict_SetItemString(solution->data,"labels",labels);
	
	return solution;
}

static void pyBalSolution_dealloc(pyBalSolution *self) {
	self->params->Destroy();
	//PyTypeList.tp_dealloc(PyDict_GetItemString(self->data,"t"));
	//printf("refcnt before: %d\n", self->data->ob_refcnt);
	Py_DECREF(self->data);
	//printf("refcnt after: %d\n", self->data->ob_refcnt);
	self->ob_type->tp_free((PyObject *) self);
}

static PyObject * pyBalSolution_getattro(pyBalSolution *self, PyObject *name) {
	Py_INCREF(name);
	char *n = PyString_AsString(name);
	PyObject *result = NULL;
	
	if (strcmp(n, "parameters") == 0) {
		int i, npar = self->params->GetNumber();
		result = (PyObject *) PyList_New((Py_ssize_t) npar);
		for(i=0; i<npar; i++)
			PyList_SET_ITEM(result,i,Py_BuildValue("d",self->params->At(i)));
	}
	else if(strcmp(n, "data") == 0) {
		Py_INCREF(self->data);
		result = self->data;
	}
	else {
		result = PyObject_GenericGetAttr((PyObject*)self, name);
	}
	Py_DECREF(name);
	return result;
}

/******************** pyBalBifurcationDiagram ********************/


static PyObject * pyBalBifurcationDiagram_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	pyBalBifurcationDiagram *self;
	self = (pyBalBifurcationDiagram *) type->tp_alloc(type,0);
	self->diagram = NULL;
	self->params = NULL;
	self->dynsys = NULL;
	self->ic = NULL;
	self->nic = -1;
	return (PyObject *) self;
}

static int pyBalBifurcationDiagram_init(pyBalBifurcationDiagram *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"system","parameters",NULL};
	PyObject *ds = NULL;
	PyObject *p = NULL;
	
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist, &ds, &p))
		return NULL;
	
	self->params = (pyBalParameters *) p;
	
	self->dynsys = (pyBalDynamicalSystem *) ds;
	if(self->params->bifparams->GetNumber() != self->dynsys->dynsys->GetNumberOfParameters()) {
		self->params = NULL;
		self->dynsys = NULL;
		return NULL;
	}
	self->dynsys->dynsys->SetParameters(self->params->bifparams);
	double *x0 = new double[self->dynsys->dynsys->GetDimension()];
	for(int i=0; i<self->dynsys->dynsys->GetDimension(); i++)
		x0[i] = 0.0;
	
	self->diagram = bal::BifurcationDiagram::Create();
	self->diagram->SetDynamicalSystem(self->dynsys->dynsys);
	self->diagram->GetODESolver()->SetIntegrationMode(bal::EVENTS);
	self->diagram->GetODESolver()->HaltAtEquilibrium(true);
	self->diagram->GetODESolver()->HaltAtCycle(false);
	self->diagram->GetODESolver()->SetTransientDuration(1e3);
	self->diagram->GetODESolver()->SetFinalTime(1e4);
	self->diagram->GetODESolver()->SetMaxNumberOfIntersections(200);
	self->diagram->GetODESolver()->SetX0(x0);
	self->diagram->RestartFromX0(true);
	self->diagram->SetNumberOfThreads(2);
	self->diagram->SetMode(bal::PARAMS);
	
	delete x0;
	return 0;
}

static void pyBalBifurcationDiagram_dealloc(pyBalBifurcationDiagram *self) {
	if(self->diagram != NULL)
		self->diagram->Destroy();
	if(self->nic > 0) {
		for(int i=0; i<self->nic; i++)
			delete [] self->ic[i];
		delete self->ic;
	}
	self->ob_type->tp_free((PyObject *) self);
}

static PyObject * pyBalBifurcationDiagram_getattro(pyBalBifurcationDiagram *self, PyObject *name) {
	Py_INCREF(name);
	char *n = PyString_AsString(name);
	PyObject *result = NULL;
		
	//printf("pyBalBifurcationDiagram_gettatro %s\n", n);
	if (strcmp(n, "system") == 0) {
		Py_INCREF(self->dynsys);
		result = (PyObject *) self->dynsys;
	}
	else if (strcmp(n, "parameters") == 0) {
		Py_INCREF(self->params);
		result = (PyObject *) self->params;
	}
	else if(strcmp(n, "dt") == 0) {
		result = Py_BuildValue("d",self->diagram->GetODESolver()->GetTimeStep());
	}
	else if(strcmp(n, "lyap_dt") == 0) {
		result = Py_BuildValue("d",self->diagram->GetODESolver()->GetLyapunovTimeStep());
	}
	else if(strcmp(n, "t0") == 0) {
		result = Py_BuildValue("d",self->diagram->GetODESolver()->GetInitialTime());
	}
	else if(strcmp(n, "tstop") == 0) {
		result = Py_BuildValue("d",self->diagram->GetODESolver()->GetFinalTime());
	}
	else if(strcmp(n, "ttran") == 0) {
		result = Py_BuildValue("d",self->diagram->GetODESolver()->GetTransientDuration());
	}
	else if(strcmp(n, "x0") == 0) {
		int i, ndim = self->dynsys->dynsys->GetDimension();
		N_Vector x0 = self->diagram->GetODESolver()->GetX0();
		result = PyList_New(ndim);
		for(i=0; i<ndim; i++)
			PyList_SET_ITEM(result,i,Py_BuildValue("d",NV_Ith_S(x0,i)));
	}
	else if(strcmp(n, "mode") == 0) {
		switch(self->diagram->GetODESolver()->GetIntegrationMode()) {
			case bal::TRAJ:
				result = Py_BuildValue("s","trajectory");
				break;
			case bal::EVENTS:
				result = Py_BuildValue("s","events");
				break;
			case bal::BOTH:
				result = Py_BuildValue("s","trajectory + events");
				break;
			case bal::LYAP:
				result = Py_BuildValue("s","lyapunov exponents");
				break;
		}
	}
	else if(strcmp(n, "diagram_mode") == 0) {
		switch(self->diagram->GetMode()) {
			case bal::PARAMS:
				result = Py_BuildValue("s","parameters");
				break;
			case bal::IC:
				result = Py_BuildValue("s","IC");
				break;
		}
	}
	else if(strcmp(n, "intersections") == 0) {
		result = Py_BuildValue("i",self->diagram->GetODESolver()->GetMaxNumberOfIntersections());
	}
	else if(strcmp(n, "outfile") == 0) {
		result = Py_BuildValue("s",self->diagram->GetFilename());
	}
	else if(strcmp(n, "resetX0") == 0) {
		result = (self->diagram->RestartsFromX0() ? Py_True : Py_False);
	}
	else if(strcmp(n, "equilbreak") == 0) {
		result = (self->diagram->GetODESolver()->HaltsAtEquilibrium() ? Py_True : Py_False);
	}
	else if(strcmp(n, "cyclebreak") == 0) {
		result = (self->diagram->GetODESolver()->HaltsAtCycle() ? Py_True : Py_False);
	}
	else if(strcmp(n, "nthreads") == 0) {
		result = Py_BuildValue("i",self->diagram->GetNumberOfThreads());
	}
	else {
		result = PyObject_GenericGetAttr((PyObject*)self, name);
	}
	Py_DECREF(name);
	return result;
}

static int pyBalBifurcationDiagram_setattro(pyBalBifurcationDiagram *self, PyObject *name, PyObject *value) {
	int err = 0;
	Py_INCREF(name);
	char *n = PyString_AsString(name);
	
	//printf("pyBalBifurcationDiagram_settatro %s\n", n);
	if (strcmp(n, "system") == 0) {
		Py_DECREF(self->dynsys);
		Py_INCREF(value);
		self->dynsys = (pyBalDynamicalSystem *) value;
		self->diagram->SetDynamicalSystem(self->dynsys->dynsys);
	}
	else if (strcmp(n, "parameters") == 0) {
		Py_DECREF(self->params);
		Py_INCREF(value);
		self->params = (pyBalParameters *) value;
		self->dynsys->dynsys->SetParameters(self->params->bifparams);
	}
	else if(strcmp(n, "dt") == 0) {
		double dt = PyFloat_AsDouble(value);
		if(dt > 0)
			self->diagram->GetODESolver()->SetTimeStep(dt);
		else
			err = -1;
	}
	else if(strcmp(n, "lyap_dt") == 0) {
		double dt = PyFloat_AsDouble(value);
		if(dt > 0)
			self->diagram->GetODESolver()->SetLyapunovTimeStep(dt);
		else
			err = -1;
	}
	else if(strcmp(n, "t0") == 0) {
		double t0 = PyFloat_AsDouble(value);
		if(t0 >= 0)
			self->diagram->GetODESolver()->SetInitialTime(t0);
		else
			err = -1;
	}
	else if(strcmp(n, "tstop") == 0) {
		double tstop = PyFloat_AsDouble(value);
		if(tstop >= 0)
			self->diagram->GetODESolver()->SetFinalTime(tstop);
		else
			err = -1;
	}
	else if(strcmp(n, "ttran") == 0) {
		double ttran = PyFloat_AsDouble(value);
		if(ttran >= 0)
			self->diagram->GetODESolver()->SetTransientDuration(ttran);
		else
			err = -1;
	}
	else if(strcmp(n, "x0") == 0) {
		if(self->diagram->GetMode() == bal::PARAMS) {
			int i, ndim = self->dynsys->dynsys->GetDimension();
			double *x0 = new double[ndim];
			for(i=0; i<ndim; i++)
				x0[i] = PyFloat_AsDouble(PyList_GET_ITEM(value,i));
			self->diagram->GetODESolver()->SetX0(x0);
			delete x0;
		}
		else if(self->diagram->GetMode() == bal::IC) {
			int i, j, ndim = self->dynsys->dynsys->GetDimension();
			PyObject *list;
			if(self->nic > 0) {
				for(i=0; i<self->nic; i++)
					delete [] self->ic[i];
				delete self->ic;
			}
			self->nic = PyList_Size(value);
			self->ic = new double*[self->nic];
			for(i=0; i<self->nic; i++) {
				self->ic[i] = new double[ndim];
				list = PyList_GET_ITEM(value,i);
				for(j=0; j<ndim; j++)
					self->ic[i][j] = PyFloat_AsDouble(PyList_GET_ITEM(list,j));
			}
			self->diagram->SetInitialConditions(self->nic, self->ic);
		}
	}
	else if(strcmp(n, "mode") == 0) {
		Py_INCREF(value);
		char *v = PyString_AsString(value);
		if(strcmp(v, "trajectory") == 0)
			self->diagram->GetODESolver()->SetIntegrationMode(bal::TRAJ);
		else if(strcmp(v, "events") == 0)
			self->diagram->GetODESolver()->SetIntegrationMode(bal::EVENTS);
		else if(strcmp(v, "trajectory + events") == 0 || strcmp(v, "both") == 0)
			self->diagram->GetODESolver()->SetIntegrationMode(bal::BOTH);
		else if(strcmp(v, "lyap") == 0 || strcmp(v, "lyapunov") == 0)
			self->diagram->GetODESolver()->SetIntegrationMode(bal::LYAP);
		else
			err = -1;
		Py_DECREF(value);
	}
	else if(strcmp(n, "diagram_mode") == 0) {
		Py_INCREF(value);
		char *v = PyString_AsString(value);
		if(strcmp(v, "parameters") == 0)
			self->diagram->SetMode(bal::PARAMS);
		else if(strcmp(v, "IC") == 0 || strcmp(v, "ic") == 0)
			self->diagram->SetMode(bal::IC);
		else
			err = -1;
		Py_DECREF(value);
	}	
	else if(strcmp(n, "intersections") == 0) {
		self->diagram->GetODESolver()->SetMaxNumberOfIntersections(PyInt_AsLong(value));
	}
	else if(strcmp(n, "outfile") == 0) {
		Py_INCREF(value);
		self->diagram->SetFilename(PyString_AsString(value));
	}
	else if(strcmp(n, "resetX0") == 0) {
		if(PyBool_Check(value))
			self->diagram->RestartFromX0(value == Py_False ? false : true);
		else if(PyInt_Check(value))
			self->diagram->RestartFromX0(PyInt_AsLong(value) == 0 ? false : true);
		else if(PyString_Check(value)) {
			Py_INCREF(value);
			char *n = PyString_AsString(value);
			if(strcmp(n,"yes") == 0)
				self->diagram->RestartFromX0(true);
			else if(strcmp(n,"no") == 0)
				self->diagram->RestartFromX0(false);
			else
				err = -1;
			Py_DECREF(value);
		}
		else
			err = -1;
	}
	else if(strcmp(n, "nthreads") == 0) {
		self->diagram->SetNumberOfThreads(PyInt_AsLong(value));
	}
	else if(strcmp(n, "equilbreak") == 0) {
		if(PyBool_Check(value))
			self->diagram->GetODESolver()->HaltAtEquilibrium(value == Py_False ? false : true);
		else if(PyInt_Check(value))
			self->diagram->GetODESolver()->HaltAtEquilibrium(PyInt_AsLong(value) == 0 ? false : true);
		else
			err = -1;
	}
	else if(strcmp(n, "cyclebreak") == 0) {
		if(PyBool_Check(value))
			self->diagram->GetODESolver()->HaltAtCycle(value == Py_False ? false : true);
		else if(PyInt_Check(value))
			self->diagram->GetODESolver()->HaltAtCycle(PyInt_AsLong(value) == 0 ? false : true);
		else
			err = -1;
	}
	else if(strcmp(n, "equiltol") == 0) {
		self->diagram->GetODESolver()->SetEquilibriumTolerance(PyFloat_AsDouble(value));
	}
	else if(strcmp(n, "reltol") == 0) {
		self->diagram->GetODESolver()->SetRelativeTolerance(PyFloat_AsDouble(value));
	}
	else if(strcmp(n, "abstol") == 0) {
		self->diagram->GetODESolver()->SetAbsoluteTolerance(PyFloat_AsDouble(value));
	}
	else {
		err = -1;
	}
	Py_DECREF(name);
	return 0;
}

static PyObject * pyBalBifurcationDiagram_compute(pyBalBifurcationDiagram *self) {
	self->diagram->ComputeDiagram();
	return Py_BuildValue("i",0);
}

static PyObject * pyBalBifurcationDiagram_summary(pyBalBifurcationDiagram *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"filename",NULL};
	char *filename = NULL;
	
	if (! PyArg_ParseTupleAndKeywords(args, kwds, "|s", kwlist, &filename))
		return NULL;

	PyObject *result;

	if(filename != NULL) {
		result = (self->diagram->SaveSummaryData(filename) ? Py_True : Py_False);
	}
	else {
		int size[2];
		double **data = self->diagram->GetSummaryData(size);
		if(data == NULL) {
			result = Py_None;
		}
		else {
			int i, j;
			PyObject *row;
			result = PyList_New((Py_ssize_t) size[0]);
			for(i=0; i<size[0]; i++) {
				row = PyList_New((Py_ssize_t) size[1]);
				for(j=0; j<size[1]-1; j++)
					PyList_SET_ITEM(row,j,Py_BuildValue("d",data[i][j]));
				PyList_SET_ITEM(row,size[1]-1,Py_BuildValue("i",(int)data[i][size[1]-1]));
				PyList_SET_ITEM(result,i,row);
			}
			for(i=0; i<size[0]; i++)
				delete data[i];
			delete data;
		}
	}
	return result;
}

static PyMethodDef pyBalBifurcationDiagram_methods[] = {
	{"run", (PyCFunction) pyBalBifurcationDiagram_compute, METH_NOARGS, "Compute the bifurcation diagram"},
	{"summary", (PyCFunction) pyBalBifurcationDiagram_summary, 
		METH_VARARGS | METH_KEYWORDS, "Return or save to file the bifurcation diagram"},
	{NULL}	/* Sentinel */
};

static PyTypeObject pyBalBifurcationDiagramType = {
	PyObject_HEAD_INIT(NULL)
	0,									/* ob_size */
	"bal.BifurcationDiagram",			/* ob_name */
	sizeof(pyBalBifurcationDiagram),	/* tp_basicsize */
	0,									/* tp_itemsize */
	(destructor)pyBalBifurcationDiagram_dealloc,/* tp_dealloc */
	0,									/* tp_print */
	0,									/* tp_getattr */
	0,									/* tp_setattr */
	0,									/* tp_compare */
	0,									/* tp_repr */
	0,									/* tp_as_number */
	0,									/* tp_as_sequence */
	0,									/* tp_as_mapping */
	0,									/* tp_hash */
	0,									/* tp_call */
	0,									/* tp_str */
	(getattrofunc)pyBalBifurcationDiagram_getattro,/* tp_getattro  */
	(setattrofunc)pyBalBifurcationDiagram_setattro,/* tp_setattro  */
	0,									/* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,/* tp_flags */
	"Bifurcation Diagram object",		/* tp_doc */
	0,									/* tp_traverse */
    0,									/* tp_clear */
    0,									/* tp_richcompare */
    0,									/* tp_weaklistoffset */
    0,									/* tp_iter */
    0,									/* tp_iternext */
    pyBalBifurcationDiagram_methods,	/* tp_methods */
    0,									/* tp_members */
    0,									/* tp_getset */
    0,									/* tp_base */
    0,									/* tp_dict */
    0,									/* tp_descr_get */
    0,									/* tp_descr_set */
    0,									/* tp_dictoffset */
	(initproc)pyBalBifurcationDiagram_init,/* tp_init */
    0,									/* tp_alloc */
    pyBalBifurcationDiagram_new			/* tp_new */
};


/******************** Module stuff ********************/

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

PyMODINIT_FUNC initbal(void) {
	PyObject *m;

	ballib = dlopen(BALLIB, RTLD_LAZY);
	if(ballib == NULL)
		return;
	printf("\n\t BAL - Bifurcation Analysis Library\n");
	printf("\n\tSuccessfully loaded the library %s.\n", BALLIB);
	ballibext = dlopen(BALLIBEXT, RTLD_LAZY);
	if(ballibext != NULL)
		printf("\tSuccessfully loaded the library %s.\n", BALLIBEXT);
	printf("\n");
	
	if (PyType_Ready(&pyBalDynamicalSystemType) < 0)
		return;
	if (PyType_Ready(&pyBalParametersType) < 0)
		return;
	if (PyType_Ready(&pyBalODESolverType) < 0)
		return;
	if (PyType_Ready(&pyBalSolutionType) < 0)
		return;
	if (PyType_Ready(&pyBalBifurcationDiagramType) < 0)
		return;

	m = Py_InitModule3("pybal.bal", module_methods, "BAL module");
	if (m == NULL)
		return;
	
	Py_INCREF(&pyBalDynamicalSystemType);
	PyModule_AddObject(m, "DynamicalSystem", (PyObject *) &pyBalDynamicalSystemType);
	Py_INCREF(&pyBalParametersType);
	PyModule_AddObject(m, "Parameters", (PyObject *) &pyBalParametersType);
	Py_INCREF(&pyBalODESolverType);
	PyModule_AddObject(m, "ODESolver", (PyObject *) &pyBalODESolverType);
	Py_INCREF(&pyBalSolutionType);
	PyModule_AddObject(m, "Solution", (PyObject *) &pyBalSolutionType);
	Py_INCREF(&pyBalBifurcationDiagramType);
	PyModule_AddObject(m, "BifurcationDiagram", (PyObject *) &pyBalBifurcationDiagramType);
}

