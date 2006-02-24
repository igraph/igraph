#ifndef PYTHON_COMMON_H
#define PYTHON_COMMON_H

#include <Python.h>

#ifdef RC_DEBUG
#define RC_ALLOC(T, P) fprintf(stderr, "[ alloc ] " T " @ %p\n", P)
#define RC_DEALLOC(T, P) fprintf(stderr, "[dealloc] " T " @ %p\n", P);
#define RC_TRAVERSE(T, P)
//#define RC_TRAVERSE(T, P) fprintf(stderr, "[ travr ] " T " @ %p\n", P);
#else
#define RC_ALLOC(T, P)
#define RC_DEALLOC(T, P)
#define RC_TRAVERSE(T, P)
#endif

/* Compatibility stuff for Python 2.3 */
#ifndef Py_RETURN_TRUE
#define Py_RETURN_TRUE { Py_INCREF(Py_True); return Py_True; }
#endif

#ifndef Py_RETURN_FALSE
#define Py_RETURN_FALSE { Py_INCREF(Py_False); return Py_False; }
#endif

#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE { Py_INCREF(Py_None); return Py_None; }
#endif


PyObject* igraphmodule_unimplemented(PyObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_resolve_graph_weakref(PyObject* ref);
#endif
