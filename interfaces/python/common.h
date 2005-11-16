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

PyObject* igraphmodule_unimplemented(PyObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_resolve_graph_weakref(PyObject* ref);
#endif
