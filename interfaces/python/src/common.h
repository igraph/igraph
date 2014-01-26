/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Tamas Nepusz <ntamas@gmail.com>
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef PYTHON_COMMON_H
#define PYTHON_COMMON_H

#include <Python.h>

#ifdef RC_DEBUG
#  define RC_ALLOC(T, P) fprintf(stderr, "[ alloc ] " T " @ %p\n", P)
#  define RC_DECREF(T, P) fprintf(stderr, "[ ref - ] " T " @ %p (was: %d)\n", P, (int)P->ob_refcnt);
#  define RC_INCREF(T, P) fprintf(stderr, "[ ref + ] " T " @ %p (was: %d)\n", P, (int)P->ob_refcnt);
#  define RC_PRINT(P) fprintf(stderr, "[refcntr] %s @ %p = %d\n", ((PyTypeObject*)P->ob_type)->tp_name, P, (int)P->ob_refcnt);
#  define RC_DEALLOC(T, P) fprintf(stderr, "[dealloc] " T " @ %p\n", P);
#  define RC_TRAVERSE(T, P)
#else
#  define RC_ALLOC(T, P)
#  define RC_DECREF(T, P)
#  define RC_INCREF(T, P)
#  define RC_PRINT(P)
#  define RC_DEALLOC(T, P)
#  define RC_TRAVERSE(T, P)
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

#ifndef Py_RETURN_NOTIMPLEMENTED
#define Py_RETURN_NOTIMPLEMENTED { Py_INCREF(Py_NotImplemented); return Py_NotImplemented; }
#endif

#ifndef Py_RETURN
#define Py_RETURN(x) { if (x) { Py_RETURN_TRUE; } else { Py_RETURN_FALSE; } }
#endif

/* Compatibility stuff for Python 2.4 */
#if (PY_MAJOR_VERSION <= 2) & (PY_MINOR_VERSION <= 4)
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizessizeargfunc intintargfunc
#define Py_ssize_t int
#endif

#define ATTRIBUTE_TYPE_VERTEX 1
#define ATTRIBUTE_TYPE_EDGE 2

PyObject* igraphmodule_unimplemented(PyObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_resolve_graph_weakref(PyObject* ref);
#endif
