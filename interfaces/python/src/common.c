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

#include "common.h"
#include "structmember.h"

/**
 * \ingroup python_interface
 * \brief Handler function for all unimplemented \c igraph.Graph methods
 * 
 * This function is called whenever an unimplemented \c igraph.Graph method
 * is called ("unimplemented" meaning that there is a method name in the
 * method table of \c igraph.Graph , but there isn't any working implementation
 * either because the underlying \c igraph API might be subject to change
 * or because the calling format from Python is not decided yet (or maybe
 * because of laziness or lack of time ;))
 * 
 * All of the parameters are ignored, they are here just to make the
 * function satisfy the requirements of \c PyCFunction, thus allowing it
 * to be included in a method table.
 * 
 * \return NULL
 */
PyObject* igraphmodule_unimplemented(PyObject* self, PyObject* args, PyObject* kwds)
{
   PyErr_SetString(PyExc_NotImplementedError, "This method is unimplemented.");
   return NULL;
}

/**
 * \ingroup python_interface
 * \brief Resolves a weak reference to an \c igraph.Graph
 * \return the \c igraph.Graph object or NULL if the weak reference is dead.
 * Sets an exception in the latter case.
 */
PyObject* igraphmodule_resolve_graph_weakref(PyObject* ref) {
  PyObject *o;
  
  if (!PyWeakref_Check(ref)) {
    PyErr_SetString(PyExc_TypeError, "weak reference expected");
    return NULL;
  }
  o=PyWeakref_GetObject(ref);
  if (o == Py_None) {
    PyErr_SetString(PyExc_TypeError, "underlying graph has already been destroyed");
    return NULL;
  }
  return o;
}
