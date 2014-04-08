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

#ifndef PYTHON_EDGEOBJECT_H
#define PYTHON_EDGEOBJECT_H

#include <Python.h>
#include "graphobject.h"
#include "py2compat.h"

/**
 * \ingroup python_interface_edge
 * \brief A structure representing an edge of a graph
 */
typedef struct {
  PyObject_HEAD
  igraphmodule_GraphObject* gref;
  igraph_integer_t idx;
  Py_hash_t hash;
} igraphmodule_EdgeObject;

int igraphmodule_Edge_clear(igraphmodule_EdgeObject *self);
void igraphmodule_Edge_dealloc(igraphmodule_EdgeObject* self);

int igraphmodule_Edge_Check(PyObject *obj);
int igraphmodule_Edge_Validate(PyObject *obj);

PyObject* igraphmodule_Edge_New(igraphmodule_GraphObject *gref, igraph_integer_t idx);
PyObject* igraphmodule_Edge_repr(igraphmodule_EdgeObject *self);
PyObject* igraphmodule_Edge_attributes(igraphmodule_EdgeObject* self);
PyObject* igraphmodule_Edge_attribute_names(igraphmodule_EdgeObject* self);
igraph_integer_t igraphmodule_Edge_get_index_igraph_integer(igraphmodule_EdgeObject* self);
long igraphmodule_Edge_get_index_long(igraphmodule_EdgeObject* self);
PyObject* igraphmodule_Edge_update_attributes(PyObject* self, PyObject* args,
    PyObject* kwds);

extern PyTypeObject igraphmodule_EdgeType;

#endif
