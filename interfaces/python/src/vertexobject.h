/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

#ifndef PYTHON_VERTEXOBJECT_H
#define PYTHON_VERTEXOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_vertex
 * \brief A structure representing a vertex of a graph
 */
typedef struct
{
  PyObject_HEAD
  igraphmodule_GraphObject* gref;
  long idx;
} igraphmodule_VertexObject;

int igraphmodule_Vertex_clear(igraphmodule_VertexObject *self);
void igraphmodule_Vertex_dealloc(igraphmodule_VertexObject* self);

PyObject* igraphmodule_Vertex_New(igraphmodule_GraphObject *gref, long idx);
PyObject* igraphmodule_Vertex_repr(igraphmodule_VertexObject *self);
PyObject* igraphmodule_Vertex_attributes(igraphmodule_VertexObject* self);
PyObject* igraphmodule_Vertex_attribute_names(igraphmodule_VertexObject* self);
long igraphmodule_Vertex_get_index_long(igraphmodule_VertexObject* self);

extern PyTypeObject igraphmodule_VertexType;

#endif
