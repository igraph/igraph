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

#ifndef PYTHON_VERTEXSEQOBJECT_H
#define PYTHON_VERTEXSEQOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_vertexseq
 * \brief A structure representing the vertex sequence of a graph
 */
typedef struct {
  PyObject_HEAD
  igraphmodule_GraphObject* gref;
  igraph_vs_t vs;
  PyObject* weakreflist;
} igraphmodule_VertexSeqObject;

PyObject* igraphmodule_VertexSeq_new(PyTypeObject *subtype,
  PyObject* args, PyObject* kwds);
int igraphmodule_VertexSeq_init(igraphmodule_VertexSeqObject* self,
  PyObject* args, PyObject* kwds);
void igraphmodule_VertexSeq_dealloc(igraphmodule_VertexSeqObject* self);

int igraphmodule_VertexSeq_sq_length(igraphmodule_VertexSeqObject *self);

PyObject* igraphmodule_VertexSeq_find(igraphmodule_VertexSeqObject *self,
  PyObject *args);
PyObject* igraphmodule_VertexSeq_select(igraphmodule_VertexSeqObject *self,
  PyObject *args);

int igraphmodule_VertexSeq_to_vector_t(igraphmodule_VertexSeqObject *self,
  igraph_vector_t *v);
PyObject* igraphmodule_VertexSeq_get_graph(igraphmodule_VertexSeqObject *self,
  void* closure);

extern PyTypeObject igraphmodule_VertexSeqType;

#endif
