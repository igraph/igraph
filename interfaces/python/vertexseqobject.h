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

#ifndef PYTHON_VERTEXSEQOBJECT_H
#define PYTHON_VERTEXSEQOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_vertexseq
 * \brief A structure representing the vertex sequence of a graph
 */
typedef struct
{
  PyObject_HEAD
  PyObject* gref;
} igraphmodule_VertexSeqObject;

PyObject* igraphmodule_VertexSeq_New(igraphmodule_GraphObject *g);
int igraphmodule_VertexSeq_traverse(igraphmodule_VertexSeqObject *self,
					   visitproc visit, void *arg);
int igraphmodule_VertexSeq_clear(igraphmodule_VertexSeqObject *self);
void igraphmodule_VertexSeq_dealloc(igraphmodule_VertexSeqObject* self);

int igraphmodule_VertexSeq_sq_length(igraphmodule_VertexSeqObject *self);

static PyTypeObject igraphmodule_VertexSeqType;

#endif
