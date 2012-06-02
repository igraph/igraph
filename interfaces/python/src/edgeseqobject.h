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

#ifndef PYTHON_EDGESEQOBJECT_H
#define PYTHON_EDGESEQOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_edgeseq
 * \brief A structure representing the edge sequence of a graph
 */
typedef struct
{
  PyObject_HEAD
  igraphmodule_GraphObject* gref;
  igraph_es_t es;
  PyObject* weakreflist;
} igraphmodule_EdgeSeqObject;

PyObject* igraphmodule_EdgeSeq_new(PyTypeObject *subtype,
  PyObject *args, PyObject *kwds);
igraphmodule_EdgeSeqObject* igraphmodule_EdgeSeq_copy(
  igraphmodule_EdgeSeqObject *o);
int igraphmodule_EdgeSeq_init(igraphmodule_EdgeSeqObject *self,
  PyObject *args, PyObject *kwds);
void igraphmodule_EdgeSeq_dealloc(igraphmodule_EdgeSeqObject* self);

int igraphmodule_EdgeSeq_sq_length(igraphmodule_EdgeSeqObject *self);

PyObject* igraphmodule_EdgeSeq_find(igraphmodule_EdgeSeqObject *self,
  PyObject *args);
PyObject* igraphmodule_EdgeSeq_select(igraphmodule_EdgeSeqObject *self,
  PyObject *args);

PyObject* igraphmodule_EdgeSeq_get_graph(igraphmodule_EdgeSeqObject *self,
  void* closure);

extern PyTypeObject igraphmodule_EdgeSeqType;

#endif
