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

#include "bfsiter.h"
#include "common.h"
#include "error.h"
#include "py2compat.h"
#include "vertexobject.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_bfsiter BFS iterator object
 */

PyTypeObject igraphmodule_BFSIterType;

/**
 * \ingroup python_interface_bfsiter
 * \brief Allocate a new BFS iterator object for a given graph and a given root
 * \param g the graph object being referenced
 * \param vid the root vertex index
 * \param advanced whether the iterator should be advanced (returning distance and parent as well)
 * \return the allocated PyObject
 */
PyObject* igraphmodule_BFSIter_new(igraphmodule_GraphObject *g, PyObject *root, igraph_neimode_t mode, igraph_bool_t advanced) {
  igraphmodule_BFSIterObject* o;
  long int no_of_nodes, r;
  
  o=PyObject_GC_New(igraphmodule_BFSIterObject, &igraphmodule_BFSIterType);
  Py_INCREF(g);
  o->gref=g;
  o->graph=&g->g;
  
  if (!PyInt_Check(root) && !PyObject_IsInstance(root, (PyObject*)&igraphmodule_VertexType)) {
    PyErr_SetString(PyExc_TypeError, "root must be integer or igraph.Vertex");
    return NULL;
  }
  
  no_of_nodes=igraph_vcount(&g->g);
  o->visited=(char*)calloc(no_of_nodes, sizeof(char));
  if (o->visited == 0) {
    PyErr_SetString(PyExc_MemoryError, "out of memory");
    return NULL;
  }
  
  if (igraph_dqueue_init(&o->queue, 100)) {
    PyErr_SetString(PyExc_MemoryError, "out of memory");
    return NULL;
  }
  if (igraph_vector_init(&o->neis, 0)) {
    PyErr_SetString(PyExc_MemoryError, "out of memory");
    igraph_dqueue_destroy(&o->queue);
    return NULL;
  }
  
  if (PyInt_Check(root)) {
    r=PyInt_AsLong(root);
  } else {
    r=((igraphmodule_VertexObject*)root)->idx;
  }
  if (igraph_dqueue_push(&o->queue, r) ||
      igraph_dqueue_push(&o->queue, 0) ||
      igraph_dqueue_push(&o->queue, -1)) {
    igraph_dqueue_destroy(&o->queue);
    igraph_vector_destroy(&o->neis);
    PyErr_SetString(PyExc_MemoryError, "out of memory");
    return NULL;
  }
  o->visited[r]=1;
  
  if (!igraph_is_directed(&g->g)) mode=IGRAPH_ALL;
  o->mode=mode;
  o->advanced=advanced;
  
  PyObject_GC_Track(o);
  
  RC_ALLOC("BFSIter", o);
  
  return (PyObject*)o;
}

/**
 * \ingroup python_interface_bfsiter
 * \brief Support for cyclic garbage collection in Python
 * 
 * This is necessary because the \c igraph.BFSIter object contains several
 * other \c PyObject pointers and they might point back to itself.
 */
int igraphmodule_BFSIter_traverse(igraphmodule_BFSIterObject *self,
				  visitproc visit, void *arg) {
  int vret;

  RC_TRAVERSE("BFSIter", self);
  
  if (self->gref) {
    vret=visit((PyObject*)self->gref, arg);
    if (vret != 0) return vret;
  }
  
  return 0;
}

/**
 * \ingroup python_interface_bfsiter
 * \brief Clears the iterator's subobject (before deallocation)
 */
int igraphmodule_BFSIter_clear(igraphmodule_BFSIterObject *self) {
  PyObject *tmp;

  PyObject_GC_UnTrack(self);
  
  tmp=(PyObject*)self->gref;
  self->gref=NULL;
  Py_XDECREF(tmp);

  igraph_dqueue_destroy(&self->queue);
  igraph_vector_destroy(&self->neis);
  free(self->visited);
  self->visited=0;
  
  return 0;
}

/**
 * \ingroup python_interface_bfsiter
 * \brief Deallocates a Python representation of a given BFS iterator object
 */
void igraphmodule_BFSIter_dealloc(igraphmodule_BFSIterObject* self) {
  igraphmodule_BFSIter_clear(self);

  RC_DEALLOC("BFSIter", self);
  
  PyObject_GC_Del(self);
}

PyObject* igraphmodule_BFSIter_iter(igraphmodule_BFSIterObject* self) {
  Py_INCREF(self);
  return (PyObject*)self;
}

PyObject* igraphmodule_BFSIter_iternext(igraphmodule_BFSIterObject* self) {
  if (!igraph_dqueue_empty(&self->queue)) {
    igraph_integer_t vid = (igraph_integer_t)igraph_dqueue_pop(&self->queue);
    igraph_integer_t dist = (igraph_integer_t)igraph_dqueue_pop(&self->queue);
    igraph_integer_t parent = (igraph_integer_t)igraph_dqueue_pop(&self->queue);
    long int i;
    
    if (igraph_neighbors(self->graph, &self->neis, vid, self->mode)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
	
    for (i=0; i<igraph_vector_size(&self->neis); i++) {
      igraph_integer_t neighbor = (igraph_integer_t)VECTOR(self->neis)[i];
      if (self->visited[neighbor]==0) {
	self->visited[neighbor]=1;
	if (igraph_dqueue_push(&self->queue, neighbor) ||
	    igraph_dqueue_push(&self->queue, dist+1) ||
	    igraph_dqueue_push(&self->queue, vid)) {
	  igraphmodule_handle_igraph_error();
	  return NULL;
	}
      }
    }

    if (self->advanced) {
      PyObject *vertexobj, *parentobj;
      vertexobj = igraphmodule_Vertex_New(self->gref, vid);
      if (!vertexobj)
        return NULL;
      if (parent >= 0) {
        parentobj = igraphmodule_Vertex_New(self->gref, parent);
        if (!parentobj)
            return NULL;
      } else {
        Py_INCREF(Py_None);
        parentobj=Py_None;
      }
      return Py_BuildValue("NlN", vertexobj, (long int)dist, parentobj);
    } else {
      return igraphmodule_Vertex_New(self->gref, vid);
    }
  } else {
    return NULL;
  }
}

/**
 * \ingroup python_interface_bfsiter
 * Method table for the \c igraph.BFSIter object
 */
PyMethodDef igraphmodule_BFSIter_methods[] = {
  {NULL}
};

/** \ingroup python_interface_bfsiter
 * Python type object referencing the methods Python calls when it performs various operations on
 * a BFS iterator of a graph
 */
PyTypeObject igraphmodule_BFSIterType =
{
  PyVarObject_HEAD_INIT(0, 0)
  "igraph.BFSIter",                         // tp_name
  sizeof(igraphmodule_BFSIterObject),       // tp_basicsize
  0,                                        // tp_itemsize
  (destructor)igraphmodule_BFSIter_dealloc, // tp_dealloc
  0,                                        // tp_print
  0,                                        // tp_getattr
  0,                                        // tp_setattr
  0,                                        /* tp_compare (2.x) / tp_reserved (3.x) */
  0,                                        // tp_repr
  0,                                        // tp_as_number
  0,                                        // tp_as_sequence
  0,                                        // tp_as_mapping
  0,                                        // tp_hash
  0,                                        // tp_call
  0,                                        // tp_str
  0,                                        // tp_getattro
  0,                                        // tp_setattro
  0,                                        // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, // tp_flags
  "igraph BFS iterator object",             // tp_doc
  (traverseproc) igraphmodule_BFSIter_traverse, /* tp_traverse */
  (inquiry) igraphmodule_BFSIter_clear,     /* tp_clear */
  0,                                        // tp_richcompare
  0,                                        // tp_weaklistoffset
  (getiterfunc)igraphmodule_BFSIter_iter,   /* tp_iter */
  (iternextfunc)igraphmodule_BFSIter_iternext, /* tp_iternext */
  0,                                        /* tp_methods */
  0,                                        /* tp_members */
  0,                                        /* tp_getset */
  0,                                        /* tp_base */
  0,                                        /* tp_dict */
  0,                                        /* tp_descr_get */
  0,                                        /* tp_descr_set */
  0,                                        /* tp_dictoffset */
  0,                                        /* tp_init */
  0,                                        /* tp_alloc */
  0,                                        /* tp_new */
  0,                                        /* tp_free */
};

