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

#include "edgeseqobject.h"
#include "edgeobject.h"
#include "common.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_edgeseq Edge sequence object
 */

PyTypeObject igraphmodule_EdgeSeqType;

/**
 * \ingroup python_interface_edgeseq
 * \brief Allocate a new edge sequence object for a given graph
 * \param g the graph object being referenced
 * \return the allocated PyObject
 */
PyObject* igraphmodule_EdgeSeq_New(igraphmodule_GraphObject *g) {
  igraphmodule_EdgeSeqObject* o;
  
  o=PyObject_GC_New(igraphmodule_EdgeSeqObject, &igraphmodule_EdgeSeqType);
  o->gref=PyWeakref_NewRef((PyObject*)g, NULL);
  //Py_INCREF(g);
  PyObject_GC_Track(o);
  
  RC_ALLOC("EdgeSeq", o);
  
  return (PyObject*)o;
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Support for cyclic garbage collection in Python
 * 
 * This is necessary because the \c igraph.EdgeSeq object contains several
 * other \c PyObject pointers and they might point back to itself.
 */
int igraphmodule_EdgeSeq_traverse(igraphmodule_EdgeSeqObject *self,
				  visitproc visit, void *arg) {
  int vret;

  RC_TRAVERSE("EdgeSeq", self);
  
  if (self->gref) {
    vret=visit(self->gref, arg);
    if (vret != 0) return vret;
  }
  
  return 0;
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Clears the graph object's subobject (before deallocation)
 */
int igraphmodule_EdgeSeq_clear(igraphmodule_EdgeSeqObject *self) {
  PyObject *tmp;

  PyObject_GC_UnTrack(self);
  
  tmp=self->gref;
  self->gref=NULL;
  Py_XDECREF(tmp);

  return 0;
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Deallocates a Python representation of a given edge sequence object
 */
void igraphmodule_EdgeSeq_dealloc(igraphmodule_EdgeSeqObject* self) {
  igraphmodule_EdgeSeq_clear(self);

  RC_DEALLOC("EdgeSeq", self);
  
  PyObject_GC_Del(self);
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Returns the length of the sequence (i.e. the number of edges in the graph)
 */
int igraphmodule_EdgeSeq_sq_length(igraphmodule_EdgeSeqObject* self) {
  igraph_t *g;
  
  g=&((igraphmodule_GraphObject*)PyWeakref_GetObject(self->gref))->g;
  if ((PyObject*)g == Py_None) return 0;
  
  return (int)igraph_ecount(g);
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Returns the item at the given index in the sequence
 */
PyObject* igraphmodule_EdgeSeq_sq_item(igraphmodule_EdgeSeqObject* self,
				       int i) {
  igraphmodule_GraphObject *o;
  igraph_t *g;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  g=&o->g;
  if (i<0 || i>=(int)igraph_ecount(g)) {
    PyErr_SetString(PyExc_IndexError, "edge index out of range");
    return NULL;
  }
  /// @todo caching
  return igraphmodule_Edge_New(self->gref, i);
}

/** \ingroup python_interface_edgeseq
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_EdgeSeq_attributes(igraphmodule_EdgeSeqObject* self) {
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  return igraphmodule_Graph_edge_attributes(o);
}

/** \ingroup python_interface_edgeseq
 * \brief Returns the list of values for a given attribute
 */
PyObject* igraphmodule_EdgeSeq_get_attribute_values(igraphmodule_EdgeSeqObject* self, PyObject* o) {
  igraphmodule_GraphObject *gr;
  PyObject* result;
  
  gr=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!gr) return NULL;
  
  result=PyDict_GetItem(((PyObject**)gr->g.attr)[2], o);
  if (result) {
    Py_INCREF(result);
    return result;
  }
  
  if (!PyErr_Occurred())
    PyErr_SetString(PyExc_KeyError, "Attribute does not exist");
  return NULL;
}

/**
 * \ingroup python_interface_edgeseq
 * Method table for the \c igraph.EdgeSeq object
 */
PyMethodDef igraphmodule_EdgeSeq_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_EdgeSeq_attributes,
      METH_NOARGS,
      "attributes() -> list\n\n"
      "Returns the attribute name list of the graph's edges\n"
  },
  {"get_attribute_values", (PyCFunction)igraphmodule_EdgeSeq_get_attribute_values,
      METH_O,
      "get_attribute_values(attrname) -> list\n\n"
      "Returns the value of a given edge attribute for all edges.\n\n"
      "@param attrname: the name of the attribute\n"
  },
  {NULL}
};

/**
 * \ingroup python_interface_edgeseq
 * This is the collection of functions necessary to implement the
 * edge sequence as a real sequence (e.g. allowing to reference
 * edges by indices)
 */
static PySequenceMethods igraphmodule_EdgeSeq_as_sequence = {
  (inquiry)igraphmodule_EdgeSeq_sq_length,
  0,               /* sq_concat */
  0,               /* sq_repeat */
  (intargfunc)igraphmodule_EdgeSeq_sq_item, /* sq_item */
  0,                                          /* sq_slice */
  0,                                          /* sq_ass_item */
  0,                                          /* sq_ass_slice */
  0,                                          /* sq_contains */
  0,                                          /* sq_inplace_concat */
  0,                                          /* sq_inplace_repeat */
};

/** \ingroup python_interface_edgeseq
 * Python type object referencing the methods Python calls when it performs various operations on
 * an edge sequence of a graph
 */
PyTypeObject igraphmodule_EdgeSeqType =
{
  PyObject_HEAD_INIT(NULL)                  //
  0,                                        // ob_size
  "igraph.EdgeSeq",                         // tp_name
  sizeof(igraphmodule_EdgeSeqObject),       // tp_basicsize
  0,                                        // tp_itemsize
  (destructor)igraphmodule_EdgeSeq_dealloc, // tp_dealloc
  0,                                        // tp_print
  0,                                        // tp_getattr
  0,                                        // tp_setattr
  0,                                        // tp_compare
  0,                                        // tp_repr
  0,                                        // tp_as_number
  &igraphmodule_EdgeSeq_as_sequence,        // tp_as_sequence
  0,                                        // tp_as_mapping
  0,                                        // tp_hash
  0,                                        // tp_call
  0,                                        // tp_str
  0,                                        // tp_getattro
  0,                                        // tp_setattro
  0,                                        // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, // tp_flags
  "Class representing the edge set of a graph.\n\n"
  "The individual edges can be accessed by indexing the edge sequence\n"
  "object. It can be used as an iterable as well, or even in a list\n"
  "comprehension:\n\n"
  "  >>> g=igraph.Graph.Full(3)\n"
  "  >>> for e in g.es:\n"
  "  ...   print e.tuple\n"
  "  ...\n"
  "  (0, 1)\n"
  "  (0, 2)\n"
  "  (1, 2)\n"
  "  >>> [max(e.tuple) for e in g.es]\n"
  "  [1, 2, 2]\n", // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  0,                                          // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  igraphmodule_EdgeSeq_methods,             // tp_methods
};

