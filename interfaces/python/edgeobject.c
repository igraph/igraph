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

#include "edgeobject.h"
#include "graphobject.h"
#include "error.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_edge Edge object
 */

/**
 * \ingroup python_interface_edge
 * \brief Allocates a new Python edge object
 * \param gref weak reference of the \c igraph.Graph being referenced by the edge
 * \param idx the index of the edge
 * 
 * \warning \c igraph references its edges by indices, so if
 * you delete some edges from the graph, the edge indices will
 * change. Since the \c igraph.Edge objects do not follow these
 * changes, your existing edge objects will point to elsewhere
 * (or they might even get invalidated).
 */
PyObject* igraphmodule_Edge_New(PyObject *gref, long idx) {
  igraphmodule_EdgeObject* self;
  self=PyObject_GC_New(igraphmodule_EdgeObject, &igraphmodule_EdgeType);
  if (self) {
    RC_ALLOC("Edge", self);
    Py_INCREF(gref);
    self->gref=gref;
    self->idx=idx;
  }
  return (PyObject*)self;
}

/**
 * \ingroup python_interface_edge
 * \brief Support for cyclic garbage collection in Python
 */
int igraphmodule_Edge_traverse(igraphmodule_EdgeObject *self,
			       visitproc visit, void *arg) {
  int vret;
  
  if (self->gref) {
    vret=visit((PyObject*)self->gref, arg);
    if (vret != 0) return vret;
  }
  
  return 0;
}

/**
 * \ingroup python_interface_edge
 * \brief Clears the graph object's subobject (before deallocation)
 */
int igraphmodule_Edge_clear(igraphmodule_EdgeObject *self) {
  PyObject *tmp;

  tmp=(PyObject*)self->gref;
  self->gref=NULL;
  Py_XDECREF(tmp);

  return 0;
}

/**
 * \ingroup python_interface_edge
 * \brief Deallocates a Python representation of a given edge object
 */
void igraphmodule_Edge_dealloc(igraphmodule_EdgeObject* self) {
  igraphmodule_Edge_clear(self);

  RC_DEALLOC("Edge", self);

  PyObject_GC_Del((PyObject*)self);
}

/** \ingroup python_interface_edge
 * \brief Formats an \c igraph.Edge object in a human-consumable format.
 * 
 * \return the formatted textual representation as a \c PyObject
 */
PyObject* igraphmodule_Edge_str(igraphmodule_EdgeObject *self)
{
  PyObject *s, *o;
  
  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  s=PyString_FromFormat("Edge #%ld of: ", self->idx);
  PyString_Concat(&s, igraphmodule_Graph_str((igraphmodule_GraphObject*)o));
  return s;
}

/** \ingroup python_interface_edge
 * \brief Returns the number of edge attributes
 */
int igraphmodule_Edge_attribute_count(igraphmodule_EdgeObject* self) {
  igraph_vector_t t;
  long result;
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return 0;
  
  if (igraph_vector_init(&t, 0)) return 0;
  if (igraph_list_edge_attributes(&o->g, NULL, &t)) {
    igraph_vector_destroy(&t);
    return 0;
  }
  
  result=igraph_vector_size(&t);
  igraph_vector_destroy(&t);
  return result;
}

/** \ingroup python_interface_edge
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_Edge_attributes(igraphmodule_EdgeObject* self) {
  igraph_vector_t t;
  igraph_vector_ptr_t ns;
  long result;
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  return igraphmodule_Graph_edge_attributes(o, NULL, NULL);
}

/** \ingroup python_interface_edge
 * \brief Returns the corresponding value to a given attribute of the edge
 * \param self the edge object
 * \param s the attribute name to be queried
 */
PyObject* igraphmodule_Edge_get_attribute(igraphmodule_EdgeObject* self,
					  PyObject* s) {
  igraph_attribute_type_t t=-1;
  void* value=NULL;
  igraphmodule_GraphObject *o;
  int result;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  if (!PyString_Check(s)) {
    PyErr_SetString(PyExc_TypeError, "Attribute name must be string");
    return NULL;
  }
  
  result=igraph_get_edge_attribute(&o->g, PyString_AsString(s),
				     self->idx, &value, &t);
  if (result == IGRAPH_EINVAL) {
    PyErr_SetString(PyExc_KeyError, "Attribute does not exist");
    return NULL;
  } else if (result) {
    return igraphmodule_handle_igraph_error();
  }
  
  if (t==IGRAPH_ATTRIBUTE_NUM)
    return PyFloat_FromDouble((double)(*(real_t*)value));
  if (t==IGRAPH_ATTRIBUTE_STR)
    return PyString_FromString((char*)value);
  return igraphmodule_handle_igraph_error();
}

/** \ingroup python_interface_edge
 * \brief Sets the corresponding value of a given attribute of the edge
 * \param self the edge object
 * \param k the attribute name to be set
 * \param v the value to be set
 * \return 0 if everything's ok, -1 in case of error
 */
int igraphmodule_Edge_set_attribute(igraphmodule_EdgeObject* self, PyObject* k, PyObject* v) {
  igraph_attribute_type_t t=-1, t2=-1;
  void* value=NULL;
  char* key;
  real_t value0;
  igraphmodule_GraphObject *o;
  int result;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return -1;
  
  if (!PyString_Check(k)) {
    PyErr_SetString(PyExc_TypeError, "Attribute name must be string");
    return -1;
  }
  if (v!=NULL && !PyString_Check(v) && !PyInt_Check(v) && !PyLong_Check(v) && !PyFloat_Check(v)) {
    PyErr_SetString(PyExc_TypeError, "Attribute value must be either numeric or string");
    return -1;
  }
  
  key=PyString_AsString(k);
  result=igraph_get_edge_attribute_type(&o->g, key, &t2);
  if (result == IGRAPH_EINVAL) {
    t2=-1;         // to indicate that there's no such attribute yet
    PyErr_Clear(); // to indicate that we handled the situation
  } else if (result) {
    igraphmodule_handle_igraph_error(); return -1;
  }
  
  if (v==NULL) {
    // we are deleting attribute
    if (igraph_remove_edge_attribute(&o->g, key)) {
      igraphmodule_handle_igraph_error(); return -1;
    }
    return 0;
  } else if (PyString_Check(v)) {
    t=IGRAPH_ATTRIBUTE_STR;
    value=(void*)PyString_AsString(v);
  } else {
    t=IGRAPH_ATTRIBUTE_NUM;
    if (PyInt_Check(v))
      value0=(real_t)(PyInt_AsLong(v));
    else if (PyLong_Check(v))
      value0=(real_t)(PyLong_AsLong(v));
    else if (PyFloat_Check(v))
      value0=(real_t)(PyFloat_AsDouble(v));
    if (PyErr_Occurred()) return -1;
    value=(void*)&value0;
  }

  if ((long)t2 == -1) {
    // Attribute does not exist yet, so we add it
    if (igraph_add_edge_attribute(&o->g, key, t)) {
      igraphmodule_handle_igraph_error(); return -1;
    }
    t2=t;
  } else {
    if (t2!=t) {
      /* The type of the existing attribute differs from the new one */
      PyErr_SetString(PyExc_TypeError, "Edge attribute type does not match the type of the given value");
      return -1;
    }
  }
  if (igraph_set_edge_attribute(&o->g, key, self->idx, value)) {
    igraphmodule_handle_igraph_error(); return -1;
  }
  return 0;
}

/**
 * \ingroup python_interface_edge
 * Returns the source node index of an edge
 */
PyObject* igraphmodule_Edge_get_from(igraphmodule_EdgeObject* self, void* closure) {
  PyObject *o;
  integer_t from, to;
  
  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  if (igraph_edge(&((igraphmodule_GraphObject*)o)->g, self->idx, &from, &to)) {
    igraphmodule_handle_igraph_error(); return NULL;
  }
  return PyInt_FromLong(from);
}

/**
 * \ingroup python_interface_edge
 * Returns the target node index of an edge
 */
PyObject* igraphmodule_Edge_get_to(igraphmodule_EdgeObject* self, void* closure) {
  PyObject *o;
  integer_t from, to;
  
  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  if (igraph_edge(&((igraphmodule_GraphObject*)o)->g, self->idx, &from, &to)) {
    igraphmodule_handle_igraph_error(); return NULL;
  }
  return PyInt_FromLong(to);
}

/**
 * \ingroup python_interface_edge
 * Method table for the \c igraph.Edge object
 */
PyMethodDef igraphmodule_Edge_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_Edge_attributes,
      METH_NOARGS,
      "Returns the attribute list of the graph's edges\n"
  },
  {NULL}
};

/**
 * \ingroup python_interface_edge
 * Getter/setter table for the \c igraph.Edge object
 */
PyGetSetDef igraphmodule_Edge_getseters[] = {
  {"source", igraphmodule_Edge_get_from, NULL,
      "Source node index of this edge", NULL
  },
  {"target", igraphmodule_Edge_get_to, NULL,
      "Target node index of this edge", NULL
  },
  {NULL}
};

/** \ingroup python_interface_edge
 * This structure is the collection of functions necessary to implement
 * the edge as a mapping (i.e. to allow the retrieval and setting of
 * igraph attributes in Python as if it were of a Python mapping type)
 */
PyMappingMethods igraphmodule_Edge_as_mapping = {
  // returns the number of edge attributes
  (inquiry)igraphmodule_Edge_attribute_count,
  // returns an attribute by name
  (binaryfunc)igraphmodule_Edge_get_attribute,
  // sets an attribute by name
  (objobjargproc)igraphmodule_Edge_set_attribute
};

/** \ingroup python_interface_edge
 * Python type object referencing the methods Python calls when it performs various operations on
 * an edge of a graph
 */
PyTypeObject igraphmodule_EdgeType =
{
  PyObject_HEAD_INIT(NULL)                    //
  0,                                          // ob_size
  "igraph.Edge",                              // tp_name
  sizeof(igraphmodule_EdgeObject),            // tp_basicsize
  0,                                          // tp_itemsize
  (destructor)igraphmodule_Edge_dealloc,      // tp_dealloc
  0,                                          // tp_print
  0,                                          // tp_getattr
  0,                                          // tp_setattr
  0,                                          // tp_compare
  0,                                          // tp_repr
  0,                                          // tp_as_number
  0,                                          // tp_as_sequence
  &igraphmodule_Edge_as_mapping,              // tp_as_mapping
  0,                                          // tp_hash
  0,                                          // tp_call
  (reprfunc)igraphmodule_Edge_str,            // tp_str
  0,                                          // tp_getattro
  0,                                          // tp_setattro
  0,                                          // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, // tp_flags
  "igraph edge object",                       // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  0,                                          // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  igraphmodule_Edge_methods,                  // tp_methods
  0,                                          // tp_members
  igraphmodule_Edge_getseters,                // tp_getset
};

