/* -*- mode: C -*-  */
/* vim:ts=2 sw=2 sts=2 et: */

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

PyTypeObject igraphmodule_EdgeType;

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
    /* RC_ALLOC("Edge", self); */
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

  /* RC_DEALLOC("Edge", self); */

  PyObject_GC_Del((PyObject*)self);
}

/** \ingroup python_interface_edge
 * \brief Formats an \c igraph.Edge object as a string
 * 
 * \return the formatted textual representation as a \c PyObject
 */
PyObject* igraphmodule_Edge_repr(igraphmodule_EdgeObject *self) {
  PyObject *o, *s, *grepr, *drepr;

  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  grepr=PyObject_Repr(o);
  if (!grepr) return NULL;
  drepr=PyObject_Repr(igraphmodule_Edge_attributes(self));
  if (!drepr) { Py_DECREF(grepr); return NULL; }
  s=PyString_FromFormat("igraph.Edge(%s,%ld,%s)", PyString_AsString(grepr),
    self->idx, PyString_AsString(drepr));
  Py_DECREF(grepr);
  Py_DECREF(drepr);
  return s;
}

/** \ingroup python_interface_edge
 * \brief Returns the number of edge attributes
 */
int igraphmodule_Edge_attribute_count(igraphmodule_EdgeObject* self) {
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return 0;
  if (!((PyObject**)o->g.attr)[1]) return 0;
  return PyDict_Size(((PyObject**)o->g.attr)[1]);
}

/** \ingroup python_interface_edge
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_Edge_attribute_names(igraphmodule_EdgeObject* self) {
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  return igraphmodule_Graph_edge_attributes(o);
}

/** \ingroup python_interface_edge
 * \brief Returns a dict with attribute names and values
 */
PyObject* igraphmodule_Edge_attributes(igraphmodule_EdgeObject* self) {
  igraphmodule_GraphObject *o;
  PyObject *names, *dict;
  long i, n;

  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  dict=PyDict_New();
  if (!dict) return NULL;

  names=igraphmodule_Graph_edge_attributes(o);
  if (!names) { Py_DECREF(dict); return NULL; }

  n=PyList_Size(names);
  for (i=0; i<n; i++) {
    PyObject *name = PyList_GetItem(names, i);
    if (name) {
      PyObject *dictit;
      dictit = PyDict_GetItem(((PyObject**)o->g.attr)[ATTRHASH_IDX_EDGE], name);
      if (dictit) {
        PyObject *value = PyList_GetItem(dictit, self->idx);
        if (value) {
          Py_INCREF(value);
          PyDict_SetItem(dict, name, value);
        }
      }
    }
  }

  return dict;
}

/** \ingroup python_interface_edge
 * \brief Returns the corresponding value to a given attribute of the edge
 * \param self the edge object
 * \param s the attribute name to be queried
 */
PyObject* igraphmodule_Edge_get_attribute(igraphmodule_EdgeObject* self,
                      PyObject* s) {
  igraphmodule_GraphObject *o;
  PyObject* result;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  result=PyDict_GetItem(((PyObject**)o->g.attr)[2], s);
  if (result) {
    /* result is a list, so get the element with index self->idx */
    if (!PyList_Check(result)) {
      PyErr_SetString(igraphmodule_InternalError, "Edge attribute dict member is not a list");
      return NULL;
    }
    result=PyList_GetItem(result, self->idx);
    Py_INCREF(result);
    return result;
  }
  
  /* result is NULL, check whether there was an error */
  if (!PyErr_Occurred())
    PyErr_SetString(PyExc_KeyError, "Attribute does not exist");
  return NULL;
}

/** \ingroup python_interface_edge
 * \brief Sets the corresponding value of a given attribute of the edge
 * \param self the edge object
 * \param k the attribute name to be set
 * \param v the value to be set
 * \return 0 if everything's ok, -1 in case of error
 */
int igraphmodule_Edge_set_attribute(igraphmodule_EdgeObject* self, PyObject* k, PyObject* v) {
  igraphmodule_GraphObject *o;
  PyObject* result;
  int r;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return -1;

  if (v==NULL)
    // we are deleting attribute
    return PyDict_DelItem(((PyObject**)o->g.attr)[2], k);
  
  result=PyDict_GetItem(((PyObject**)o->g.attr)[2], k);
  if (result) {
    /* result is a list, so set the element with index self->idx */
    if (!PyList_Check(result)) {
      PyErr_SetString(igraphmodule_InternalError, "Vertex attribute dict member is not a list");
      return -1;
    }
    /* we actually don't own a reference here to v, so we must increase
     * its reference count, because PyList_SetItem will "steal" a reference!
     * It took me 1.5 hours between London and Manchester to figure it out */
    Py_INCREF(v);
    r=PyList_SetItem(result, self->idx, v);
    if (r == -1) { Py_DECREF(v); }
    return r;
  }
  
  /* result is NULL, check whether there was an error */
  if (!PyErr_Occurred()) {
    /* no, there wasn't, so we must simply add the attribute */
    int n=(int)igraph_ecount(&o->g), i;
    result=PyList_New(n);
    for (i=0; i<n; i++) {
      if (i != self->idx) {
    Py_INCREF(Py_None);
    if (PyList_SetItem(result, i, Py_None) == -1) {
      Py_DECREF(Py_None);
      Py_DECREF(result);
      return -1;
    }
      } else {
    /* Same game with the reference count here */
    Py_INCREF(v);
    if (PyList_SetItem(result, i, v) == -1) {
      Py_DECREF(v);
      Py_DECREF(result);
      return -1;
    }
      }
    }
    if (PyDict_SetItem(((PyObject**)o->g.attr)[2], k, result) == -1) {
      Py_DECREF(result); /* TODO: is it needed here? maybe not! */
      return -1;
    }
    return 0;
  }
  
  return -1;
}

/**
 * \ingroup python_interface_edge
 * Returns the source node index of an edge
 */
PyObject* igraphmodule_Edge_get_from(igraphmodule_EdgeObject* self, void* closure) {
  PyObject *o;
  igraph_integer_t from, to;
  
  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  if (igraph_edge(&((igraphmodule_GraphObject*)o)->g, self->idx, &from, &to)) {
    igraphmodule_handle_igraph_error(); return NULL;
  }
  return PyInt_FromLong((long)from);
}

/**
 * \ingroup python_interface_edge
 * Returns the target node index of an edge
 */
PyObject* igraphmodule_Edge_get_to(igraphmodule_EdgeObject* self, void* closure) {
  PyObject *o;
  igraph_integer_t from, to;
  
  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  if (igraph_edge(&((igraphmodule_GraphObject*)o)->g, self->idx, &from, &to)) {
    igraphmodule_handle_igraph_error(); return NULL;
  }
  return PyInt_FromLong((long)to);
}

/**
 * \ingroup python_interface_edge
 * Returns the target node index of an edge
 */
PyObject* igraphmodule_Edge_get_tuple(igraphmodule_EdgeObject* self, void* closure) {
  PyObject *o;
  igraph_integer_t from, to;
  
  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  if (igraph_edge(&((igraphmodule_GraphObject*)o)->g, self->idx, &from, &to)) {
    igraphmodule_handle_igraph_error(); return NULL;
  }
  return Py_BuildValue("(ii)", (long)from, (long)to);
}

/**
 * \ingroup python_interface_edge
 * Method table for the \c igraph.Edge object
 */
PyMethodDef igraphmodule_Edge_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_Edge_attributes,
    METH_NOARGS,
    "attributes() -> list\n\n"
    "Returns a dict of attribute names and values for the edge\n"
  },
  {"attribute_names", (PyCFunction)igraphmodule_Edge_attribute_names,
      METH_NOARGS,
      "attribute_names() -> list\n\n"
      "Returns the list of edge attribute names\n"
  },
  {NULL}
};

/**
 * \ingroup python_interface_edge
 * Getter/setter table for the \c igraph.Edge object
 */
PyGetSetDef igraphmodule_Edge_getseters[] = {
  {"source", (getter)igraphmodule_Edge_get_from, NULL,
      "Source node index of this edge", NULL
  },
  {"target", (getter)igraphmodule_Edge_get_to, NULL,
      "Target node index of this edge", NULL
  },
  {"tuple", (getter)igraphmodule_Edge_get_tuple, NULL,
      "Source and target node index of this edge as a tuple", NULL
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
  (reprfunc)igraphmodule_Edge_repr,           // tp_repr
  0,                                          // tp_as_number
  0,                                          // tp_as_sequence
  &igraphmodule_Edge_as_mapping,              // tp_as_mapping
  0,                                          // tp_hash
  0,                                          // tp_call
  0,                                          // tp_str
  0,                                          // tp_getattro
  0,                                          // tp_setattro
  0,                                          // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, // tp_flags
  "Class representing a single edge in a graph.\n\n"
  "The edge is referenced by its index, so if the underlying graph\n"
  "changes, the semantics of the edge object might change as well\n"
  "(if the edge indices are altered in the original graph).\n\n"
  "The attributes of the edge can be accessed by using the edge\n"
  "as a hash:\n\n"
  "  >>> e[\"weight\"] = 2\n"
  "  >>> print e[\"weight\"]\n"
  "  2\n", // tp_doc
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

