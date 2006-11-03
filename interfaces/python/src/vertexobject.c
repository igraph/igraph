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

#include "vertexobject.h"
#include "graphobject.h"
#include "error.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_vertex Vertex object
 */

PyTypeObject igraphmodule_VertexType;

/**
 * \ingroup python_interface_vertex
 * \brief Allocates a new Python vertex object
 * \param gref weak reference of the \c igraph.Graph being referenced by the vertex
 * \param idx the index of the vertex
 * 
 * \warning \c igraph references its vertices by indices, so if
 * you delete some vertices from the graph, the vertex indices will
 * change. Since the \c igraph.Vertex objects do not follow these
 * changes, your existing vertex objects will point to elsewhere
 * (or they might even get invalidated).
 */
PyObject* igraphmodule_Vertex_New(PyObject *gref, long idx) {
  igraphmodule_VertexObject* self;
  self=PyObject_GC_New(igraphmodule_VertexObject, &igraphmodule_VertexType);
  if (self) {
    RC_ALLOC("Vertex", self);
    Py_INCREF(gref);
    self->gref=gref;
    self->idx=idx;
  }
  return (PyObject*)self;
}

/**
 * \ingroup python_interface_vertex
 * \brief Support for cyclic garbage collection in Python
 */
int igraphmodule_Vertex_traverse(igraphmodule_VertexObject *self,
				    visitproc visit, void *arg) {
  int vret;
  
  if (self->gref) {
    vret=visit((PyObject*)self->gref, arg);
    if (vret != 0) return vret;
  }
  
  return 0;
}

/**
 * \ingroup python_interface_vertex
 * \brief Clears the graph object's subobject (before deallocation)
 */
int igraphmodule_Vertex_clear(igraphmodule_VertexObject *self) {
  PyObject *tmp;

  tmp=(PyObject*)self->gref;
  self->gref=NULL;
  Py_XDECREF(tmp);

  return 0;
}

/**
 * \ingroup python_interface_vertex
 * \brief Deallocates a Python representation of a given vertex object
 */
void igraphmodule_Vertex_dealloc(igraphmodule_VertexObject* self) {
  igraphmodule_Vertex_clear(self);

  RC_DEALLOC("Vertex", self);

  PyObject_GC_Del((PyObject*)self);
}

/** \ingroup python_interface_vertex
 * \brief Formats an \c igraph.Vertex object in a human-consumable format.
 * 
 * \return the formatted textual representation as a \c PyObject
 */
PyObject* igraphmodule_Vertex_str(igraphmodule_VertexObject *self)
{
  PyObject *s, *o;
  
  o=igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  s=PyString_FromFormat("Vertex #%ld of: ", self->idx);
  PyString_Concat(&s, igraphmodule_Graph_str((igraphmodule_GraphObject*)o));
  return s;
}

/** \ingroup python_interface_vertex
 * \brief Returns the number of vertex attributes
 */
int igraphmodule_Vertex_attribute_count(igraphmodule_VertexObject* self) {
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return 0;
  if (!((PyObject**)o->g.attr)[1]) return 0;
  return PyDict_Size(((PyObject**)o->g.attr)[1]);
}

/** \ingroup python_interface_vertex
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_Vertex_attributes(igraphmodule_VertexObject* self) {
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  return igraphmodule_Graph_vertex_attributes(o);
}

/** \ingroup python_interface_vertex
 * \brief Returns the corresponding value to a given attribute of the vertex
 * \param self the vertex object
 * \param s the attribute name to be queried
 */
PyObject* igraphmodule_Vertex_get_attribute(igraphmodule_VertexObject* self,
					   PyObject* s) {
  igraphmodule_GraphObject *o;
  PyObject* result;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  result=PyDict_GetItem(((PyObject**)o->g.attr)[1], s);
  if (result) {
    /* result is a list, so get the element with index self->idx */
    if (!PyList_Check(result)) {
      PyErr_SetString(igraphmodule_InternalError, "Vertex attribute dict member is not a list");
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

/** \ingroup python_interface_vertex
 * \brief Sets the corresponding value of a given attribute of the vertex
 * \param self the vertex object
 * \param k the attribute name to be set
 * \param v the value to be set
 * \return 0 if everything's ok, -1 in case of error
 */
int igraphmodule_Vertex_set_attribute(igraphmodule_VertexObject* self, PyObject* k, PyObject* v) {
  igraphmodule_GraphObject *o;
  PyObject* result;
  int r;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return -1;

  if (v==NULL)
    // we are deleting attribute
    return PyDict_DelItem(((PyObject**)o->g.attr)[1], k);
  
  result=PyDict_GetItem(((PyObject**)o->g.attr)[1], k);
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
    if (r == -1) Py_DECREF(v);
    return r;
  }
  
  /* result is NULL, check whether there was an error */
  if (!PyErr_Occurred()) {
    /* no, there wasn't, so we must simply add the attribute */
    int n=(int)igraph_vcount(&o->g), i;
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
    if (PyDict_SetItem(((PyObject**)o->g.attr)[1], k, result) == -1) {
      Py_DECREF(result);
      return -1;
    }
    return 0;
  }
  
  return -1;
}

/**
 * \ingroup python_interface_vertex
 * Returns the vertex index
 */
PyObject* igraphmodule_Vertex_get_index(igraphmodule_VertexObject* self, void* closure) {
  return PyInt_FromLong((long)self->idx);
}

/**
 * \ingroup python_interface_vertex
 * Method table for the \c igraph.Vertex object
 */
PyMethodDef igraphmodule_Vertex_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_Vertex_attributes,
      METH_NOARGS,
      "attributes() -> list\n\n"
      "Returns the attribute list of the graph's vertices\n"
  },
  {NULL}
};

/** \ingroup python_interface_vertex
 * This structure is the collection of functions necessary to implement
 * the vertex as a mapping (i.e. to allow the retrieval and setting of
 * igraph attributes in Python as if it were of a Python mapping type)
 */
PyMappingMethods igraphmodule_Vertex_as_mapping = {
  // returns the number of vertex attributes
  (inquiry)igraphmodule_Vertex_attribute_count,
  // returns an attribute by name
  (binaryfunc)igraphmodule_Vertex_get_attribute,
  // sets an attribute by name
  (objobjargproc)igraphmodule_Vertex_set_attribute
};

/**
 * \ingroup python_interface_vertex
 * Getter/setter table for the \c igraph.Vertex object
 */
PyGetSetDef igraphmodule_Vertex_getseters[] = {
  {"index", (getter)igraphmodule_Vertex_get_index, NULL,
      "Index of the vertex", NULL
  },
  {NULL}
};

/** \ingroup python_interface_vertex
 * Python type object referencing the methods Python calls when it performs various operations on
 * a vertex of a graph
 */
PyTypeObject igraphmodule_VertexType =
{
  PyObject_HEAD_INIT(NULL)                    //
  0,                                          // ob_size
  "igraph.Vertex",                            // tp_name
  sizeof(igraphmodule_VertexObject),          // tp_basicsize
  0,                                          // tp_itemsize
  (destructor)igraphmodule_Vertex_dealloc,    // tp_dealloc
  0,                                          // tp_print
  0,                                          // tp_getattr
  0,                                          // tp_setattr
  0,                                          // tp_compare
  0,                                          // tp_repr
  0,                                          // tp_as_number
  0,                                          // tp_as_sequence
  &igraphmodule_Vertex_as_mapping,            // tp_as_mapping
  0,                                          // tp_hash
  0,                                          // tp_call
  (reprfunc)igraphmodule_Vertex_str,          // tp_str
  0,                                          // tp_getattro
  0,                                          // tp_setattro
  0,                                          // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, // tp_flags
  "Class representing a single vertex in a graph.\n\n"
  "The vertex is referenced by its index, so if the underlying graph\n"
  "changes, the semantics of the vertex object might change as well\n"
  "(if the vertex indices are altered in the original graph).\n\n"
  "The attributes of the vertex can be accessed by using the vertex\n"
  "as a hash:\n\n"
  "  >>> v[\"color\"] = \"red\"\n"
  "  >>> print v[\"color\"]\n"
  "  red\n", // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  0,                                          // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  igraphmodule_Vertex_methods,                // tp_methods
  0,                                          // tp_members
  igraphmodule_Vertex_getseters,              // tp_getset
};

