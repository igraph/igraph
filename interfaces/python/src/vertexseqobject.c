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

#include "vertexseqobject.h"
#include "vertexobject.h"
#include "common.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_vertexseq Vertex sequence object
 */

PyTypeObject igraphmodule_VertexSeqType;

/**
 * \ingroup python_interface_vertexseq
 * \brief Allocate a new vertex sequence object for a given graph
 * \param g the graph object being referenced
 * \return the allocated PyObject
 */
PyObject* igraphmodule_VertexSeq_New(igraphmodule_GraphObject *g) {
  igraphmodule_VertexSeqObject* o;
  
  o=PyObject_GC_New(igraphmodule_VertexSeqObject, &igraphmodule_VertexSeqType);
  o->gref=PyWeakref_NewRef((PyObject*)g, NULL);
  //Py_INCREF(g);
  PyObject_GC_Track(o);
  
  RC_ALLOC("VertexSeq", o);
  
  return (PyObject*)o;
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Support for cyclic garbage collection in Python
 * 
 * This is necessary because the \c igraph.VertexSeq object contains several
 * other \c PyObject pointers and they might point back to itself.
 */
int igraphmodule_VertexSeq_traverse(igraphmodule_VertexSeqObject *self,
                    visitproc visit, void *arg) {
  int vret;

  RC_TRAVERSE("VertexSeq", self);
  
  if (self->gref) {
    vret=visit(self->gref, arg);
    if (vret != 0) return vret;
  }
  
  return 0;
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Clears the graph object's subobject (before deallocation)
 */
int igraphmodule_VertexSeq_clear(igraphmodule_VertexSeqObject *self) {
  PyObject *tmp;

  PyObject_GC_UnTrack(self);
  
  tmp=self->gref;
  self->gref=NULL;
  Py_XDECREF(tmp);

  return 0;
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Deallocates a Python representation of a given vertex sequence object
 */
void igraphmodule_VertexSeq_dealloc(igraphmodule_VertexSeqObject* self) {
  igraphmodule_VertexSeq_clear(self);

  RC_DEALLOC("VertexSeq", self);
  
  PyObject_GC_Del(self);
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Returns the length of the sequence (i.e. the number of vertices in the graph)
 */
int igraphmodule_VertexSeq_sq_length(igraphmodule_VertexSeqObject* self) {
  igraph_t *g;
  
  g=&((igraphmodule_GraphObject*)PyWeakref_GetObject(self->gref))->g;
  if ((PyObject*)g == Py_None) return 0;
  
  return (int)igraph_vcount(g);
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Returns the item at the given index in the sequence
 */
PyObject* igraphmodule_VertexSeq_sq_item(igraphmodule_VertexSeqObject* self,
                     int i) {
  igraphmodule_GraphObject *o;
  igraph_t *g;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  g=&o->g;
  if (i<0 || i>=(int)igraph_vcount(g)) {
    PyErr_SetString(PyExc_IndexError, "vertex index out of range");
    return NULL;
  }
  /// @todo caching
  return igraphmodule_Vertex_New(self->gref, i);
}

/** \ingroup python_interface_vertexseq
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_VertexSeq_attributes(igraphmodule_VertexSeqObject* self) {
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  return igraphmodule_Graph_vertex_attributes(o);
}

/** \ingroup python_interface_vertexseq
 * \brief Returns the count of attribute names
 */
PyObject* igraphmodule_VertexSeq_attribute_count(igraphmodule_VertexSeqObject* self) {
  PyObject *list;
    long int size;
    igraphmodule_GraphObject *o;
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  list=igraphmodule_Graph_vertex_attributes(o);
  size=PySequence_Size(list);
  Py_DECREF(list);

  return Py_BuildValue("i", size);
}

/** \ingroup python_interface_vertexseq
 * \brief Returns the list of values for a given attribute
 */
PyObject* igraphmodule_VertexSeq_get_attribute_values(igraphmodule_VertexSeqObject* self, PyObject* o) {
  igraphmodule_GraphObject *gr;
  PyObject *result, *values, *item;
  
  gr=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!gr) return NULL;
  
  values=PyDict_GetItem(((PyObject**)gr->g.attr)[ATTRHASH_IDX_VERTEX], o);
  if (values) {
    long int i, n;
    n = PyList_Size(values);
    result = PyList_New(n);
    if (!result) return result;

    for (i=0; i<n; i++) {
        item = PyList_GET_ITEM(values, i);
        Py_INCREF(item);
        PyList_SET_ITEM(result, i, item);
    }
    return result;
  }
  
  if (!PyErr_Occurred())
    PyErr_SetString(PyExc_KeyError, "Attribute does not exist");
  return NULL;
}

PyObject* igraphmodule_VertexSeq_get_attribute_values_mapping(igraphmodule_VertexSeqObject *self, PyObject *o) {
  /* Handle integer indices according to the sequence protocol */
  if (PyInt_Check(o)) return igraphmodule_VertexSeq_sq_item(self, PyInt_AsLong(o));
  return igraphmodule_VertexSeq_get_attribute_values(self, o);
}

/** \ingroup python_interface_vertexseq
 * \brief Sets the list of values for a given attribute
 */
int igraphmodule_VertexSeq_set_attribute_values_mapping(igraphmodule_VertexSeqObject* self, PyObject* attrname, PyObject* values) {
  PyObject *dict, *list, *item;
  igraphmodule_GraphObject *gr;
  long i, n;

  gr=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!gr) return -1;

  dict = ((PyObject**)gr->g.attr)[ATTRHASH_IDX_VERTEX];
  
  if (values == 0) return PyDict_DelItem(dict, attrname);

  n=PySequence_Size(values);
  if (n<0) return -1;
  if (n != (long)igraph_vcount(&gr->g)) {
    PyErr_SetString(PyExc_ValueError, "value list length must be equal to the number of vertices in the graph");
    return -1;
  }

  /* Check if we already have attributes with the given name */
  list = PyDict_GetItem(dict, attrname);
  if (list != 0) {
    /* Yes, we have. Modify its items to the items found in values */
    for (i=0; i<n; i++) {
      item = PyList_GetItem(values, i);
      if (item == 0) return -1;
      Py_INCREF(item);
      if (PyList_SetItem(list, i, item)) {
        Py_DECREF(item);
        return -1;
      } /* PyList_SetItem stole a reference to the item automatically */ 
    }
  } else if (values != 0) {
    /* We don't have attributes with the given name yet. Create an entry
     * in the dict, create a new list and copy everything */
    list = PyList_New(n);
    if (list == 0) return -1;
    for (i=0; i<n; i++) {
      item = PyList_GetItem(values, i);
      if (item == 0) {
        Py_DECREF(list);
        return -1;
      }
      Py_INCREF(item);
      PyList_SET_ITEM(list, i, item);
    }
    if (PyDict_SetItem(dict, attrname, list)) {
      Py_DECREF(list);
      return -1;
    }
  }

  return 0;
}

PyObject* igraphmodule_VertexSeq_set_attribute_values(igraphmodule_VertexSeqObject *self, PyObject *args, PyObject *kwds) {
  static char* kwlist[] = { "attrname", "values", NULL };
  PyObject *attrname, *values;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist,
				   &attrname, &values))
	return NULL;

  if (igraphmodule_VertexSeq_set_attribute_values_mapping(self, attrname, values))
	return NULL;

  Py_RETURN_NONE;
}

/**
 * \ingroup python_interface_vertexseq
 * Method table for the \c igraph.VertexSeq object
 */
PyMethodDef igraphmodule_VertexSeq_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_VertexSeq_attributes,
   METH_NOARGS,
   "attributes() -> list\n\n"
   "Returns the attribute name list of the graph's vertices\n"
  },
  {"get_attribute_values", (PyCFunction)igraphmodule_VertexSeq_get_attribute_values,
   METH_O,
   "get_attribute_values(attrname) -> list\n"
   "Returns the value of a given vertex attribute for all vertices in a list.\n\n"
   "The values stored in the list are exactly the same objects that are stored\n"
   "in the vertex attribute, meaning that in the case of mutable objects,\n"
   "the modification of the list element does affect the attribute stored in\n"
   "the vertex. In the case of immutable objects, modification of the list\n"
   "does not affect the attribute values.\n\n"
   "@param attrname: the name of the attribute\n"
  },
  {"set_attribute_values", (PyCFunction)igraphmodule_VertexSeq_set_attribute_values,
   METH_VARARGS | METH_KEYWORDS,
   "set_attribute_values(attrname, values) -> list\n"
   "Sets the value of a given vertex attribute for all vertices\n"
   "@param attrname: the name of the attribute\n"
   "@param values: the new attribute values in a list\n"
  },
  {NULL}
};

/**
 * \ingroup python_interface_vertexseq
 * This is the collection of functions necessary to implement the
 * vertex sequence as a real sequence (e.g. allowing to reference
 * vertices by indices)
 */
static PySequenceMethods igraphmodule_VertexSeq_as_sequence = {
  (inquiry)igraphmodule_VertexSeq_sq_length,
  0,               /* sq_concat */
  0,               /* sq_repeat */
  (intargfunc)igraphmodule_VertexSeq_sq_item, /* sq_item */
  0,                                          /* sq_slice */
  0,                                          /* sq_ass_item */
  0,                                          /* sq_ass_slice */
  0,                                          /* sq_contains */
  0,                                          /* sq_inplace_concat */
  0,                                          /* sq_inplace_repeat */
};

/**
 * \ingroup python_interface_vertexseq
 * This is the collection of functions necessary to implement the
 * vertex sequence as a mapping (which maps attribute names to values)
 */
static PyMappingMethods igraphmodule_VertexSeq_as_mapping = {
	/* returns the number of vertex attributes */
	(inquiry) igraphmodule_VertexSeq_attribute_count,
	/* returns the values of an attribute by name */
	(binaryfunc) igraphmodule_VertexSeq_get_attribute_values_mapping,
	/* sets the values of an attribute by name */
	(objobjargproc) igraphmodule_VertexSeq_set_attribute_values_mapping,
};

/** \ingroup python_interface_vertexseq
 * Python type object referencing the methods Python calls when it performs various operations on
 * a vertex sequence of a graph
 */
PyTypeObject igraphmodule_VertexSeqType =
{
  PyObject_HEAD_INIT(NULL)                    /* */
  0,                                          /* ob_size */
  "igraph.VertexSeq",                         /* tp_name */
  sizeof(igraphmodule_VertexSeqObject),       /* tp_basicsize */
  0,                                          /* tp_itemsize */
  (destructor)igraphmodule_VertexSeq_dealloc, /* tp_dealloc */
  0,                                          /* tp_print */
  0,                                          /* tp_getattr */
  0,                                          /* tp_setattr */
  0,                                          /* tp_compare */
  0,                                          /* tp_repr */
  0,                                          /* tp_as_number */
  &igraphmodule_VertexSeq_as_sequence,        /* tp_as_sequence */
  &igraphmodule_VertexSeq_as_mapping,         /* tp_as_mapping */
  0,                                          /* tp_hash */
  0,                                          /* tp_call */
  0,                                          /* tp_str */
  0,                                          /* tp_getattro */
  0,                                          /* tp_setattro */
  0,                                          /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /* tp_flags */
  "Class representing the vertex set of a graph.\n\n"
  "The individual vertices can be accessed by indexing the vertex sequence\n"
  "object. It can be used as an iterable as well, or even in a list\n"
  "comprehension:\n\n"
  "  >>> g=igraph.Graph.Full(3)\n"
  "  >>> for v in g.vs:\n"
  "  ...   v[\"value\"] = v.index ** 2\n"
  "  ...\n"
  "  >>> [v[\"value\"] ** 0.5 for v in g.vs]\n"
  "  [0.0, 1.0, 2.0]\n\n"
  "The vertex set can also be used as a dictionary where the keys are the\n"
  "attribute names. The values corresponding to the keys are the values\n"
  "of the given attribute of every vertex in the graph:\n\n"
  "  >>> g=igraph.Graph.Full(3)\n"
  "  >>> for idx, v in enumerate(g.vs):\n"
  "  ...   v[\"weight\"] = idx*(idx+1)\n"
  "  ...\n"
  "  >>> g.vs[\"weight\"]\n"
  "  [0, 2, 6]\n"
  "  >>> g.vs[\"weight\"] = range(3)\n"
  "  >>> g.vs[\"weight\"]\n"
  "  [0, 1, 2]\n", /* tp_doc */
  0,                                          /* tp_traverse */
  0,                                          /* tp_clear */
  0,                                          /* tp_richcompare */
  0,                                          /* tp_weaklistoffset */
  0,                                          /* tp_iter */
  0,                                          /* tp_iternext */
  igraphmodule_VertexSeq_methods,             /* tp_methods */
};

