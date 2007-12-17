/* -*- mode: C -*-  */
/* vim:ts=2 sts=2 sw=2 et: */

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
  igraph_vs_all(&o->vs);
  
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

  igraph_vs_destroy(&self->vs);

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
 * \brief Returns the length of the sequence
 */
int igraphmodule_VertexSeq_sq_length(igraphmodule_VertexSeqObject* self) {
  igraph_t *g;
  igraph_integer_t result;
  g=&((igraphmodule_GraphObject*)PyWeakref_GetObject(self->gref))->g;
  if (igraph_vs_size(g, &self->vs, &result)) return -1;
  return (int)result;
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Returns the item at the given index in the sequence
 */
PyObject* igraphmodule_VertexSeq_sq_item(igraphmodule_VertexSeqObject* self,
                     int i) {
  igraphmodule_GraphObject *o;
  igraph_t *g;
  long idx = -1;

  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  g=&o->g;
  switch (igraph_vs_type(&self->vs)) {
    case IGRAPH_VS_ALL:
      if (i >= 0 && i < (long)igraph_vcount(g)) idx = i;
      break;
    case IGRAPH_VS_VECTOR:
    case IGRAPH_VS_VECTORPTR:
      if (i >= 0 && i < igraph_vector_size(self->vs.data.vecptr))
        idx = (long)VECTOR(*self->vs.data.vecptr)[i];
      break;
    case IGRAPH_VS_1:
      if (i == 0) idx = (long)self->vs.data.vid;
      break;
    case IGRAPH_VS_SEQ:
      if (i >= 0 && i < self->vs.data.seq.to - self->vs.data.seq.from)
        idx = (long)(self->vs.data.seq.from + i);
      break;
    /* TODO: IGRAPH_VS_ADJ, IGRAPH_VS_NONADJ - someday :) */
  }

  if (idx < 0) {
    PyErr_SetString(PyExc_IndexError, "vertex index out of range");
    return NULL;
  }

  return igraphmodule_Vertex_New(self->gref, idx);
}

/** \ingroup python_interface_vertexseq
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_VertexSeq_attribute_names(igraphmodule_VertexSeqObject* self) {
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
  PyObject *result=0, *values, *item;
  long int i, n;

  gr=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!gr) return NULL;
  
  values=PyDict_GetItem(((PyObject**)gr->g.attr)[ATTRHASH_IDX_VERTEX], o);
  if (!values) {
    PyErr_SetString(PyExc_KeyError, "Attribute does not exist");
    return NULL;
  } else if (PyErr_Occurred()) return NULL;

  switch (igraph_vs_type(&self->vs)) {
    case IGRAPH_VS_NONE:
      n = 0;
      result = PyList_New(0);
      break;

    case IGRAPH_VS_ALL:
      n = PyList_Size(values);
      result = PyList_New(n);
      if (!result) return 0;
      
      for (i=0; i<n; i++) {
        item = PyList_GET_ITEM(values, i);
        Py_INCREF(item);
        PyList_SET_ITEM(result, i, item);
      }
      break;

    case IGRAPH_VS_VECTOR:
    case IGRAPH_VS_VECTORPTR:
      n = igraph_vector_size(self->vs.data.vecptr);
      result = PyList_New(n);
      if (!result) return 0;

      for (i=0; i<n; i++) {
        item = PyList_GET_ITEM(values, (long)VECTOR(*self->vs.data.vecptr)[i]);
        Py_INCREF(item);
        PyList_SET_ITEM(result, i, item);
      }
      break;

    case IGRAPH_VS_SEQ:
      n = self->vs.data.seq.to - self->vs.data.seq.from;
      result = PyList_New(n);
      if (!result) return 0;

      for (i=0; i<n; i++) {
        item = PyList_GET_ITEM(values, (long)self->vs.data.seq.from+i);
        Py_INCREF(item);
        PyList_SET_ITEM(result, i, item);
      }
      break;

    default:
      PyErr_SetString(PyExc_RuntimeError, "invalid vertex selector");
  }

  return result;
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
  igraph_vector_t vs;
  long i, n;

  gr=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!gr) return -1;

  dict = ((PyObject**)gr->g.attr)[ATTRHASH_IDX_VERTEX];
  
  if (values == 0) {
    if (igraph_vs_type(&self->vs) == IGRAPH_VS_ALL)
      return PyDict_DelItem(dict, attrname);
    PyErr_SetString(PyExc_TypeError, "can't delete attribute from a vertex sequence not representing the whole graph");
    return -1;
  }

  n=PySequence_Size(values);
  if (n<0) return -1;

  if (igraph_vs_type(&self->vs) == IGRAPH_VS_ALL) {
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
        item = PyList_GET_ITEM(values, i);
        if (item == 0) { Py_DECREF(list); return -1; }
        Py_INCREF(item);
        PyList_SET_ITEM(list, i, item);
      }
      if (PyDict_SetItem(dict, attrname, list)) {
        Py_DECREF(list);
        return -1;
      }
    }
  } else {
    /* We are working with a subset of the graph. Convert the sequence to a
     * vector and loop through it */
    if (igraph_vector_init(&vs, 0)) {
      igraphmodule_handle_igraph_error();
      return -1;
    } 
    if (igraph_vs_as_vector(&gr->g, self->vs, &vs)) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(&vs);
      return -1;
    }
    if (n != (long)igraph_vector_size(&vs)) {
      PyErr_SetString(PyExc_ValueError, "value list length must be equal to the number of vertices in the vertex set");
      return -1;
    }
    /* Check if we already have attributes with the given name */
    list = PyDict_GetItem(dict, attrname);
    if (list != 0) {
      /* Yes, we have. Modify its items to the items found in values */
      for (i=0; i<n; i++) {
        item = PyList_GetItem(values, VECTOR(vs)[i]);
        if (item == 0) continue;
        Py_INCREF(item);
        if (PyList_SetItem(list, i, item)) {
          Py_DECREF(item);
          return -1;
        } /* PyList_SetItem stole a reference to the item automatically */ 
      }
    } else if (values != 0) {
      /* We don't have attributes with the given name yet. Create an entry
       * in the dict, create a new list, fill with None for vertices not in the
       * sequence and copy the rest */
      long n2 = igraph_vcount(&gr->g);
      list = PyList_New(n2);
      if (list == 0) return -1;
      for (i=0; i<n2; i++) {
        Py_INCREF(Py_None);
        PyList_SET_ITEM(list, i, Py_None);
      }
      for (i=0; i<n; i++) {
        item = PyList_GET_ITEM(values, (long)VECTOR(vs)[i]);
        if (item == 0) { Py_DECREF(list); return -1; }
        Py_INCREF(item);
        PyList_SET_ITEM(list, i, item);
      }
      if (PyDict_SetItem(dict, attrname, list)) {
        Py_DECREF(list);
        return -1;
      }
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
 * \brief Selects a subset of the vertex sequence based on some criteria
 */
PyObject* igraphmodule_VertexSeq_select(igraphmodule_VertexSeqObject *self,
  PyObject *args, PyObject *kwds) {
  igraphmodule_VertexSeqObject *result;
  igraphmodule_GraphObject *gr;
  long i, j, n, m;
  gr=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  result=(igraphmodule_VertexSeqObject*)igraphmodule_VertexSeq_New(gr);
  if (result==0) return NULL;

  result->vs = self->vs;
  /* First, filter by positional arguments */
  n = PyTuple_Size(args);
  for (i=0; i<n; i++) {
    PyObject *item = PyTuple_GET_ITEM(args, i);
    if (item == Py_None) {
      /* None means: select nothing */
      igraph_vs_destroy(&result->vs);
      igraph_vs_none(&result->vs);
      /* We can simply bail out here */
      return (PyObject*)result;
    } else if (PyCallable_Check(item)) {
      /* Call the callable for every vertex in the current sequence to
       * determine what's up */
      igraph_bool_t was_excluded = 0;
      igraph_vector_t v;

      if (igraph_vector_init(&v, 0)) {
        igraphmodule_handle_igraph_error();
        return 0;
      }

      m = PySequence_Size((PyObject*)result);
      for (j=0; j<m; j++) {
        PyObject *vertex = PySequence_GetItem((PyObject*)result, j);
        PyObject *call_result;
        if (vertex == 0) {
          Py_DECREF(result);
          igraph_vector_destroy(&v);
          return NULL;
        }
        call_result = PyObject_CallFunctionObjArgs(item, vertex, NULL);
        if (call_result == 0) {
          Py_DECREF(vertex); Py_DECREF(result);
          igraph_vector_destroy(&v);
          return NULL;
        }
        if (PyObject_IsTrue(call_result))
          igraph_vector_push_back(&v,
            igraphmodule_Vertex_get_index_long((igraphmodule_VertexObject*)vertex));
        else was_excluded=1;
        Py_DECREF(call_result);
        Py_DECREF(vertex);
      }

      if (was_excluded) {
        igraph_vs_destroy(&result->vs);
        if (igraph_vs_vector_copy(&result->vs, &v)) {
          Py_DECREF(result);
          igraph_vector_destroy(&v);
          igraphmodule_handle_igraph_error();
          return NULL;
        }
      }

      igraph_vector_destroy(&v);
    } else if (PyInt_Check(item)) {
      /* Integers are treated specially: from now on, all remaining items
       * in the argument list must be integers and they will be used together
       * to restrict the vertex set. Integers are interpreted as indices on the
       * vertex set and NOT on the original, untouched vertex sequence of the
       * graph */
      igraph_vector_t v, v2;
      if (igraph_vector_init(&v, 0)) {
        igraphmodule_handle_igraph_error();
        return 0;
      }
      if (igraph_vector_init(&v2, 0)) {
        igraph_vector_destroy(&v);
        igraphmodule_handle_igraph_error();
        return 0;
      }
      if (igraph_vs_as_vector(&gr->g, result->vs, &v2)) {
        igraph_vector_destroy(&v);
        igraph_vector_destroy(&v2);
        igraphmodule_handle_igraph_error();
        return 0;
      }
      for (; i<n; i++) {
        PyObject *item2 = PyTuple_GET_ITEM(args, i);
        long idx;
        if (!PyInt_Check(item2)) {
          Py_DECREF(result);
          PyErr_SetString(PyExc_TypeError, "vertex indices expected");
          igraph_vector_destroy(&v);
          igraph_vector_destroy(&v2);
          return NULL;
        }
        idx = PyInt_AsLong(item2);
        if (igraph_vector_push_back(&v, VECTOR(v2)[idx])) {
          Py_DECREF(result);
          igraphmodule_handle_igraph_error();
          igraph_vector_destroy(&v);
          igraph_vector_destroy(&v2);
          return NULL;
        }
      }
      igraph_vector_destroy(&v2);
      igraph_vs_destroy(&result->vs);
      if (igraph_vs_vector_copy(&result->vs, &v)) {
        Py_DECREF(result);
        igraphmodule_handle_igraph_error();
        igraph_vector_destroy(&v);
        return NULL;
      }
      igraph_vector_destroy(&v);
    } else {
      /* Iterators and everything that was not handled directly */
      PyObject *iter, *item2;
      igraph_vector_t v, v2;
      
      iter = PyObject_GetIter(item);
      if (iter == 0) {
        PyErr_SetString(PyExc_TypeError, "invalid vertex filter among positional arguments");
        Py_DECREF(result);
        return 0;
      }
      /* Allocate stuff */
      if (igraph_vector_init(&v, 0)) {
        Py_DECREF(iter);
        igraphmodule_handle_igraph_error();
        return 0;
      }
      if (igraph_vector_init(&v2, 0)) {
        Py_DECREF(iter);
        igraph_vector_destroy(&v);
        igraphmodule_handle_igraph_error();
        return 0;
      }
      if (igraph_vs_as_vector(&gr->g, result->vs, &v2)) {
        Py_DECREF(iter);
        igraph_vector_destroy(&v);
        igraph_vector_destroy(&v2);
        igraphmodule_handle_igraph_error();
        return 0;
      }
      /* Do the iteration */
      while ((item2=PyIter_Next(iter)) != 0) {
        if (PyInt_Check(item2)) {
          long idx = PyInt_AsLong(item2);
          Py_DECREF(item2);
          if (igraph_vector_push_back(&v, VECTOR(v2)[idx])) {
            Py_DECREF(result);
            Py_DECREF(iter);
            igraphmodule_handle_igraph_error();
            igraph_vector_destroy(&v);
            igraph_vector_destroy(&v2);
            return NULL;
          }
        } else {
          /* We simply ignore elements that we don't know */
          Py_DECREF(item2);
        }
      }
      /* Deallocate stuff */
      igraph_vector_destroy(&v2);
      Py_DECREF(iter);
      if (PyErr_Occurred()) {
        igraph_vector_destroy(&v);
        Py_DECREF(result);
        return 0;
      }
      igraph_vs_destroy(&result->vs);
      if (igraph_vs_vector_copy(&result->vs, &v)) {
        Py_DECREF(result);
        igraphmodule_handle_igraph_error();
        igraph_vector_destroy(&v);
        return NULL;
      }
      igraph_vector_destroy(&v);
    }
  }

  return (PyObject*)result;
}

/**
 * \ingroup python_interface_vertexseq
 * Method table for the \c igraph.VertexSeq object
 */
PyMethodDef igraphmodule_VertexSeq_methods[] = {
  {"attribute_names", (PyCFunction)igraphmodule_VertexSeq_attribute_names,
   METH_NOARGS,
   "attribute_names() -> list\n\n"
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
   "Sets the value of a given vertex attribute for all vertices\n\n"
   "@param attrname: the name of the attribute\n"
   "@param values: the new attribute values in a list\n"
  },
  {"select", (PyCFunction)igraphmodule_VertexSeq_select,
   METH_VARARGS | METH_KEYWORDS,
   "select(...) -> VertexSeq\n"
   "Selects a subset of the vertex sequence based on some criteria\n\n"
   "The selection criteria can be specified by the positional and the keyword\n"
   "arguments.\n\n"
   "  - If the first positional argument is C{None}, an empty sequence is\n"
   "    returned.\n\n"
   "  - If the first positional argument is a callable object, the object\n"
   "    will be called for every vertex in the sequence. If it returns\n"
   "    C{True}, the vertex will be included, otherwise it will be excluded.\n\n"
   "  - If the first positional argument is an iterable, it must return\n"
   "    integers and they will be considered as indices of the current vertex\n"
   "    set (NOT the vertex set of the graph -- the difference matters when\n"
   "    one filters a vertex set that has already been filtered by a previous\n"
   "    invocation of L{VertexSet.select()}. In this case, the indices do not\n"
   "    refer directly to the vertices of the graph but to the elements of the\n"
   "    vertex sequence.\n\n"
   "  - If the first positional argument is an integer, all remaining arguments\n"
   "    are expected to be integers. They are considered as indices of the\n"
   "    current vertex set again.\n\n"
   "\n"
   "@return: the new, filtered vertex sequence\n"
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
  "Class representing a sequence of vertices in the graph.\n\n"
  "This class is most easily accessed by the C{vs} field of the L{Graph} object\n"
  "which returns an ordered sequence of all vertices in the graph. The vertex\n"
  "sequence can be refined by invoking the L{VertexSeq.select()} method.\n\n"
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

