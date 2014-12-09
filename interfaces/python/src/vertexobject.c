/* -*- mode: C -*-  */
/* vim: set ts=2 sw=2 sts=2 et: */

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

#include "attributes.h"
#include "convert.h"
#include "error.h"
#include "graphobject.h"
#include "vertexobject.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_vertex Vertex object
 */

PyTypeObject igraphmodule_VertexType;

/**
 * \ingroup python_interface_vertex
 * \brief Checks whether the given Python object is a vertex
 */
int igraphmodule_Vertex_Check(PyObject* obj) {
  if (!obj)
    return 0;
  
  return PyObject_IsInstance(obj, (PyObject*)(&igraphmodule_VertexType));
}

/**
 * \ingroup python_interface_vertex
 * \brief Checks whether the index in the given vertex object is a valid one.
 * \return nonzero if the vertex object is valid. Raises an appropriate Python
 *   exception and returns zero if the vertex object is invalid.
 */
int igraphmodule_Vertex_Validate(PyObject* obj) {
  igraph_integer_t n;
  igraphmodule_VertexObject *self;
  igraphmodule_GraphObject *graph;

  if (!igraphmodule_Vertex_Check(obj)) {
    PyErr_SetString(PyExc_TypeError, "object is not a Vertex");
    return 0;
  }

  self = (igraphmodule_VertexObject*)obj;
  graph = self->gref;

  if (graph == 0) {
    PyErr_SetString(PyExc_ValueError, "Vertex object refers to a null graph");
    return 0;
  }

  if (self->idx < 0) {
    PyErr_SetString(PyExc_ValueError, "Vertex object refers to a negative vertex index");
    return 0;
  }

  n = igraph_vcount(&graph->g);

  if (self->idx >= n) {
    PyErr_SetString(PyExc_ValueError, "Vertex object refers to a nonexistent vertex");
    return 0;
  }

  return 1;
}

/**
 * \ingroup python_interface_vertex
 * \brief Allocates a new Python vertex object
 * \param gref the \c igraph.Graph being referenced by the vertex
 * \param idx the index of the vertex
 * 
 * \warning \c igraph references its vertices by indices, so if
 * you delete some vertices from the graph, the vertex indices will
 * change. Since the \c igraph.Vertex objects do not follow these
 * changes, your existing vertex objects will point to elsewhere
 * (or they might even get invalidated).
 */
PyObject* igraphmodule_Vertex_New(igraphmodule_GraphObject *gref, igraph_integer_t idx) {
  igraphmodule_VertexObject* self;
  self=PyObject_New(igraphmodule_VertexObject, &igraphmodule_VertexType);
  if (self) {
    RC_ALLOC("Vertex", self);
    Py_INCREF(gref);
    self->gref=gref;
    self->idx=idx;
    self->hash=-1;
  }
  return (PyObject*)self;
}

/**
 * \ingroup python_interface_vertex
 * \brief Clears the vertex's subobject (before deallocation)
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

  PyObject_Del((PyObject*)self);
}

/** \ingroup python_interface_vertex
 * \brief Formats an \c igraph.Vertex object to a string
 * 
 * \return the formatted textual representation as a \c PyObject
 */
PyObject* igraphmodule_Vertex_repr(igraphmodule_VertexObject *self) {
  PyObject *s;
  PyObject *attrs;
#ifndef IGRAPH_PYTHON3
  PyObject *grepr, *drepr;
#endif

  attrs = igraphmodule_Vertex_attributes(self);
  if (attrs == 0)
    return NULL;

#ifdef IGRAPH_PYTHON3
  s = PyUnicode_FromFormat("igraph.Vertex(%R, %ld, %R)",
      (PyObject*)self->gref, (long int)self->idx, attrs);
  Py_DECREF(attrs);
#else
  grepr=PyObject_Repr((PyObject*)self->gref);
  drepr=PyObject_Repr(igraphmodule_Vertex_attributes(self));
  Py_DECREF(attrs);
  if (!grepr || !drepr) {
    Py_XDECREF(grepr);
    Py_XDECREF(drepr);
    return NULL;
  }
  s=PyString_FromFormat("igraph.Vertex(%s,%ld,%s)", PyString_AsString(grepr),
    (long int)self->idx, PyString_AsString(drepr));
  Py_DECREF(grepr);
  Py_DECREF(drepr);
#endif
  return s;
}

/** \ingroup python_interface_vertex
 * \brief Returns the hash code of the vertex
 */
Py_hash_t igraphmodule_Vertex_hash(igraphmodule_VertexObject* self) {
  Py_hash_t hash_graph;
  Py_hash_t hash_index;
  Py_hash_t result;
  PyObject* index_o;

  if (self->hash != -1)
    return self->hash;

  index_o = PyInt_FromLong((long int)self->idx);
  if (index_o == 0)
    return -1;

  hash_index = PyObject_Hash(index_o);
  Py_DECREF(index_o);

  if (hash_index == -1)
    return -1;

  hash_graph = PyObject_Hash((PyObject*)self->gref);
  if (hash_graph == -1)
    return -1;

  result = hash_graph ^ hash_index;
  if (result == -1)
    result = 590923713U;

  self->hash = result;

  return result;
}

/** \ingroup python_interface_vertex
 * \brief Rich comparison of a vertex with another
 */
PyObject* igraphmodule_Vertex_richcompare(igraphmodule_VertexObject *a,
    PyObject *b, int op) {

  igraphmodule_VertexObject* self = a;
  igraphmodule_VertexObject* other;

  if (!igraphmodule_Vertex_Check(b))
    Py_RETURN_NOTIMPLEMENTED;

  other = (igraphmodule_VertexObject*)b;

  if (self->gref != other->gref)
    Py_RETURN_FALSE;

  switch (op) {
    case Py_EQ:
      Py_RETURN(self->idx == other->idx);

    case Py_NE:
      Py_RETURN(self->idx != other->idx);

    case Py_LE:
      Py_RETURN(self->idx <= other->idx);

    case Py_LT:
      Py_RETURN(self->idx <  other->idx);

    case Py_GE:
      Py_RETURN(self->idx >= other->idx);

    case Py_GT:
      Py_RETURN(self->idx >  other->idx);

    default:
      Py_RETURN_NOTIMPLEMENTED;
  }
}

/** \ingroup python_interface_vertex
 * \brief Returns the number of vertex attributes
 */
Py_ssize_t igraphmodule_Vertex_attribute_count(igraphmodule_VertexObject* self) {
  igraphmodule_GraphObject *o = self->gref;
  
  if (!o) return 0;
  if (!((PyObject**)o->g.attr)[1]) return 0;
  return PyDict_Size(((PyObject**)o->g.attr)[1]);
}

/** \ingroup python_interface_vertex
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_Vertex_attribute_names(igraphmodule_VertexObject* self) {
  if (!self->gref) return NULL;
  return igraphmodule_Graph_vertex_attributes(self->gref);
}

/** \ingroup python_interface_vertex
 * \brief Returns a dict with attribue names and values
 */
PyObject* igraphmodule_Vertex_attributes(igraphmodule_VertexObject* self) {
  igraphmodule_GraphObject *o = self->gref;
  PyObject *names, *dict;
  long i, n;

  if (!igraphmodule_Vertex_Validate((PyObject*)self))
    return 0;

  dict=PyDict_New();
  if (!dict) return NULL;

  names=igraphmodule_Graph_vertex_attributes(o);
  if (!names) {
    Py_DECREF(dict);
    return NULL;
  }

  n=PyList_Size(names);
  for (i=0; i<n; i++) {
    PyObject *name = PyList_GetItem(names, i);
    if (name) {
      PyObject *dictit;
      dictit = PyDict_GetItem(((PyObject**)o->g.attr)[ATTRHASH_IDX_VERTEX], name);
      if (dictit) {
        PyObject *value = PyList_GetItem(dictit, self->idx);
        if (value) {
          /* No need to Py_INCREF, PyDict_SetItem will do that */
          PyDict_SetItem(dict, name, value);
        }
      }
    }
  }

  Py_DECREF(names);

  return dict;
}

/**
 * \ingroup python_interface_vertex
 * \brief Updates some attributes of a vertex
 * 
 * Incidentally, this method is re-used intact in edgeobject.c for edges.
 *
 * \param self the vertex object
 * \param args positional arguments
 * \param kwds keyword arguments
 */
PyObject* igraphmodule_Vertex_update_attributes(PyObject* self, PyObject* args,
    PyObject* kwds) {
  PyObject* items[] = { Py_None, kwds, 0 };
  PyObject** pObj;
  PyObject *key, *value, *it, *item, *keys;

  igraph_bool_t ok = 1;

  if (!PyArg_ParseTuple(args, "|O", &items[0]))
    return NULL;

  pObj = items;
  for (pObj = items; ok && *pObj != 0; pObj++) {
    PyObject* obj = *pObj;
    PyObject* keys_func;

    if (obj == Py_None)
      continue;

    keys_func = PyObject_GetAttrString(obj, "keys");
    if (keys_func == 0)
      PyErr_Clear();

    if (keys_func != 0 && PyCallable_Check(keys_func)) {
      /* Object has a "keys" method, so we iterate over the keys */
      keys = PyObject_CallObject(keys_func, 0);
      if (keys == 0) {
        ok = 0;
      } else {
        /* Iterate over the keys */
        it = PyObject_GetIter(keys);
        if (it == 0) {
          ok = 0;
        } else {
          while (ok && ((key = PyIter_Next(it)) != 0)) {
            value = PyObject_GetItem(obj, key);
            if (value == 0) {
              ok = 0;
            } else {
              PyObject_SetItem((PyObject*)self, key, value);
              Py_DECREF(value);
            }
            Py_DECREF(key);
          }
          Py_DECREF(it);
          if (PyErr_Occurred())
            ok = 0;
        }
        Py_DECREF(keys);
      }
    } else {
      /* Object does not have a "keys" method; assume that it
       * yields tuples when treated as an iterator */
      it = PyObject_GetIter(obj);
      if (!it) {
        ok = 0;
      } else {
        while (ok && ((item = PyIter_Next(it)) != 0)) {
          if (!PySequence_Check(item) || PyBaseString_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "cannot convert update sequence element to a sequence");
            ok = 0;
          } else {
            key = PySequence_GetItem(item, 0);
            if (key == 0) {
              ok = 0;
            } else {
              value = PySequence_GetItem(item, 1);
              if (value == 0) {
                ok = 0;
              } else {
                PyObject_SetItem((PyObject*)self, key, value);
                Py_DECREF(value);
              }
              Py_DECREF(key);
            }
          }
          Py_DECREF(item);
        }
        Py_DECREF(it);
        if (PyErr_Occurred())
          ok = 0;
      }
    }

    if (keys_func != 0) {
      Py_DECREF(keys_func);
    }
  }

  if (ok)
    Py_RETURN_NONE;
  return 0;
}

/** \ingroup python_interface_vertex
 * \brief Returns the corresponding value to a given attribute of the vertex
 * \param self the vertex object
 * \param s the attribute name to be queried
 */
PyObject* igraphmodule_Vertex_get_attribute(igraphmodule_VertexObject* self,
                       PyObject* s) {
  igraphmodule_GraphObject *o = self->gref;
  PyObject* result;

  if (!igraphmodule_Vertex_Validate((PyObject*)self))
    return 0;

  if (!igraphmodule_attribute_name_check(s))
    return 0;

  result=PyDict_GetItem(((PyObject**)o->g.attr)[ATTRHASH_IDX_VERTEX], s);
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
  igraphmodule_GraphObject *o=self->gref;
  PyObject* result;
  int r;
  
  if (!igraphmodule_Vertex_Validate((PyObject*)self))
    return -1;

  if (!igraphmodule_attribute_name_check(k))
    return -1;

  if (PyString_IsEqualToASCIIString(k, "name"))
    igraphmodule_invalidate_vertex_name_index(&o->g);

  if (v==NULL)
    // we are deleting attribute
    return PyDict_DelItem(((PyObject**)o->g.attr)[ATTRHASH_IDX_VERTEX], k);
  
  result=PyDict_GetItem(((PyObject**)o->g.attr)[ATTRHASH_IDX_VERTEX], k);
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
    Py_DECREF(result);   /* compensating for PyDict_SetItem */
    return 0;
  }
  
  return -1;
}

/**
 * \ingroup python_interface_vertex
 * Returns the vertex index
 */
PyObject* igraphmodule_Vertex_get_index(igraphmodule_VertexObject* self, void* closure) {
  return PyInt_FromLong((long int)self->idx);
}

/**
 * \ingroup python_interface_vertex
 * Returns the vertex index as an igraph_integer_t
 */
igraph_integer_t igraphmodule_Vertex_get_index_igraph_integer(igraphmodule_VertexObject* self) {
  return self->idx;
}

/**
 * \ingroup python_interface_vertex
 * Returns the vertex index as an ordinary C long
 */
long igraphmodule_Vertex_get_index_long(igraphmodule_VertexObject* self) {
  return (long)self->idx;
}

/**
 * \ingroup python_interface_vertexseq
 * Returns the graph where the vertex belongs
 */
PyObject* igraphmodule_Vertex_get_graph(igraphmodule_VertexObject* self,
  void* closure) {
  Py_INCREF(self->gref);
  return (PyObject*)self->gref;
}

/**************************************************************************/
/* Implementing proxy method in Vertex that just forward the call to the
 * appropriate Graph method.
 * 
 * These methods may also execute a postprocessing function on the result
 * of the Graph method; for instance, this mechanism is used to turn the
 * result of Graph.neighbors() (which is a list of vertex indices) into a
 * list of Vertex objects.
 */

/* Dummy postprocessing function that does nothing. */
static PyObject* _identity(igraphmodule_VertexObject* vertex, PyObject* obj) {
  Py_INCREF(obj);
  return obj;
}

/* Postprocessing function that converts a Python list of integers into a
 * list of vertices in-place. */
static PyObject* _convert_to_vertex_list(igraphmodule_VertexObject* vertex, PyObject* obj) {
  Py_ssize_t i, n;

  if (!PyList_Check(obj)) {
    PyErr_SetString(PyExc_TypeError, "_convert_to_vertex_list expected list of integers");
    return NULL;
  }

  n = PyList_Size(obj);
  for (i = 0; i < n; i++) {
    PyObject* idx = PyList_GET_ITEM(obj, i);
    PyObject* v;
    int idx_int;

    if (!PyInt_Check(idx)) {
      PyErr_SetString(PyExc_TypeError, "_convert_to_vertex_list expected list of integers");
      return NULL;
    }

    if (PyInt_AsInt(idx, &idx_int))
      return NULL;

    v = igraphmodule_Vertex_New(vertex->gref, idx_int);
    PyList_SetItem(obj, i, v);   /* reference to v stolen, reference to idx discarded */
  }

  Py_INCREF(obj);
  return obj;
}

#define GRAPH_PROXY_METHOD_PP(FUNC, METHODNAME, POSTPROCESS) \
    PyObject* igraphmodule_Vertex_##FUNC(igraphmodule_VertexObject* self, PyObject* args, PyObject* kwds) { \
      PyObject *new_args, *item, *result;                     \
      long int i, num_args = args ? PyTuple_Size(args)+1 : 1; \
                                                              \
      /* Prepend ourselves to args */                         \
      new_args = PyTuple_New(num_args);                       \
      Py_INCREF(self); PyTuple_SET_ITEM(new_args, 0, (PyObject*)self);   \
      for (i = 1; i < num_args; i++) {                        \
        item = PyTuple_GET_ITEM(args, i-1);                   \
        Py_INCREF(item); PyTuple_SET_ITEM(new_args, i, item); \
      }                                                       \
                                                              \
      /* Get the method instance */                           \
      item = PyObject_GetAttrString((PyObject*)(self->gref), METHODNAME);  \
      result = PyObject_Call(item, new_args, kwds);           \
      Py_DECREF(item);                                        \
      Py_DECREF(new_args);                                    \
                                                              \
      /* Optional postprocessing */                           \
      if (result) {                                           \
        PyObject* pp_result = POSTPROCESS(self, result);      \
        Py_DECREF(result);                                    \
        return pp_result;                                     \
      }                                                       \
      return NULL;                                            \
    }

#define GRAPH_PROXY_METHOD(FUNC, METHODNAME) \
        GRAPH_PROXY_METHOD_PP(FUNC, METHODNAME, _identity)

GRAPH_PROXY_METHOD(betweenness, "betweenness");
GRAPH_PROXY_METHOD(closeness, "closeness");
GRAPH_PROXY_METHOD(constraint, "constraint");
GRAPH_PROXY_METHOD(degree, "degree");
GRAPH_PROXY_METHOD(delete, "delete_vertices");
GRAPH_PROXY_METHOD(diversity, "diversity");
GRAPH_PROXY_METHOD(eccentricity, "eccentricity");
GRAPH_PROXY_METHOD(get_shortest_paths, "get_shortest_paths");
GRAPH_PROXY_METHOD(indegree, "indegree");
GRAPH_PROXY_METHOD(is_minimal_separator, "is_minimal_separator");
GRAPH_PROXY_METHOD(is_separator, "is_separator");
GRAPH_PROXY_METHOD_PP(neighbors, "neighbors", _convert_to_vertex_list);
GRAPH_PROXY_METHOD(outdegree, "outdegree");
GRAPH_PROXY_METHOD(pagerank, "pagerank");
GRAPH_PROXY_METHOD_PP(predecessors, "predecessors", _convert_to_vertex_list);
GRAPH_PROXY_METHOD(personalized_pagerank, "personalized_pagerank");
GRAPH_PROXY_METHOD(shortest_paths, "shortest_paths");
GRAPH_PROXY_METHOD(strength, "strength");
GRAPH_PROXY_METHOD_PP(successors, "successors", _convert_to_vertex_list);

#undef GRAPH_PROXY_METHOD

#define GRAPH_PROXY_METHOD_SPEC(FUNC, METHODNAME) \
  {METHODNAME, (PyCFunction)igraphmodule_Vertex_##FUNC, METH_VARARGS | METH_KEYWORDS, \
    "Proxy method to L{Graph." METHODNAME "()}\n\n"              \
    "This method calls the " METHODNAME " method of the L{Graph} class " \
    "with this vertex as the first argument, and returns the result.\n\n"\
    "@see: Graph." METHODNAME "() for details."}
#define GRAPH_PROXY_METHOD_SPEC_2(FUNC, METHODNAME, METHODNAME_IN_GRAPH) \
  {METHODNAME, (PyCFunction)igraphmodule_Vertex_##FUNC, METH_VARARGS | METH_KEYWORDS, \
    "Proxy method to L{Graph." METHODNAME_IN_GRAPH "()}\n\n"              \
    "This method calls the " METHODNAME_IN_GRAPH " method of the L{Graph} class " \
    "with this vertex as the first argument, and returns the result.\n\n"\
    "@see: Graph." METHODNAME_IN_GRAPH "() for details."}

/**
 * \ingroup python_interface_vertex
 * Method table for the \c igraph.Vertex object
 */
PyMethodDef igraphmodule_Vertex_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_Vertex_attributes,
    METH_NOARGS,
    "attributes() -> dict\n\n"
    "Returns a dict of attribute names and values for the vertex\n"
  },
  {"attribute_names", (PyCFunction)igraphmodule_Vertex_attribute_names,
    METH_NOARGS,
    "attribute_names() -> list\n\n"
    "Returns the list of vertex attribute names\n"
  },
  {"update_attributes", (PyCFunction)igraphmodule_Vertex_update_attributes,
    METH_VARARGS | METH_KEYWORDS,
    "update_attributes(E, **F) -> None\n\n"
    "Updates the attributes of the vertex from dict/iterable E and F.\n\n"
    "If E has a C{keys()} method, it does: C{for k in E: self[k] = E[k]}.\n"
    "If E lacks a C{keys()} method, it does: C{for (k, v) in E: self[k] = v}.\n"
    "In either case, this is followed by: C{for k in F: self[k] = F[k]}.\n\n"
    "This method thus behaves similarly to the C{update()} method of Python\n"
    "dictionaries."
  },
  GRAPH_PROXY_METHOD_SPEC(betweenness, "betweenness"),
  GRAPH_PROXY_METHOD_SPEC(closeness, "closeness"),
  GRAPH_PROXY_METHOD_SPEC(constraint, "constraint"),
  GRAPH_PROXY_METHOD_SPEC(degree, "degree"),
  GRAPH_PROXY_METHOD_SPEC_2(delete, "delete", "delete_vertices"),
  GRAPH_PROXY_METHOD_SPEC(diversity, "diversity"),
  GRAPH_PROXY_METHOD_SPEC(eccentricity, "eccentricity"),
  GRAPH_PROXY_METHOD_SPEC(get_shortest_paths, "get_shortest_paths"),
  GRAPH_PROXY_METHOD_SPEC(indegree, "indegree"),
  GRAPH_PROXY_METHOD_SPEC(is_minimal_separator, "is_minimal_separator"),
  GRAPH_PROXY_METHOD_SPEC(is_separator, "is_separator"),
  GRAPH_PROXY_METHOD_SPEC(neighbors, "neighbors"),
  GRAPH_PROXY_METHOD_SPEC(outdegree, "outdegree"),
  GRAPH_PROXY_METHOD_SPEC(pagerank, "pagerank"),
  GRAPH_PROXY_METHOD_SPEC(predecessors, "predecessors"),
  GRAPH_PROXY_METHOD_SPEC(personalized_pagerank, "personalized_pagerank"),
  GRAPH_PROXY_METHOD_SPEC(shortest_paths, "shortest_paths"),
  GRAPH_PROXY_METHOD_SPEC(strength, "strength"),
  GRAPH_PROXY_METHOD_SPEC(successors, "successors"),
  {NULL}
};

#undef GRAPH_PROXY_METHOD_SPEC
#undef GRAPH_PROXY_METHOD_SPEC_2

/** \ingroup python_interface_vertex
 * This structure is the collection of functions necessary to implement
 * the vertex as a mapping (i.e. to allow the retrieval and setting of
 * igraph attributes in Python as if it were of a Python mapping type)
 */
PyMappingMethods igraphmodule_Vertex_as_mapping = {
  // returns the number of vertex attributes
  (lenfunc)igraphmodule_Vertex_attribute_count,
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
  {"graph", (getter)igraphmodule_Vertex_get_graph, NULL,
      "The graph the vertex belongs to", NULL
  },
  {NULL}
};

/** \ingroup python_interface_vertex
 * Python type object referencing the methods Python calls when it performs various operations on
 * a vertex of a graph
 */
PyTypeObject igraphmodule_VertexType =
{
  PyVarObject_HEAD_INIT(0, 0)
  "igraph.Vertex",                            /* tp_name */
  sizeof(igraphmodule_VertexObject),          /* tp_basicsize */
  0,                                          /* tp_itemsize */
  (destructor)igraphmodule_Vertex_dealloc,    /* tp_dealloc */
  0,                                          /* tp_print */
  0,                                          /* tp_getattr */
  0,                                          /* tp_setattr */
  0,                                          /* tp_compare (2.x) / tp_reserved (3.x) */
  (reprfunc)igraphmodule_Vertex_repr,         /* tp_repr */
  0,                                          /* tp_as_number */
  0,                                          /* tp_as_sequence */
  &igraphmodule_Vertex_as_mapping,            /* tp_as_mapping */
  (hashfunc)igraphmodule_Vertex_hash,         /* tp_hash */
  0,                                          /* tp_call */
  0,                                          /* tp_str */
  0,                                          /* tp_getattro */
  0,                                          /* tp_setattro */
  0,                                          /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
  "Class representing a single vertex in a graph.\n\n"
  "The vertex is referenced by its index, so if the underlying graph\n"
  "changes, the semantics of the vertex object might change as well\n"
  "(if the vertex indices are altered in the original graph).\n\n"
  "The attributes of the vertex can be accessed by using the vertex\n"
  "as a hash:\n\n"
  "  >>> v[\"color\"] = \"red\"                  #doctest: +SKIP\n"
  "  >>> print v[\"color\"]                      #doctest: +SKIP\n"
  "  red\n", /* tp_doc */
  0,                                          /* tp_traverse */
  0,                                          /* tp_clear */
  (richcmpfunc)igraphmodule_Vertex_richcompare, /* tp_richcompare */
  0,                                          /* tp_weaklistoffset */
  0,                                          /* tp_iter */
  0,                                          /* tp_iternext */
  igraphmodule_Vertex_methods,                /* tp_methods */
  0,                                          /* tp_members */
  igraphmodule_Vertex_getseters,              /* tp_getset */
};

