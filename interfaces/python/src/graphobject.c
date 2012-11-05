/* vim:set ts=4 sw=2 sts=2 et:  */
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
#include "arpackobject.h"
#include "bfsiter.h"
#include "common.h"
#include "convert.h"
#include "edgeseqobject.h"
#include "error.h"
#include "filehandle.h"
#include "graphobject.h"
#include "indexing.h"
#include "memory.h"
#include "py2compat.h"
#include "pyhelpers.h"
#include "vertexseqobject.h"
#include <float.h>

PyTypeObject igraphmodule_GraphType;

#define CREATE_GRAPH(py_graph, c_graph) { \
  py_graph = (igraphmodule_GraphObject *) Py_TYPE(self)->tp_alloc(Py_TYPE(self), 0); \
  if (py_graph != NULL) { \
    igraphmodule_Graph_init_internal(py_graph); \
    py_graph->g = (c_graph); \
  } \
  RC_ALLOC("Graph", py_graph); \
}
#define CREATE_GRAPH_FROM_TYPE(py_graph, c_graph, py_type) { \
  py_graph = (igraphmodule_GraphObject *) py_type->tp_alloc(py_type, 0); \
  if (py_graph != NULL) { \
    igraphmodule_Graph_init_internal(py_graph); \
    py_graph->g = (c_graph); \
  } \
  RC_ALLOC("Graph", py_graph); \
}

/**********************************************************************
 * Basic implementation of igraph.Graph                               *
 **********************************************************************/

/** \defgroup python_interface_graph Graph object
 * \ingroup python_interface */

/**
 * \ingroup python_interface_internal
 * \brief Initializes the internal structures in an \c igraph.Graph object's
 * C representation.
 * 
 * This function must be called whenever we create a new Graph object with
 * \c tp_alloc
 */
void igraphmodule_Graph_init_internal(igraphmodule_GraphObject * self)
{
  if (!self) return;
  self->destructor = NULL;
  self->weakreflist = NULL;
}

/**
 * \ingroup python_interface_graph
 * \brief Creates a new igraph object in Python
 * 
 * This function is called whenever a new \c igraph.Graph object is created in
 * Python. An optional \c n parameter can be passed from Python,
 * representing the number of vertices in the graph. If it is omitted,
 * the default value is 0.
 * 
 * <b>Example call from Python:</b>
\verbatim
g = igraph.Graph(5);
\endverbatim
 *
 * In fact, the parameters are processed by \c igraphmodule_Graph_init
 * 
 * \return the new \c igraph.Graph object or NULL if an error occurred.
 * 
 * \sa igraphmodule_Graph_init
 * \sa igraph_empty
 */
PyObject *igraphmodule_Graph_new(PyTypeObject * type, PyObject * args,
                                 PyObject * kwds)
{
  igraphmodule_GraphObject *self;

  self = (igraphmodule_GraphObject *) type->tp_alloc(type, 0);
  RC_ALLOC("Graph", self);

  /* don't need it, the constructor will do it */
  /*if (self != NULL) {
     igraph_empty(&self->g, 1, 0);
     } */
  igraphmodule_Graph_init_internal(self);

  return (PyObject *) self;
}

/**
 * \ingroup python_interface_graph
 * \brief Clears the graph object's subobject (before deallocation)
 */
int igraphmodule_Graph_clear(igraphmodule_GraphObject * self)
{
  PyObject *tmp;
  PyObject_GC_UnTrack(self);

  tmp = self->destructor;
  self->destructor = NULL;
  Py_XDECREF(tmp);

  return 0;
}

/**
 * \ingroup python_interface_graph
 * \brief Support for cyclic garbage collection in Python
 * 
 * This is necessary because the \c igraph.Graph object contains several
 * other \c PyObject pointers and they might point back to itself.
 */
int igraphmodule_Graph_traverse(igraphmodule_GraphObject * self,
                                visitproc visit, void *arg)
{
  int vret, i;

  RC_TRAVERSE("Graph", self);

  if (self->destructor) {
    vret = visit(self->destructor, arg);
    if (vret != 0)
      return vret;
  }

  if (self->g.attr) {
    for (i = 0; i < 3; i++) {
      vret = visit(((PyObject **) (self->g.attr))[i], arg);
      if (vret != 0)
        return vret;
    }
  }

  return 0;
}

/**
 * \ingroup python_interface_graph
 * \brief Deallocates a Python representation of a given igraph object
 */
void igraphmodule_Graph_dealloc(igraphmodule_GraphObject * self)
{
  PyObject *r;

  /* Clear weak references */
  if (self->weakreflist != NULL)
     PyObject_ClearWeakRefs((PyObject *) self);

  igraph_destroy(&self->g);

  if (PyCallable_Check(self->destructor)) {
    r = PyObject_CallObject(self->destructor, NULL);
    if (r) {
      Py_DECREF(r);
    }
  }

  igraphmodule_Graph_clear(self);

  RC_DEALLOC("Graph", self);

  Py_TYPE(self)->tp_free((PyObject*)self);
}

/**
 * \ingroup python_interface_graph
 * \brief Initializes a new \c igraph object in Python
 * 
 * This function is called whenever a new \c igraph.Graph object is initialized in
 * Python (note that initializing is not equal to creating: an object might
 * be created but not initialized when it is being recovered from a serialized
 * state).
 * 
 * Throws \c AssertionError in Python if \c vcount is less than or equal to zero.
 * \return the new \c igraph.Graph object or NULL if an error occurred.
 * 
 * \sa igraphmodule_Graph_new
 * \sa igraph_empty
 * \sa igraph_create
 */
int igraphmodule_Graph_init(igraphmodule_GraphObject * self,
                            PyObject * args, PyObject * kwds) {
  static char *kwlist[] = { "n", "edges", "directed", NULL };
  long int n = 0;
  PyObject *edges = NULL, *dir = Py_False;
  igraph_vector_t edges_vector;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lO!O", kwlist,
                                   &n, &PyList_Type, &edges, &dir))
    return -1;

  if (edges && PyList_Check(edges)) {
    /* Caller specified an edge list, so we use igraph_create */
    /* We have to convert the Python list to a igraph_vector_t */
    if (igraphmodule_PyObject_to_edgelist(edges, &edges_vector, 0)) {
      igraphmodule_handle_igraph_error();
      return -1;
    }

    if (igraph_create
        (&self->g, &edges_vector, (igraph_integer_t) n, PyObject_IsTrue(dir))) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(&edges_vector);
      return -1;
    }

    igraph_vector_destroy(&edges_vector);
  }
  else {
    /* No edge list was specified, let's use igraph_empty */
    if (igraph_empty(&self->g, (igraph_integer_t) n, PyObject_IsTrue(dir))) {
      igraphmodule_handle_igraph_error();
      return -1;
    }
  }

  return 0;
}

/** \ingroup python_interface_graph
 * \brief Creates an \c igraph.Graph object from an existing \c igraph_t
 *
 * The newly created \c igraph.Graph object will take ownership of the
 * given \c igraph_t. This function is not accessible from Python the
 * normal way, but it is exposed via the C API of the Python module.
 * See \c api.h for more details.
 */
PyObject* igraphmodule_Graph_from_igraph_t(igraph_t *graph) {
  igraphmodule_GraphObject* result;
  PyTypeObject* type = &igraphmodule_GraphType;

  CREATE_GRAPH_FROM_TYPE(result, *graph, type);

  return (PyObject*)result;
}

/** \ingroup python_interface_graph
 * \brief Formats an \c igraph.Graph object in a human-readable format.
 * 
 * This function is rather simple now, it returns the number of vertices
 * and edges in a string.
 * 
 * \return the formatted textual representation as a \c PyObject
 */
PyObject *igraphmodule_Graph_str(igraphmodule_GraphObject * self)
{
  if (igraph_is_directed(&self->g))
    return PyString_FromFormat("Directed graph (|V| = %ld, |E| = %ld)",
                               (long)igraph_vcount(&self->g),
                               (long)igraph_ecount(&self->g));
  else
    return PyString_FromFormat("Undirected graph (|V| = %ld, |E| = %ld)",
                               (long)igraph_vcount(&self->g),
                               (long)igraph_ecount(&self->g));
}

/** \ingroup python_interface_copy
 * \brief Creates an exact deep copy of the graph
 * \return the copy of the graph
 */
PyObject *igraphmodule_Graph_copy(igraphmodule_GraphObject * self)
{
  igraphmodule_GraphObject *result;
  igraph_t g;

  if (igraph_copy(&g, &self->g)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/**********************************************************************
 * The most basic igraph interface                                    *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Returns the number of vertices in an \c igraph.Graph object.
 * \return the number of vertices as a \c PyObject
 * \sa igraph_vcount
 */
PyObject *igraphmodule_Graph_vcount(igraphmodule_GraphObject * self)
{
  PyObject *result;
  result = Py_BuildValue("l", (long)igraph_vcount(&self->g));
  return result;
}

/** \ingroup python_interface_graph
 * \brief Returns the number of edges in an \c igraph.Graph object.
 * \return the number of edges as a \c PyObject
 * \sa igraph_ecount
 */
PyObject *igraphmodule_Graph_ecount(igraphmodule_GraphObject * self)
{
  PyObject *result;
  result = Py_BuildValue("l", (long)igraph_ecount(&self->g));
  return result;
}

/** \ingroup python_interface_graph
 * \brief Checks whether an \c igraph.Graph object is a DAG.
 * \return \c True if the graph is directed, \c False otherwise.
 * \sa igraph_is_dag
 */
PyObject *igraphmodule_Graph_is_dag(igraphmodule_GraphObject * self)
{
  igraph_bool_t res;

  if (igraph_is_dag(&self->g, &res)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (res)
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Checks whether an \c igraph.Graph object is directed.
 * \return \c True if the graph is directed, \c False otherwise.
 * \sa igraph_is_directed
 */
PyObject *igraphmodule_Graph_is_directed(igraphmodule_GraphObject * self)
{
  if (igraph_is_directed(&self->g))
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/**
 * \ingroup python_interface_graph
 * \brief Checks whether a matching is valid in the context of an \c igraph.Graph
 * object.
 * \sa igraph_is_matching
 */
PyObject *igraphmodule_Graph_is_matching(igraphmodule_GraphObject* self,
    PyObject* args, PyObject* kwds) {
  static char* kwlist[] = { "matching", "types", NULL };
  PyObject *matching_o, *types_o = Py_None;
  igraph_vector_long_t* matching = 0;
  igraph_vector_bool_t* types = 0;
  igraph_bool_t result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &matching_o,
        &types_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_long_t(matching_o, self, &matching,
        ATTRIBUTE_TYPE_VERTEX))
	return NULL;

  if (igraphmodule_attrib_to_vector_bool_t(types_o, self, &types, ATTRIBUTE_TYPE_VERTEX)) {
    if (matching != 0) { igraph_vector_long_destroy(matching); free(matching); }
	return NULL;
  }

  if (igraph_is_matching(&self->g, types, matching, &result)) {
    if (matching != 0) { igraph_vector_long_destroy(matching); free(matching); }
    if (types != 0) { igraph_vector_bool_destroy(types); free(types); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (matching != 0) { igraph_vector_long_destroy(matching); free(matching); }
  if (types != 0) { igraph_vector_bool_destroy(types); free(types); }

  if (result)
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/**
 * \ingroup python_interface_graph
 * \brief Checks whether a matching is valid and maximal in the context of an
 *        \c igraph.Graph object.
 * \sa igraph_is_maximal_matching
 */
PyObject *igraphmodule_Graph_is_maximal_matching(igraphmodule_GraphObject* self,
    PyObject* args, PyObject* kwds) {
  static char* kwlist[] = { "matching", "types", NULL };
  PyObject *matching_o, *types_o = Py_None;
  igraph_vector_long_t* matching = 0;
  igraph_vector_bool_t* types = 0;
  igraph_bool_t result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &matching_o,
        &types_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_long_t(matching_o, self, &matching,
        ATTRIBUTE_TYPE_VERTEX))
	return NULL;

  if (igraphmodule_attrib_to_vector_bool_t(types_o, self, &types, ATTRIBUTE_TYPE_VERTEX)) {
    if (matching != 0) { igraph_vector_long_destroy(matching); free(matching); }
	return NULL;
  }

  if (igraph_is_maximal_matching(&self->g, types, matching, &result)) {
    if (matching != 0) { igraph_vector_long_destroy(matching); free(matching); }
    if (types != 0) { igraph_vector_bool_destroy(types); free(types); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (matching != 0) { igraph_vector_long_destroy(matching); free(matching); }
  if (types != 0) { igraph_vector_bool_destroy(types); free(types); }

  if (result)
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}


/** \ingroup python_interface_graph
 * \brief Checks whether an \c igraph.Graph object is simple.
 * \return \c True if the graph is simple, \c False otherwise.
 * \sa igraph_is_simple
 */
PyObject *igraphmodule_Graph_is_simple(igraphmodule_GraphObject *self) {
  igraph_bool_t res;

  if (igraph_is_simple(&self->g, &res)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (res)
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Adds vertices to an \c igraph.Graph
 * \return the extended \c igraph.Graph object
 * \sa igraph_add_vertices
 */
PyObject *igraphmodule_Graph_add_vertices(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds) {
  long n;

  if (!PyArg_ParseTuple(args, "l", &n))
    return NULL;

  if (igraph_add_vertices(&self->g, (igraph_integer_t) n, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Removes vertices from an \c igraph.Graph
 * \return the modified \c igraph.Graph object
 * 
 * \todo Need more error checking on vertex IDs. (igraph fails when an
 * invalid vertex ID is given)
 * \sa igraph_delete_vertices
 */
PyObject *igraphmodule_Graph_delete_vertices(igraphmodule_GraphObject * self,
                                             PyObject * args, PyObject * kwds) {
  PyObject *list;
  igraph_vs_t vs;

  if (!PyArg_ParseTuple(args, "O", &list)) return NULL;
  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, 0, 0)) return NULL;

  if (igraph_delete_vertices(&self->g, vs)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    return NULL;
  }

  igraph_vs_destroy(&vs);
  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Adds edges to an \c igraph.Graph
 * \return the extended \c igraph.Graph object
 * 
 * \todo Need more error checking on vertex IDs. (igraph fails when an
 * invalid vertex ID is given)
 * \sa igraph_add_edges
 */
PyObject *igraphmodule_Graph_add_edges(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  PyObject *list;
  igraph_vector_t v;

  if (!PyArg_ParseTuple(args, "O", &list))
    return NULL;

  if (igraphmodule_PyObject_to_edgelist(list, &v, &self->g))
    return NULL;

  /* do the hard work :) */
  if (igraph_add_edges(&self->g, &v, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&v);
    return NULL;
  }

  igraph_vector_destroy(&v);
  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Deletes edges from an \c igraph.Graph
 * \return the extended \c igraph.Graph object
 * 
 * \todo Need more error checking on vertex IDs. (igraph fails when an
 * invalid vertex ID is given)
 * \sa igraph_delete_edges
 */
PyObject *igraphmodule_Graph_delete_edges(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  PyObject *list;
  igraph_es_t es;
  static char *kwlist[] = { "edges", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &list))
    return NULL;
	
  if (igraphmodule_PyObject_to_es_t(list, &es, &self->g, 0)) {
    /* something bad happened during conversion, return immediately */
    return NULL;
  }

  if (igraph_delete_edges(&self->g, es)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    return NULL;
  }

  igraph_es_destroy(&es);
  Py_RETURN_NONE;
}

/**********************************************************************
 * tructural properties                                              *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief The degree of some vertices in an \c igraph.Graph
 * \return the degree list as a Python object
 * \sa igraph_degree
 */
PyObject *igraphmodule_Graph_degree(igraphmodule_GraphObject * self,
                                    PyObject * args, PyObject * kwds)
{
  PyObject *list = Py_None;
  PyObject *loops = Py_True;
  PyObject *dtype_o = Py_None;
  PyObject *dmode_o = Py_None;
  igraph_neimode_t dmode = IGRAPH_ALL;
  igraph_vector_t result;
  igraph_vs_t vs;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "vertices", "mode", "loops", "type", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO", kwlist,
                                   &list, &dmode_o, &loops, &dtype_o))
    return NULL;

  if (dmode_o == Py_None && dtype_o != Py_None) {
    dmode_o = dtype_o;
    PY_IGRAPH_DEPRECATED("type=... keyword argument is deprecated since igraph 0.6, use mode=... instead");
  }

  if (igraphmodule_PyObject_to_neimode_t(dmode_o, &dmode))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, &return_single, 0)) {
    return NULL;
  }

  if (igraph_vector_init(&result, 0)) {
    igraph_vs_destroy(&vs);
    return NULL;
  }

  if (igraph_degree(&self->g, &result, vs,
                    dmode, PyObject_IsTrue(loops))) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  else
    list = PyInt_FromLong((long int)VECTOR(result)[0]);

  igraph_vector_destroy(&result);
  igraph_vs_destroy(&vs);

  return list;
}

/**
 * \ingroup python_interface_graph
 * \brief Structural diversity index of some vertices in an \c igraph.Graph
 * \sa igraph_diversity
 */
PyObject *igraphmodule_Graph_diversity(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds) {
  PyObject *list = Py_None;
  PyObject *weights_o = Py_None;
  igraph_vector_t result, *weights = 0;
  igraph_vs_t vs;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "vertices", "weights", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist,
                                   &list, &weights_o))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&result, 0)) {
    igraph_vs_destroy(&vs);
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) {
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (igraph_diversity(&self->g, weights, &result, vs)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&result);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_FLOAT);
  else
    list = PyFloat_FromDouble(VECTOR(result)[0]);

  igraph_vector_destroy(&result);
  igraph_vs_destroy(&vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief The strength (weighted degree) of some vertices in an \c igraph.Graph
 * \return the strength list as a Python object
 * \sa igraph_strength
 */
PyObject *igraphmodule_Graph_strength(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  PyObject *list = Py_None;
  PyObject *loops = Py_True;
  PyObject *dtype_o = Py_None;
  PyObject *dmode_o = Py_None;
  PyObject *weights_o = Py_None;
  igraph_neimode_t dmode = IGRAPH_ALL;
  igraph_vector_t result, *weights = 0;
  igraph_vs_t vs;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "vertices", "mode", "loops", "weights",
    "type", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOOO", kwlist,
                                   &list, &dmode_o, &loops, &weights_o,
                                   &dtype_o))
    return NULL;

  if (dmode_o == Py_None && dtype_o != Py_None) {
    dmode_o = dtype_o;
    PY_IGRAPH_DEPRECATED("type=... keyword argument is deprecated since igraph 0.6, use mode=... instead");
  }

  if (igraphmodule_PyObject_to_neimode_t(dmode_o, &dmode))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&result, 0)) {
    igraph_vs_destroy(&vs);
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) {
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (igraph_strength(&self->g, &result, vs, dmode,
        PyObject_IsTrue(loops), weights)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&result);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_FLOAT);
  else
    list = PyFloat_FromDouble(VECTOR(result)[0]);

  igraph_vector_destroy(&result);
  igraph_vs_destroy(&vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the graph density
 * \return the density
 * \sa igraph_density
 */
PyObject *igraphmodule_Graph_density(igraphmodule_GraphObject * self,
                                     PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "loops", NULL };
  igraph_real_t result;
  PyObject *loops = Py_False;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &loops))
    return NULL;

  if (igraph_density(&self->g, &result, PyObject_IsTrue(loops))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  return Py_BuildValue("d", (double)result);
}

/** \ingroup python_interface_graph
 * \brief The maximum degree of some vertices in an \c igraph.Graph
 * \return the maxium degree as a Python object
 * \sa igraph_maxdegree
 */
PyObject *igraphmodule_Graph_maxdegree(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  PyObject *list = Py_None;
  igraph_neimode_t dmode = IGRAPH_ALL;
  PyObject *dtype_o = Py_None;
  PyObject *dmode_o = Py_None;
  PyObject *loops = Py_False;
  igraph_integer_t result;
  igraph_vs_t vs;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "vertices", "mode", "loops", "type", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO", kwlist,
                                   &list, &dmode_o, &loops, &dtype_o))
    return NULL;

  if (dmode_o == Py_None && dtype_o != Py_None) {
    dmode_o = dtype_o;
    PY_IGRAPH_DEPRECATED("type=... keyword argument is deprecated since igraph 0.6, use mode=... instead");
  }

  if (igraphmodule_PyObject_to_neimode_t(dmode_o, &dmode)) return NULL;

  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_maxdegree(&self->g, &result, vs,
                       dmode, PyObject_IsTrue(loops))) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    return NULL;
  }

  igraph_vs_destroy(&vs);

  return PyInt_FromLong((long)result);
}

/** \ingroup python_interface_graph
 * \brief Checks whether an edge is a loop edge 
 * \return a boolean or a list of booleans
 * \sa igraph_is_loop
 */
PyObject *igraphmodule_Graph_is_loop(igraphmodule_GraphObject *self,
                                     PyObject *args, PyObject *kwds) {
  PyObject *list = Py_None;
  igraph_vector_bool_t result;
  igraph_es_t es;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "edges", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &list))
    return NULL;

  if (igraphmodule_PyObject_to_es_t(list, &es, &self->g, &return_single)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_bool_init(&result, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    return NULL;
  }

  if (igraph_is_loop(&self->g, &result, es)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    igraph_vector_bool_destroy(&result);
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_bool_t_to_PyList(&result);
  else {
    list = (VECTOR(result)[0]) ? Py_True : Py_False;
    Py_INCREF(list);
  }

  igraph_vector_bool_destroy(&result);
  igraph_es_destroy(&es);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Checks whether an edge is a multiple edge 
 * \return a boolean or a list of booleans
 * \sa igraph_is_multiple
 */
PyObject *igraphmodule_Graph_is_multiple(igraphmodule_GraphObject *self,
                                         PyObject *args, PyObject *kwds) {
  PyObject *list = Py_None;
  igraph_vector_bool_t result;
  igraph_es_t es;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "edges", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &list))
    return NULL;

  if (igraphmodule_PyObject_to_es_t(list, &es, &self->g, &return_single)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_bool_init(&result, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    return NULL;
  }

  if (igraph_is_multiple(&self->g, &result, es)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    igraph_vector_bool_destroy(&result);
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_bool_t_to_PyList(&result);
  else {
    list = (VECTOR(result)[0]) ? Py_True : Py_False;
    Py_INCREF(list);
  }

  igraph_vector_bool_destroy(&result);
  igraph_es_destroy(&es);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Checks whether an edge is mutual 
 * \return a boolean or a list of booleans
 * \sa igraph_is_mutual
 */
PyObject *igraphmodule_Graph_is_mutual(igraphmodule_GraphObject *self,
                                       PyObject *args, PyObject *kwds) {
  PyObject *list = Py_None;
  igraph_vector_bool_t result;
  igraph_es_t es;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "edges", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &list))
    return NULL;

  if (igraphmodule_PyObject_to_es_t(list, &es, &self->g, &return_single)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_bool_init(&result, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    return NULL;
  }

  if (igraph_is_mutual(&self->g, &result, es)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    igraph_vector_bool_destroy(&result);
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_bool_t_to_PyList(&result);
  else {
    list = (VECTOR(result)[0]) ? Py_True : Py_False;
    Py_INCREF(list);
  }

  igraph_vector_bool_destroy(&result);
  igraph_es_destroy(&es);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Checks whether an \c igraph.Graph object has multiple edges.
 * \return \c True if the graph has multiple edges, \c False otherwise.
 * \sa igraph_has_multiple
 */
PyObject *igraphmodule_Graph_has_multiple(igraphmodule_GraphObject *self) {
  igraph_bool_t res;

  if (igraph_has_multiple(&self->g, &res)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (res)
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Checks the multiplicity of the edges 
 * \return the edge multiplicities as a Python list
 * \sa igraph_count_multiple
 */
PyObject *igraphmodule_Graph_count_multiple(igraphmodule_GraphObject *self,
                                            PyObject *args, PyObject *kwds) {
  PyObject *list = Py_None;
  igraph_vector_t result;
  igraph_es_t es;
  igraph_bool_t return_single = 0;

  static char *kwlist[] = { "edges", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &list))
    return NULL;

  if (igraphmodule_PyObject_to_es_t(list, &es, &self->g, &return_single)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&result, 0)) {
    igraph_es_destroy(&es);
    return NULL;
  }

  if (igraph_count_multiple(&self->g, &result, es)) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  else
    list = PyInt_FromLong((long int)VECTOR(result)[0]);

  igraph_vector_destroy(&result);
  igraph_es_destroy(&es);

  return list;
}

/** \ingroup python_interface_graph
 * \brief The neighbors of a given vertex in an \c igraph.Graph
 * This method accepts a single vertex ID as a parameter, and returns the
 * neighbors of the given vertex in the form of an integer list. A
 * second argument may be passed as well, meaning the type of neighbors to
 * be returned (\c OUT for successors, \c IN for predecessors or \c ALL
 * for both of them). This argument is ignored for undirected graphs.
 * 
 * \return the neighbor list as a Python list object
 * \sa igraph_neighbors
 */
PyObject *igraphmodule_Graph_neighbors(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  PyObject *list, *dtype_o=Py_None, *dmode_o=Py_None, *index_o;
  igraph_neimode_t dmode = IGRAPH_ALL;
  igraph_integer_t idx;
  igraph_vector_t result;

  static char *kwlist[] = { "vertex", "mode", "type", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OO", kwlist,
        &index_o, &dmode_o, &dtype_o))
    return NULL;

  if (dmode_o == Py_None && dtype_o != Py_None) {
    dmode_o = dtype_o;
    PY_IGRAPH_DEPRECATED("type=... keyword argument is deprecated since igraph 0.6, use mode=... instead");
  }

  if (igraphmodule_PyObject_to_neimode_t(dmode_o, &dmode))
    return NULL;

  if (igraphmodule_PyObject_to_vid(index_o, &idx, &self->g))
    return NULL;

  if (igraph_vector_init(&result, 1))
    return igraphmodule_handle_igraph_error();

  if (igraph_neighbors(&self->g, &result, idx, dmode)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief The incident edges of a given vertex in an \c igraph.Graph
 * This method accepts a single vertex ID as a parameter, and returns the
 * IDs of the incident edges of the given vertex in the form of an integer list.
 * A second argument may be passed as well, meaning the type of neighbors to
 * be returned (\c OUT for successors, \c IN for predecessors or \c ALL
 * for both of them). This argument is ignored for undirected graphs.
 * 
 * \return the adjacency list as a Python list object
 * \sa igraph_incident
 */
PyObject *igraphmodule_Graph_incident(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  PyObject *list, *dmode_o = Py_None, *dtype_o = Py_None, *index_o;
  igraph_neimode_t dmode = IGRAPH_OUT;
  igraph_integer_t idx;
  igraph_vector_t result;

  static char *kwlist[] = { "vertex", "mode", "type", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OO", kwlist,
        &index_o, &dmode_o, &dtype_o))
    return NULL;

  if (dmode_o == Py_None && dtype_o != Py_None) {
    dmode_o = dtype_o;
    PY_IGRAPH_DEPRECATED("type=... keyword argument is deprecated since igraph 0.6, use mode=... instead");
  }

  if (igraphmodule_PyObject_to_neimode_t(dmode_o, &dmode))
    return NULL;

  if (igraphmodule_PyObject_to_vid(index_o, &idx, &self->g))
    return NULL;

  igraph_vector_init(&result, 1);
  if (igraph_incident(&self->g, &result, idx, dmode)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the graph reciprocity
 * \return the reciprocity
 * \sa igraph_reciprocity
 */
PyObject *igraphmodule_Graph_reciprocity(igraphmodule_GraphObject * self,
                                         PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "ignore_loops", "mode", NULL };
  igraph_real_t result;
  igraph_reciprocity_t mode = IGRAPH_RECIPROCITY_DEFAULT;
  PyObject *ignore_loops = Py_True, *mode_o = Py_None;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &ignore_loops, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_reciprocity_t(mode_o, &mode))
    return NULL;

  if (igraph_reciprocity(&self->g, &result, PyObject_IsTrue(ignore_loops), mode)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  return Py_BuildValue("d", (double)result);
}

/** \ingroup python_interface_graph
 * \brief The successors of a given vertex in an \c igraph.Graph
 * This method accepts a single vertex ID as a parameter, and returns the
 * successors of the given vertex in the form of an integer list. It
 * is equivalent to calling \c igraph.Graph.neighbors with \c type=OUT
 * 
 * \return the successor list as a Python list object
 * \sa igraph_neighbors
 */
PyObject *igraphmodule_Graph_successors(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  PyObject *list, *index_o;
  igraph_integer_t idx;
  igraph_vector_t result;

  static char *kwlist[] = { "vertex", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &index_o))
    return NULL;

  if (igraphmodule_PyObject_to_vid(index_o, &idx, &self->g))
    return NULL;

  igraph_vector_init(&result, 1);
  if (igraph_neighbors(&self->g, &result, idx, IGRAPH_OUT)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief The predecessors of a given vertex in an \c igraph.Graph
 * This method accepts a single vertex ID as a parameter, and returns the
 * predecessors of the given vertex in the form of an integer list. It
 * is equivalent to calling \c igraph.Graph.neighbors with \c type=IN
 * 
 * \return the predecessor list as a Python list object
 * \sa igraph_neighbors
 */
PyObject *igraphmodule_Graph_predecessors(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  PyObject *list, *index_o;
  igraph_integer_t idx;
  igraph_vector_t result;

  static char *kwlist[] = { "vertex", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &index_o))
    return NULL;

  if (igraphmodule_PyObject_to_vid(index_o, &idx, &self->g))
    return NULL;

  igraph_vector_init(&result, 1);
  if (igraph_neighbors(&self->g, &result, idx, IGRAPH_IN)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Decides whether a graph is connected.
 * \return Py_True if the graph is connected, Py_False otherwise
 * \sa igraph_is_connected
 */
PyObject *igraphmodule_Graph_is_connected(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "mode", NULL };
  igraph_connectedness_t mode = IGRAPH_STRONG;
  igraph_bool_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|l", kwlist, &mode))
    return NULL;

  if (mode != IGRAPH_STRONG && mode != IGRAPH_WEAK) {
    PyErr_SetString(PyExc_ValueError, "mode must be either STRONG or WEAK");
    return NULL;
  }

  if (igraph_is_connected(&self->g, &res, mode)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (res)
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Decides whether there is an edge from a given vertex to an other one.
 * \return Py_True if the vertices are directly connected, Py_False otherwise
 * \sa igraph_are_connected
 */
PyObject *igraphmodule_Graph_are_connected(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "v1", "v2", NULL };
  PyObject *v1, *v2;
  igraph_integer_t idx1, idx2;
  igraph_bool_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist, &v1, &v2))
    return NULL;

  if (igraphmodule_PyObject_to_vid(v1, &idx1, &self->g))
    return NULL;

  if (igraphmodule_PyObject_to_vid(v2, &idx2, &self->g))
    return NULL;

  if (igraph_are_connected(&self->g, idx1, idx2, &res))
    return NULL;

  if (res)
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Returns the ID of an arbitrary edge between the given two vertices
 * \sa igraph_get_eid
 */
PyObject *igraphmodule_Graph_get_eid(igraphmodule_GraphObject * self,
                                     PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "v1", "v2", "directed", "error", NULL };
  long int v1, v2;
  igraph_integer_t result;
  PyObject *directed = Py_True;
  PyObject *error = Py_True;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|OO", kwlist, &v1, &v2,
                                   &directed, &error))
    return NULL;
  if (igraph_get_eid(&self->g, &result, (igraph_integer_t) v1, (igraph_integer_t) v2,
        PyObject_IsTrue(directed), PyObject_IsTrue(error)))
    return igraphmodule_handle_igraph_error();

  return Py_BuildValue("l", (long)result);
}

/** \ingroup python_interface_graph
 * \brief Returns the IDs of some edges between some vertices
 * \sa igraph_get_eids
 */
PyObject *igraphmodule_Graph_get_eids(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "pairs", "path", "directed", "error", NULL };
  PyObject *pairs_o = Py_None, *path_o = Py_None;
  PyObject *directed = Py_True;
  PyObject *error = Py_True;
  PyObject *result = NULL;
  igraph_vector_t pairs, path, res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO", kwlist,
                                   &pairs_o, &path_o, &directed,
                                   &error))
    return NULL;

  if (igraph_vector_init(&res, 0))
    return igraphmodule_handle_igraph_error();

  if (pairs_o != Py_None) {
    if (igraphmodule_PyObject_to_edgelist(pairs_o, &pairs, &self->g)) {
      igraph_vector_destroy(&res);
      return NULL;
    }
  }

  if (path_o != Py_None) {
    if (igraphmodule_PyObject_to_vector_t(path_o, &path, 1, 0)) {
      igraph_vector_destroy(&res);
      if (pairs_o != Py_None) igraph_vector_destroy(&pairs);
      return NULL;
    }
  }

  if (igraph_get_eids(&self->g, &res,
        pairs_o == Py_None ? 0 : &pairs,
        path_o  == Py_None ? 0 : &path,
        PyObject_IsTrue(directed),
        PyObject_IsTrue(error))) {
    if (pairs_o != Py_None) igraph_vector_destroy(&pairs);
    if (path_o != Py_None)  igraph_vector_destroy(&path);
    igraph_vector_destroy(&res);
    return igraphmodule_handle_igraph_error();
  }

  if (pairs_o != Py_None) igraph_vector_destroy(&pairs);
  if (path_o != Py_None)  igraph_vector_destroy(&path);

  result = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&res);
  return result;
}

/** \ingroup python_interface_graph
 * \brief Calculates the diameter of an \c igraph.Graph
 * This method accepts two optional parameters: the first one is
 * a boolean meaning whether to consider directed paths (and is
 * ignored for undirected graphs). The second one is only meaningful
 * in unconnected graphs: it is \c True if the longest geodesic
 * within a component should be returned and \c False if the number of
 * vertices should be returned. They both have a default value of \c False.
 * 
 * \return the diameter as a Python integer
 * \sa igraph_diameter
 */
PyObject *igraphmodule_Graph_diameter(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  PyObject *dir = Py_True, *vcount_if_unconnected = Py_True;
  PyObject *weights_o = Py_None;
  igraph_vector_t *weights = 0;

  static char *kwlist[] = {
    "directed", "unconn", "weights", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist,
                                   &dir, &vcount_if_unconnected,
                                   &weights_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  if (weights) {
    igraph_real_t i;
    if (igraph_diameter_dijkstra(&self->g, weights, &i, 0, 0, 0,
          PyObject_IsTrue(dir), PyObject_IsTrue(vcount_if_unconnected))) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(weights); free(weights);
      return NULL;
    }
    igraph_vector_destroy(weights); free(weights);
    return PyFloat_FromDouble((double)i);
  } else {
    igraph_integer_t i;
    if (igraph_diameter(&self->g, &i, 0, 0, 0, PyObject_IsTrue(dir),
                        PyObject_IsTrue(vcount_if_unconnected))) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
    return PyInt_FromLong((long)i);
  }
}

/** \ingroup python_interface_graph
 * \brief Returns a path of the actual diameter of the graph
 * \sa igraph_diameter
 */
PyObject *igraphmodule_Graph_get_diameter(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  PyObject *dir = Py_True, *vcount_if_unconnected = Py_True, *result;
  PyObject *weights_o = Py_None;
  igraph_vector_t *weights = 0;
  igraph_vector_t res;

  static char *kwlist[] = { "directed", "unconn", "weights", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist,
                                   &dir, &vcount_if_unconnected,
                                   &weights_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  igraph_vector_init(&res, 0);
  if (weights) {
    if (igraph_diameter_dijkstra(&self->g, weights, 0, 0, 0, &res,
          PyObject_IsTrue(dir), PyObject_IsTrue(vcount_if_unconnected))) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(weights); free(weights);
      igraph_vector_destroy(&res);
      return NULL;
    }
    igraph_vector_destroy(weights); free(weights);
  } else {
    if (igraph_diameter(&self->g, 0, 0, 0, &res, PyObject_IsTrue(dir),
                        PyObject_IsTrue(vcount_if_unconnected))) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  result = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&res);
  return result;
}

/** \ingroup python_interface_graph
 * \brief Returns the farthest points and their distances in the graph
 * \sa igraph_distance
 */
PyObject *igraphmodule_Graph_farthest_points(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  PyObject *dir = Py_True, *vcount_if_unconnected = Py_True;
  PyObject *weights_o = Py_None;
  igraph_vector_t *weights = 0;
  igraph_integer_t from, to, len;
  igraph_real_t len_real;

  static char *kwlist[] = { "directed", "unconn", "weights", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist,
                                   &dir, &vcount_if_unconnected,
                                   &weights_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  if (weights) {
    if (igraph_diameter_dijkstra(&self->g, weights, &len_real, &from, &to, 0,
          PyObject_IsTrue(dir), PyObject_IsTrue(vcount_if_unconnected))) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(weights); free(weights);
      return NULL;
    }
    igraph_vector_destroy(weights); free(weights);
    if (from >= 0)
      return Py_BuildValue("lld", (long)from, (long)to, (double)len_real);
    return Py_BuildValue("OOd", Py_None, Py_None, (double)len_real);
  } else {
    if (igraph_diameter(&self->g, &len, &from, &to, 0, PyObject_IsTrue(dir),
                        PyObject_IsTrue(vcount_if_unconnected))) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (from >= 0)
      return Py_BuildValue("lll", (long)from, (long)to, (long)len);
    return Py_BuildValue("OOl", Py_None, Py_None, (long)len);
  }
}

/**
 * \ingroup python_interface_graph
 * \brief Calculates the girth of an \c igraph.Graph
 */
PyObject *igraphmodule_Graph_girth(igraphmodule_GraphObject *self,
                                   PyObject *args, PyObject *kwds)
{
  PyObject *sc = Py_False;
  static char *kwlist[] = { "return_shortest_circle", NULL };
  igraph_integer_t girth;
  igraph_vector_t vids;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &sc))
  return NULL;

  igraph_vector_init(&vids, 0);
  if (igraph_girth(&self->g, &girth, &vids)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&vids);
    return NULL;
  }

  if (PyObject_IsTrue(sc)) {
    PyObject* o;
    o=igraphmodule_vector_t_to_PyList(&vids, IGRAPHMODULE_TYPE_INT);
    igraph_vector_destroy(&vids);
    return o;
  }
  return PyInt_FromLong((long)girth);
}

/**
 * \ingroup python_interface_graph
 * \brief Calculates the convergence degree of the edges in a graph 
 */
PyObject *igraphmodule_Graph_convergence_degree(igraphmodule_GraphObject *self) {
  igraph_vector_t result;
	PyObject *o;

  igraph_vector_init(&result, 0);
  if (igraph_convergence_degree(&self->g, &result, 0, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  o=igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&result);
  return o;
}

/**
 * \ingroup python_interface_graph
 * \brief Calculates the sizes of the convergence fields in a graph 
 */
PyObject *igraphmodule_Graph_convergence_field_size(igraphmodule_GraphObject *self) {
  igraph_vector_t ins, outs;
  PyObject *o1, *o2;

  igraph_vector_init(&ins, 0);
  igraph_vector_init(&outs, 0);
  if (igraph_convergence_degree(&self->g, 0, &ins, &outs)) {
	igraphmodule_handle_igraph_error();
	igraph_vector_destroy(&ins);
	igraph_vector_destroy(&outs);
	return NULL;
  }

  o1=igraphmodule_vector_t_to_PyList(&ins, IGRAPHMODULE_TYPE_INT);
  o2=igraphmodule_vector_t_to_PyList(&outs, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&ins);
  igraph_vector_destroy(&outs);
  return Py_BuildValue("NN", o1, o2);
}

/**
 * \ingroup python_interface_graph
 * \brief Calculates the average nearest neighbor degree of the vertices
 *   of a \c igraph.Graph
 */
PyObject *igraphmodule_Graph_knn(igraphmodule_GraphObject *self,
                                 PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "vids", "weights", NULL };
  PyObject *vids_obj = Py_None, *weights_obj = Py_None;
  PyObject *knn_obj, *knnk_obj;
  igraph_vector_t *weights = 0;
  igraph_vector_t knn, knnk;
  igraph_vs_t vids;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &vids_obj, &weights_obj)) {
    return NULL;
  }

  if (igraph_vector_init(&knn, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&knnk, 0)) {
    igraph_vector_destroy(&knn);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_PyObject_to_vs_t(vids_obj, &vids, &self->g, 0, 0)) {
    igraph_vector_destroy(&knn);
    igraph_vector_destroy(&knnk);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_obj, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) {
    igraph_vs_destroy(&vids);
    igraph_vector_destroy(&knn);
    igraph_vector_destroy(&knnk);
    return NULL;
  }

  if (igraph_avg_nearest_neighbor_degree(&self->g, vids, &knn, &knnk, weights)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vids);
    igraph_vector_destroy(&knn);
    igraph_vector_destroy(&knnk);
    if (weights) {
      igraph_vector_destroy(weights);
      free(weights);
    }
    return NULL;
  }

  igraph_vs_destroy(&vids);
  if (weights) {
    igraph_vector_destroy(weights);
    free(weights);
  }

  knn_obj = igraphmodule_vector_t_to_PyList(&knn, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&knn);
  if (!knn_obj) {
    igraph_vector_destroy(&knnk);
    return NULL;
  }

  knnk_obj = igraphmodule_vector_t_to_PyList(&knnk, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&knnk);
  if (!knnk_obj) {
    Py_DECREF(knn_obj);
    return NULL;
  }

  return Py_BuildValue("NN", knn_obj, knnk_obj);
}

/** \ingroup python_interface_graph
 * \brief Calculates the radius of an \c igraph.Graph
 * 
 * \return the radius as a Python integer
 * \sa igraph_radius
 */
PyObject *igraphmodule_Graph_radius(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  PyObject *mode_o = Py_None;
  igraph_neimode_t mode = IGRAPH_OUT;
  igraph_real_t radius;

  static char *kwlist[] = { "mode", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist,
                                   &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (igraph_radius(&self->g, &radius, mode)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  return PyFloat_FromDouble((double)radius);
}

/**********************************************************************
 * Deterministic and non-deterministic graph generators               *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Generates a graph from its adjacency matrix
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_adjacency
 */
PyObject *igraphmodule_Graph_Adjacency(PyTypeObject * type,
                                       PyObject * args, PyObject * kwds) {
  igraphmodule_GraphObject *self;
  igraph_t g;
  igraph_matrix_t m;
  PyObject *matrix, *mode_o = Py_None;
  igraph_adjacency_t mode = IGRAPH_ADJ_DIRECTED;

  static char *kwlist[] = { "matrix", "mode", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O", kwlist,
                                   &PyList_Type, &matrix, &mode_o))
    return NULL;
  if (igraphmodule_PyObject_to_adjacency_t(mode_o, &mode)) return NULL;

  if (igraphmodule_PyList_to_matrix_t(matrix, &m)) {
    PyErr_SetString(PyExc_TypeError,
                    "Error while converting adjacency matrix");
    return NULL;
  }

  if (igraph_adjacency(&g, &m, mode)) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&m);
    return NULL;
  }

  igraph_matrix_destroy(&m);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph from the Graph Atlas
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_atlas
 */
PyObject *igraphmodule_Graph_Atlas(PyTypeObject * type, PyObject * args)
{
  long n;
  igraphmodule_GraphObject *self;
  igraph_t g;

  if (!PyArg_ParseTuple(args, "l", &n))
    return NULL;

  if (igraph_atlas(&g, (igraph_integer_t) n)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on the Barabasi-Albert model
 * This is intended to be a class method in Python, so the first argument
 * is the type object and not the Python igraph object (because we have
 * to allocate that in this method).
 * 
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_barabasi_game
 */
PyObject *igraphmodule_Graph_Barabasi(PyTypeObject * type,
                                      PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n, m = 1;
  float power = 1.0f, zero_appeal = 1.0f;
  igraph_vector_t outseq;
  igraph_t *start_from = 0;
  igraph_barabasi_algorithm_t algo = IGRAPH_BARABASI_PSUMTREE;
  PyObject *m_obj = 0, *outpref = Py_False, *directed = Py_False;
  PyObject *implementation_o = Py_None;
  PyObject *start_from_o = Py_None;

  static char *kwlist[] =
    { "n", "m", "outpref", "directed", "power", "zero_appeal",
      "implementation", "start_from", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|OOOffOO", kwlist,
                                   &n, &m_obj, &outpref, &directed, &power,
                                   &zero_appeal, &implementation_o,
                                   &start_from_o))
    return NULL;

  if (igraphmodule_PyObject_to_barabasi_algorithm_t(implementation_o,
        &algo))
    return NULL;

  if (igraphmodule_PyObject_to_igraph_t(start_from_o, &start_from))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (m_obj == 0) {
    igraph_vector_init(&outseq, 0);
    m = 1;
  } else if (m_obj != 0) {
    /* let's check whether we have a constant out-degree or a list */
    if (PyInt_Check(m_obj)) {
      m = PyInt_AsLong(m_obj);
      igraph_vector_init(&outseq, 0);
    } else if (PyList_Check(m_obj)) {
      if (igraphmodule_PyObject_to_vector_t(m_obj, &outseq, 1, 0)) {
        /* something bad happened during conversion */
       return NULL;
      }
    } else {
      PyErr_SetString(PyExc_TypeError, "m must be an integer or a list of integers");
      return NULL;
    }
  }

  if (igraph_barabasi_game(&g, (igraph_integer_t) n,
                           (igraph_real_t) power,
                           (igraph_integer_t) m,
                           &outseq, PyObject_IsTrue(outpref),
                           (igraph_real_t) zero_appeal,
                           PyObject_IsTrue(directed), algo,
                           start_from)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&outseq);
    return NULL;
  }

  igraph_vector_destroy(&outseq);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a bipartite graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_barabasi_game
 */
PyObject *igraphmodule_Graph_Bipartite(PyTypeObject * type,
                                       PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  igraph_vector_bool_t types;
  igraph_vector_t edges;
  PyObject *types_o, *edges_o, *directed = Py_False;

  static char *kwlist[] = { "types", "edges", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|O", kwlist,
                                   &types_o, &edges_o, &directed))
    return NULL;

  if (igraphmodule_PyObject_to_vector_bool_t(types_o, &types))
    return NULL;

  if (igraphmodule_PyObject_to_edgelist(edges_o, &edges, 0)) {
    igraph_vector_bool_destroy(&types);
    return NULL;
  }

  if (igraph_create_bipartite(&g, &types, &edges, PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&edges);
    igraph_vector_bool_destroy(&types);
    return NULL;
  }

  igraph_vector_destroy(&edges);
  igraph_vector_bool_destroy(&types);
  
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a De Bruijn graph
 * \sa igraph_kautz
 */
PyObject *igraphmodule_Graph_De_Bruijn(PyTypeObject *type, PyObject *args,
  PyObject *kwds) {
  long int m, n;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = {"m", "n", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll", kwlist, &m, &n))
    return NULL;

  if (igraph_de_bruijn(&g, (igraph_integer_t) m, (igraph_integer_t) n)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject*)self;
}

/** \ingroup python_interface_graph
 * \brief Generates a random graph with a given degree sequence
 * This is intended to be a class method in Python, so the first argument
 * is the type object and not the Python igraph object (because we have
 * to allocate that in this method).
 * 
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_degree_sequence_game
 */
PyObject *igraphmodule_Graph_Degree_Sequence(PyTypeObject * type,
                                             PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  igraph_vector_t outseq, inseq;
  igraph_degseq_t meth = IGRAPH_DEGSEQ_SIMPLE;
  igraph_bool_t has_inseq = 0;
  PyObject *outdeg = NULL, *indeg = NULL, *method = NULL;

  static char *kwlist[] = { "out", "in", "method", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O!O", kwlist,
                                   &PyList_Type, &outdeg,
                                   &PyList_Type, &indeg,
								   &method))
    return NULL;

  if (igraphmodule_PyObject_to_degseq_t(method, &meth)) return NULL;
  if (igraphmodule_PyObject_to_vector_t(outdeg, &outseq, 1, 0)) return NULL;
  if (indeg) {
    if (igraphmodule_PyObject_to_vector_t(indeg, &inseq, 1, 0)) {
      igraph_vector_destroy(&outseq);
      return NULL;
    }
    has_inseq=1;
  }

  if (igraph_degree_sequence_game(&g, &outseq, has_inseq ? &inseq : 0, meth)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&outseq);
    if (has_inseq)
      igraph_vector_destroy(&inseq);
    return NULL;
  }

  igraph_vector_destroy(&outseq);
  if (has_inseq)
    igraph_vector_destroy(&inseq);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on the Erdõs-Rényi model
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_erdos_renyi_game
 */
PyObject *igraphmodule_Graph_Erdos_Renyi(PyTypeObject * type,
                                         PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n, m = -1;
  double p = -1.0;
  igraph_erdos_renyi_t t;
  PyObject *loops = NULL, *directed = NULL;

  static char *kwlist[] = { "n", "p", "m", "directed", "loops", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|dlO!O!", kwlist,
                                   &n, &p, &m,
                                   &PyBool_Type, &directed,
                                   &PyBool_Type, &loops))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (m == -1 && p == -1.0) {
    // no density parameters were given, throw exception
    PyErr_SetString(PyExc_TypeError, "Either m or p must be given.");
    return NULL;
  }
  if (m != -1 && p != -1.0) {
    // both density parameters were given, throw exception
    PyErr_SetString(PyExc_TypeError, "Only one must be given from m and p.");
    return NULL;
  }

  t = (m == -1) ? IGRAPH_ERDOS_RENYI_GNP : IGRAPH_ERDOS_RENYI_GNM;

  if (t == IGRAPH_ERDOS_RENYI_GNP) {
    if (p < 0.0 || p > 1.0) {
      // Invalid probability was given, throw exception
      PyErr_SetString(PyExc_ValueError, "p must be between 0 and 1.");
      return NULL;
    }
  } else {
    if (m < 0 || ((double)m)/n > n ) {
      // Invalid edge count was given, throw exception
      PyErr_SetString(PyExc_ValueError, "m must be between 0 and n^2.");
      return NULL;
    }
  }

  if (igraph_erdos_renyi_game(&g, t, (igraph_integer_t) n,
                              (igraph_real_t) ((t == IGRAPH_ERDOS_RENYI_GNM) ? m : p),
							  (directed == Py_True), (loops == Py_True))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on a simple growing model with vertex types
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_establishment_game
 */
PyObject *igraphmodule_Graph_Establishment(PyTypeObject * type,
                                           PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n, types, k;
  PyObject *type_dist, *pref_matrix;
  PyObject *directed = Py_False;
  igraph_matrix_t pm;
  igraph_vector_t td;

  char *kwlist[] = { "n", "k", "type_dist", "pref_matrix", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "llO!O!|O", kwlist,
                                   &n, &k, &PyList_Type, &type_dist,
                                   &PyList_Type, &pref_matrix, &directed))
    return NULL;

  if (n <= 0 || k <= 0) {
    PyErr_SetString(PyExc_ValueError,
                    "Number of vertices and the amount of connection trials per step must be positive.");
    return NULL;
  }
  types = PyList_Size(type_dist);

  if (igraphmodule_PyList_to_matrix_t(pref_matrix, &pm)) {
    PyErr_SetString(PyExc_TypeError,
                    "Error while converting preference matrix");
    return NULL;
  }
  if (igraph_matrix_nrow(&pm) != igraph_matrix_ncol(&pm) ||
      igraph_matrix_nrow(&pm) != types) {
    PyErr_SetString(PyExc_ValueError,
                    "Preference matrix must have exactly the same rows and columns as the number of types");
    igraph_matrix_destroy(&pm);
    return NULL;
  }
  if (igraphmodule_PyObject_to_vector_t(type_dist, &td, 1, 0)) {
    PyErr_SetString(PyExc_ValueError,
                    "Error while converting type distribution vector");
    igraph_matrix_destroy(&pm);
    return NULL;
  }

  if (igraph_establishment_game(&g, (igraph_integer_t) n,
                                (igraph_integer_t) types,
                                (igraph_integer_t) k, &td, &pm,
                                PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&pm);
    igraph_vector_destroy(&td);
    return NULL;
  }

  igraph_matrix_destroy(&pm);
  igraph_vector_destroy(&td);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a famous graph by name
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_famous
 */
PyObject *igraphmodule_Graph_Famous(PyTypeObject * type,
                                    PyObject * args, PyObject * kwds) {
  igraphmodule_GraphObject *self;
  igraph_t g;
  const char* name;

  static char *kwlist[] = { "name", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &name))
    return NULL;

  if (igraph_famous(&g, name)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}


/** \ingroup python_interface_graph
 * \brief Generates a graph based on the forest fire model
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_forest_fire_game
 */
PyObject *igraphmodule_Graph_Forest_Fire(PyTypeObject * type,
  PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n, ambs=1;
  double fw_prob, bw_factor=0.0;
  PyObject *directed = Py_False;

  static char *kwlist[] = {"n", "fw_prob", "bw_factor", "ambs", "directed", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ld|dlO", kwlist,
                                   &n, &fw_prob, &bw_factor, &ambs, &directed))
    return NULL;

  if (igraph_forest_fire_game(&g, (igraph_integer_t)n,
                              (igraph_real_t)fw_prob, (igraph_real_t)bw_factor,
                              (igraph_integer_t)ambs,
							  (igraph_bool_t)(PyObject_IsTrue(directed)))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}


/** \ingroup python_interface_graph
 * \brief Generates a full graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_full
 */
PyObject *igraphmodule_Graph_Full(PyTypeObject * type,
                                  PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n;
  PyObject *loops = Py_False, *directed = Py_False;

  char *kwlist[] = { "n", "directed", "loops", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|OO", kwlist, &n,
                                   &directed, &loops))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (igraph_full(&g, (igraph_integer_t) n, PyObject_IsTrue(directed),
			      PyObject_IsTrue(loops))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a full bipartite graph
 * \sa igraph_full_bipartite
 */
PyObject *igraphmodule_Graph_Full_Bipartite(PyTypeObject * type,
                                            PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  igraph_vector_bool_t vertex_types;
  igraph_neimode_t mode = IGRAPH_ALL;
  long int n1, n2;
  PyObject *mode_o = Py_None, *directed = Py_False, *vertex_types_o = 0;

  static char *kwlist[] = { "n1", "n2", "directed", "mode", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|OO", kwlist, &n1, &n2,
                                   &directed, &mode_o))
    return NULL;

  if (n1 < 0 || n2 < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
	  return NULL;

  if (igraph_vector_bool_init(&vertex_types, n1+n2)) {
    igraphmodule_handle_igraph_error();
	return NULL;
  }

  if (igraph_full_bipartite(&g, &vertex_types,
        (igraph_integer_t) n1, (igraph_integer_t) n2,
        PyObject_IsTrue(directed), mode)) {
	igraph_vector_bool_destroy(&vertex_types);
    igraphmodule_handle_igraph_error();
	return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  vertex_types_o = igraphmodule_vector_bool_t_to_PyList(&vertex_types);
  igraph_vector_bool_destroy(&vertex_types);
  if (vertex_types_o == 0) return NULL;
  return Py_BuildValue("NN", (PyObject *) self, vertex_types_o);
}

/** \ingroup python_interface_graph
 * \brief Generates a full citation graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_full
 */
PyObject *igraphmodule_Graph_Full_Citation(PyTypeObject *type,
    PyObject *args, PyObject *kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n;
  PyObject *directed = Py_False;

  char *kwlist[] = { "n", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|O", kwlist, &n, &directed))
    return NULL;

  if (igraph_full_citation(&g, (igraph_integer_t) n,
                           (igraph_bool_t) PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on the geometric random model
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_grg_game
 */
PyObject *igraphmodule_Graph_GRG(PyTypeObject * type,
                                 PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n;
  double r;
  PyObject *torus = Py_False;
  PyObject *o_xs, *o_ys;
  igraph_vector_t xs, ys;

  static char *kwlist[] = { "n", "radius", "torus", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ld|O", kwlist,
                                   &n, &r, &torus))
    return NULL;

  if (igraph_vector_init(&xs, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  } else if (igraph_vector_init(&ys, 0)) {
    igraph_vector_destroy(&xs);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_grg_game(&g, (igraph_integer_t) n, (igraph_real_t) r,
                      PyObject_IsTrue(torus), &xs, &ys)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&xs);
    igraph_vector_destroy(&ys);
    return NULL;
  }

  o_xs = igraphmodule_vector_t_to_PyList(&xs, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&xs);
  if (!o_xs) {
    igraph_destroy(&g);
    igraph_vector_destroy(&ys);
    return NULL;
  }
  o_ys = igraphmodule_vector_t_to_PyList(&ys, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&ys);
  if (!o_ys) {
    igraph_destroy(&g);
    Py_DECREF(o_xs);
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);
  return Py_BuildValue("NNN", (PyObject*)self, o_xs, o_ys);
}

/** \ingroup python_interface_graph
 * \brief Generates a growing random graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_growing_random_game
 */
PyObject *igraphmodule_Graph_Growing_Random(PyTypeObject * type,
                                            PyObject * args, PyObject * kwds)
{
  long n, m;
  PyObject *directed = NULL, *citation = NULL;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "n", "m", "directed", "citation", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|O!O!", kwlist, &n, &m,
                                   &PyBool_Type, &directed,
                                   &PyBool_Type, &citation))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (m < 0) {
    PyErr_SetString(PyExc_ValueError,
                    "Number of new edges per iteration must be positive.");
    return NULL;
  }

  if (igraph_growing_random_game(&g, (igraph_integer_t) n,
                                 (igraph_integer_t) m,
                                 (directed == Py_True),
                                 (citation == Py_True))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a bipartite graph from an incidence matrix
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_incidence
 */
PyObject *igraphmodule_Graph_Incidence(PyTypeObject * type,
                                       PyObject * args, PyObject * kwds) {
  igraphmodule_GraphObject *self;
  igraph_matrix_t matrix;
  igraph_vector_bool_t vertex_types;
  igraph_t g;
  PyObject *matrix_o, *vertex_types_o;
  PyObject *mode_o = Py_None, *directed = Py_False, *multiple = Py_False;
  igraph_neimode_t mode = IGRAPH_OUT;

  static char *kwlist[] = { "matrix", "directed", "mode", "multiple", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|OOO", kwlist, &PyList_Type, &matrix_o,
			  &directed, &mode_o, &multiple))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;

  if (igraph_vector_bool_init(&vertex_types, 0)) {
    igraphmodule_handle_igraph_error();
	return NULL;
  }

  if (igraphmodule_PyList_to_matrix_t(matrix_o, &matrix)) {
	igraph_vector_bool_destroy(&vertex_types);
    PyErr_SetString(PyExc_TypeError,
                    "Error while converting incidence matrix");
    return NULL;
  }

  if (igraph_incidence(&g, &vertex_types, &matrix,
			  PyObject_IsTrue(directed), mode, PyObject_IsTrue(multiple))) {
	igraphmodule_handle_igraph_error();
	igraph_matrix_destroy(&matrix);
	igraph_vector_bool_destroy(&vertex_types);
	return NULL;
  }

  igraph_matrix_destroy(&matrix);
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  vertex_types_o = igraphmodule_vector_bool_t_to_PyList(&vertex_types);
  igraph_vector_bool_destroy(&vertex_types);
  if (vertex_types_o == 0) return NULL;
  return Py_BuildValue("NN", (PyObject *) self, vertex_types_o);
}

/** \ingroup python_interface_graph
 * \brief Generates a graph with a given isomorphy class
 * This is intended to be a class method in Python, so the first argument
 * is the type object and not the Python igraph object (because we have
 * to allocate that in this method).
 * 
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_isoclass_create
 */
PyObject *igraphmodule_Graph_Isoclass(PyTypeObject * type,
                                      PyObject * args, PyObject * kwds)
{
  long int n, isoclass;
  PyObject *directed = Py_False;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "n", "class", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|O", kwlist,
                                   &n, &isoclass, &directed))
    return NULL;

  if (n < 3 || n > 4) {
    PyErr_SetString(PyExc_ValueError,
                    "Only graphs with 3 or 4 vertices are supported");
    return NULL;
  }

  if (igraph_isoclass_create(&g, (igraph_integer_t) n,
        (igraph_integer_t) isoclass, PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a Kautz graph
 * \sa igraph_kautz
 */
PyObject *igraphmodule_Graph_Kautz(PyTypeObject *type, PyObject *args,
  PyObject *kwds) {
  long int m, n;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = {"m", "n", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll", kwlist, &m, &n))
    return NULL;

  if (igraph_kautz(&g, (igraph_integer_t) m, (igraph_integer_t) n)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject*)self;
}

/** \ingroup python_interface_graph
 * \brief Generates a k-regular random graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_k_regular_game
 */
PyObject *igraphmodule_Graph_K_Regular(PyTypeObject * type,
                                 PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long int n, k;
  PyObject *directed_o = Py_False, *multiple_o = Py_False;

  static char *kwlist[] = { "n", "k", "directed", "multiple", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|OO", kwlist,
                                   &n, &k, &directed_o, &multiple_o))
    return NULL;

  if (igraph_k_regular_game(&g, (igraph_integer_t) n, (igraph_integer_t) k,
                      PyObject_IsTrue(directed_o), PyObject_IsTrue(multiple_o))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject*)self;
}

/** \ingroup python_interface_graph
 * \brief Generates a regular lattice
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_lattice
 */
PyObject *igraphmodule_Graph_Lattice(PyTypeObject * type,
                                     PyObject * args, PyObject * kwds)
{
  igraph_vector_t dimvector;
  long int nei = 1;
  igraph_bool_t directed;
  igraph_bool_t mutual;
  igraph_bool_t circular;
  PyObject *o_directed = Py_False, *o_mutual = Py_True, *o_circular = Py_True;
  PyObject *o_dimvector = Py_None;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "dim", "nei", "directed", "mutual", "circular", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|lOOO", kwlist,
                                   &PyList_Type, &o_dimvector,
                                   &nei, &o_directed, &o_mutual, &o_circular))
    return NULL;

  directed = PyObject_IsTrue(o_directed);
  mutual = PyObject_IsTrue(o_mutual);
  circular = PyObject_IsTrue(o_circular);

  if (igraphmodule_PyObject_to_vector_t(o_dimvector, &dimvector, 1, 0))
    return NULL;

  if (igraph_lattice(&g, &dimvector, (igraph_integer_t) nei,
        directed, mutual, circular)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&dimvector);
    return NULL;
  }

  igraph_vector_destroy(&dimvector);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a 3-regular Hamiltonian graph from LCF notation
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_lattice
 */
PyObject *igraphmodule_Graph_LCF(PyTypeObject *type,
                                 PyObject *args, PyObject *kwds) {
  igraph_vector_t shifts;
  long int repeats, n;
  PyObject *o_shifts;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "n", "shifts", "repeats", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "lOl", kwlist,
                                   &n, &o_shifts, &repeats))
    return NULL;

  if (igraphmodule_PyObject_to_vector_t(o_shifts, &shifts, 0, 0))
    return NULL;

  if (igraph_lcf_vector(&g, (igraph_integer_t) n, &shifts, (igraph_integer_t) repeats)) {
    igraph_vector_destroy(&shifts);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  igraph_vector_destroy(&shifts);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on vertex types and connection preferences
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_preference_game
 */
PyObject *igraphmodule_Graph_Preference(PyTypeObject * type,
                                        PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n, types;
  PyObject *type_dist, *pref_matrix;
  PyObject *directed = Py_False;
  PyObject *loops = Py_False;
  igraph_matrix_t pm;
  igraph_vector_t td;
  igraph_vector_t type_vec;
  PyObject *type_vec_o;
  PyObject *attribute_key = Py_None;
  igraph_bool_t store_attribs;

  char *kwlist[] =
    { "n", "type_dist", "pref_matrix", "attribute", "directed", "loops",
NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "lO!O!|OOO", kwlist,
                                   &n, &PyList_Type, &type_dist,
                                   &PyList_Type, &pref_matrix,
                                   &attribute_key, &directed, &loops))
    return NULL;

  types = PyList_Size(type_dist);

  if (igraphmodule_PyList_to_matrix_t(pref_matrix, &pm)) return NULL;
  if (igraphmodule_PyObject_float_to_vector_t(type_dist, &td)) {
    igraph_matrix_destroy(&pm);
    return NULL;
  }

  store_attribs = (attribute_key && attribute_key != Py_None);
  if (store_attribs && igraph_vector_init(&type_vec, (igraph_integer_t) n)) {
    igraph_matrix_destroy(&pm);
    igraph_vector_destroy(&td);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_preference_game(&g, (igraph_integer_t) n,
                             (igraph_integer_t) types, &td, 0, &pm,
                             store_attribs ? &type_vec : 0,
                             PyObject_IsTrue(directed),
                             PyObject_IsTrue(loops))) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&pm);
    igraph_vector_destroy(&td);
    if (store_attribs)
      igraph_vector_destroy(&type_vec);
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  if (store_attribs) {
    type_vec_o = igraphmodule_vector_t_to_PyList(&type_vec, IGRAPHMODULE_TYPE_INT);
    if (type_vec_o == 0) {
      igraph_matrix_destroy(&pm);
      igraph_vector_destroy(&td);
      igraph_vector_destroy(&type_vec);
      Py_DECREF(self);
      return NULL;
    }
    if (attribute_key != Py_None && attribute_key != 0) {
      if (PyDict_SetItem(ATTR_STRUCT_DICT(&self->g)[ATTRHASH_IDX_VERTEX],
                         attribute_key, type_vec_o) == -1) {
        Py_DECREF(type_vec_o);
        igraph_matrix_destroy(&pm);
        igraph_vector_destroy(&td);
        igraph_vector_destroy(&type_vec);
        Py_DECREF(self);
        return NULL;
      }
    }

    Py_DECREF(type_vec_o);
    igraph_vector_destroy(&type_vec);
  }

  igraph_matrix_destroy(&pm);
  igraph_vector_destroy(&td);
  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on asymmetric vertex types and connection preferences
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_asymmetric_preference_game
 */
PyObject *igraphmodule_Graph_Asymmetric_Preference(PyTypeObject * type,
                                                   PyObject * args,
                                                   PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n, types;
  PyObject *type_dist_matrix, *pref_matrix;
  PyObject *loops = Py_False;
  igraph_matrix_t pm;
  igraph_matrix_t td;
  igraph_vector_t in_type_vec, out_type_vec;
  PyObject *type_vec_o;
  PyObject *attribute_key = Py_None;
  igraph_bool_t store_attribs;

  char *kwlist[] =
    { "n", "type_dist_matrix", "pref_matrix", "attribute", "loops", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "lO!O!|OO", kwlist,
                                   &n, &PyList_Type, &type_dist_matrix,
                                   &PyList_Type, &pref_matrix,
                                   &attribute_key, &loops))
    return NULL;

  types = PyList_Size(type_dist_matrix);
  if (igraphmodule_PyList_to_matrix_t(pref_matrix, &pm)) return NULL;
  if (igraphmodule_PyList_to_matrix_t(type_dist_matrix, &td)) {
    igraph_matrix_destroy(&pm);
    return NULL;
  }

  store_attribs = (attribute_key && attribute_key != Py_None);
  if (store_attribs) {
    if (igraph_vector_init(&in_type_vec, (igraph_integer_t) n)) {
      igraph_matrix_destroy(&pm);
      igraph_matrix_destroy(&td);
      igraphmodule_handle_igraph_error();
      return NULL;
    }
    if (igraph_vector_init(&out_type_vec, (igraph_integer_t) n)) {
      igraph_matrix_destroy(&pm);
      igraph_matrix_destroy(&td);
      igraph_vector_destroy(&in_type_vec);
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  if (igraph_asymmetric_preference_game(&g, (igraph_integer_t) n,
                                        (igraph_integer_t) types, &td, &pm,
                                        store_attribs ? &in_type_vec : 0,
                                        store_attribs ? &out_type_vec : 0,
                                        PyObject_IsTrue(loops))) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&in_type_vec);
    igraph_vector_destroy(&out_type_vec);
    igraph_matrix_destroy(&pm);
    igraph_matrix_destroy(&td);
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  if (store_attribs) {
    type_vec_o = igraphmodule_vector_t_pair_to_PyList(&in_type_vec,
                                                      &out_type_vec);
    if (type_vec_o == NULL) {
      igraph_matrix_destroy(&pm);
      igraph_matrix_destroy(&td);
      igraph_vector_destroy(&in_type_vec);
      igraph_vector_destroy(&out_type_vec);
      Py_DECREF(self);
      return NULL;
    }
    if (attribute_key != Py_None && attribute_key != 0) {
      if (PyDict_SetItem(ATTR_STRUCT_DICT(&self->g)[ATTRHASH_IDX_VERTEX],
                         attribute_key, type_vec_o) == -1) {
        Py_DECREF(type_vec_o);
        igraph_matrix_destroy(&pm);
        igraph_matrix_destroy(&td);
        igraph_vector_destroy(&in_type_vec);
        igraph_vector_destroy(&out_type_vec);
        Py_DECREF(self);
        return NULL;
      }
    }

    Py_DECREF(type_vec_o);
    igraph_vector_destroy(&in_type_vec);
    igraph_vector_destroy(&out_type_vec);
  }

  igraph_matrix_destroy(&pm);
  igraph_matrix_destroy(&td);
  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on sort of a "windowed" Barabasi-Albert model
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_recent_degree_game
 */
PyObject *igraphmodule_Graph_Recent_Degree(PyTypeObject * type,
                                           PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  long n, m = 0, window = 0;
  float power = 0.0f, zero_appeal = 0.0f;
  igraph_vector_t outseq;
  PyObject *m_obj, *outpref = Py_False, *directed = Py_False;

  char *kwlist[] =
    { "n", "m", "window", "outpref", "directed", "power", "zero_appeal",
NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "lOl|OOff", kwlist,
                                   &n, &m_obj, &window, &outpref, &directed,
                                   &power, &zero_appeal))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  // let's check whether we have a constant out-degree or a list
  if (PyInt_Check(m_obj)) {
    m = PyInt_AsLong(m_obj);
    igraph_vector_init(&outseq, 0);
  }
  else if (PyList_Check(m_obj)) {
    if (igraphmodule_PyObject_to_vector_t(m_obj, &outseq, 1, 0)) {
      // something bad happened during conversion
      return NULL;
    }
  }

  if (igraph_recent_degree_game(&g, (igraph_integer_t) n,
                                (igraph_real_t) power,
                                (igraph_integer_t) window,
                                (igraph_integer_t) m, &outseq,
                                PyObject_IsTrue(outpref),
                                (igraph_real_t) zero_appeal,
                                PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&outseq);
    return NULL;
  }

  igraph_vector_destroy(&outseq);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a ring-shaped graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_ring
 */
PyObject *igraphmodule_Graph_Ring(PyTypeObject * type,
                                  PyObject * args, PyObject * kwds)
{
  long n;
  PyObject *directed = Py_False, *mutual = Py_False, *circular = Py_True;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "n", "directed", "mutual", "circular", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|O!O!O!", kwlist, &n,
                                   &PyBool_Type, &directed,
                                   &PyBool_Type, &mutual,
                                   &PyBool_Type, &circular))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (igraph_ring(&g, (igraph_integer_t) n, (directed == Py_True),
                  (mutual == Py_True), (circular == Py_True))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a star graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_star
 */
PyObject *igraphmodule_Graph_Star(PyTypeObject * type,
                                  PyObject * args, PyObject * kwds)
{
  long n, center = 0;
  igraph_star_mode_t mode = IGRAPH_STAR_UNDIRECTED;
  PyObject* mode_o = Py_None;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "n", "mode", "center", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|Ol", kwlist,
                                   &n, &mode_o, &center))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (center >= n || center < 0) {
    PyErr_SetString(PyExc_ValueError,
                    "Central vertex ID should be between 0 and n-1");
    return NULL;
  }

  if (igraphmodule_PyObject_to_star_mode_t(mode_o, &mode)) {
    PyErr_SetString(PyExc_ValueError,
                    "Mode should be either \"in\", \"out\", \"mutual\" or \"undirected\"");
    return NULL;
  }

  if (igraph_star(&g, (igraph_integer_t) n, mode, (igraph_integer_t) center)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a non-growing random graph with edge probabilities
 *        proportional to node fitnesses.
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_static_fitness_game
 */
PyObject *igraphmodule_Graph_Static_Fitness(PyTypeObject *type,
    PyObject* args, PyObject* kwds) {
  igraphmodule_GraphObject *self;
  igraph_t g;
  long int m;
  PyObject *fitness_out_o = Py_None, *fitness_in_o = Py_None;
  PyObject *fitness_o = Py_None;
  PyObject *multiple = Py_False, *loops = Py_False;
  igraph_vector_t fitness_out, fitness_in;

  static char *kwlist[] = { "m", "fitness_out", "fitness_in",
    "loops", "multiple", "fitness", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|OOOOO", kwlist,
                                   &m, &fitness_out_o, &fitness_in_o,
                                   &loops, &multiple, &fitness_o))
    return NULL;

  /* This trickery allows us to use "fitness" or "fitness_out" as
   * keyword argument, with "fitness_out" taking precedence over
   * "fitness" */
  if (fitness_out_o == Py_None)
    fitness_out_o = fitness_o;
  if (fitness_out_o == Py_None) {
    PyErr_SetString(PyExc_TypeError,
        "Required argument 'fitness_out' (pos 2) not found");
    return NULL;
  }

  if (igraphmodule_PyObject_float_to_vector_t(fitness_out_o, &fitness_out))
    return NULL;

  if (fitness_in_o != Py_None) {
    if (igraphmodule_PyObject_float_to_vector_t(fitness_in_o, &fitness_in)) {
      igraph_vector_destroy(&fitness_out);
      return NULL;
    }
  }

  if (igraph_static_fitness_game(&g, (igraph_integer_t) m, &fitness_out,
        fitness_in_o == Py_None ? 0 : &fitness_in,
        PyObject_IsTrue(loops), PyObject_IsTrue(multiple))) {
    igraph_vector_destroy(&fitness_out);
    if (fitness_in_o != Py_None)
      igraph_vector_destroy(&fitness_in);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  igraph_vector_destroy(&fitness_out);
  if (fitness_in_o != Py_None)
    igraph_vector_destroy(&fitness_in);

  CREATE_GRAPH_FROM_TYPE(self, g, type);
  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a non-growing random graph with prescribed power-law
 *        degree distributions.
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_static_power_law_game
 */
PyObject *igraphmodule_Graph_Static_Power_Law(PyTypeObject *type,
    PyObject* args, PyObject* kwds) {
  igraphmodule_GraphObject *self;
  igraph_t g;
  long int n, m;
  float exponent_out = -1.0f, exponent_in = -1.0f, exponent = -1.0f;
  PyObject *multiple = Py_False, *loops = Py_False;
  PyObject *finite_size_correction = Py_True;

  static char *kwlist[] = { "n", "m", "exponent_out", "exponent_in",
    "loops", "multiple", "finite_size_correction", "exponent", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|ffOOOf", kwlist,
                                   &n, &m, &exponent_out, &exponent_in,
                                   &loops, &multiple, &finite_size_correction,
                                   &exponent))
    return NULL;

  /* This trickery allows us to use "exponent" or "exponent_out" as
   * keyword argument, with "exponent_out" taking precedence over
   * "exponent" */
  if (exponent_out == -1.0)
    exponent_out = exponent;
  if (exponent_out == -1.0) {
    PyErr_SetString(PyExc_TypeError,
        "Required argument 'exponent_out' (pos 3) not found");
    return NULL;
  }

  if (igraph_static_power_law_game(&g, (igraph_integer_t) n, (igraph_integer_t) m,
        exponent_out, exponent_in,
        PyObject_IsTrue(loops), PyObject_IsTrue(multiple),
        PyObject_IsTrue(finite_size_correction))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);
  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a tree graph where almost all vertices have an equal number of children
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_tree
 */
PyObject *igraphmodule_Graph_Tree(PyTypeObject * type,
                                  PyObject * args, PyObject * kwds)
{
  long int n, children;
  PyObject *tree_mode_o = Py_None, *tree_type_o = Py_None;
  igraph_tree_mode_t mode = IGRAPH_TREE_UNDIRECTED;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "n", "children", "mode", "type", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|OO", kwlist,
                                   &n, &children,
                                   &tree_mode_o, &tree_type_o))
    return NULL;

  if (n < 0) {
    PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
    return NULL;
  }

  if (tree_mode_o == Py_None && tree_type_o != Py_None) {
    tree_mode_o = tree_type_o;
    PY_IGRAPH_DEPRECATED("type=... keyword argument is deprecated since igraph 0.6, use mode=... instead");
  }

  if (igraphmodule_PyObject_to_tree_mode_t(tree_mode_o, &mode)) {
    return NULL;
  }

  if (igraph_tree(&g, (igraph_integer_t) n, (igraph_integer_t) children, mode)) {
      igraphmodule_handle_igraph_error();
      return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph based on the Watts-Strogatz model
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_watts_strogatz_game
 */
PyObject *igraphmodule_Graph_Watts_Strogatz(PyTypeObject * type,
                                            PyObject * args, PyObject * kwds)
{
  long int nei = 1, dim, size;
  double p;
  PyObject* loops = Py_False;
  PyObject* multiple = Py_False;
  igraphmodule_GraphObject *self;
  igraph_t g;

  static char *kwlist[] = { "dim", "size", "nei", "p", "loops", "multiple", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "llld|OO", kwlist,
                                   &dim, &size, &nei, &p, &loops, &multiple))
    return NULL;

  if (igraph_watts_strogatz_game(&g, (igraph_integer_t) dim,
        (igraph_integer_t) size, (igraph_integer_t) nei, p,
        PyObject_IsTrue(loops), PyObject_IsTrue(multiple))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Generates a graph from its weighted adjacency matrix
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_weighted_adjacency
 */
PyObject *igraphmodule_Graph_Weighted_Adjacency(PyTypeObject * type,
                                       PyObject * args, PyObject * kwds) {
  igraphmodule_GraphObject *self;
  igraph_t g;
  igraph_matrix_t m;
  PyObject *matrix, *mode_o = Py_None, *attr_o = Py_None, *s = 0;
  PyObject *loops = Py_True;
  char* attr = 0;
  igraph_adjacency_t mode = IGRAPH_ADJ_DIRECTED;

  static char *kwlist[] = { "matrix", "mode", "attr", "loops", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|OOO", kwlist,
                                   &PyList_Type, &matrix, &mode_o, &attr_o,
                                   &loops))
    return NULL;

  if (igraphmodule_PyObject_to_adjacency_t(mode_o, &mode))
    return NULL;

  if (attr_o != Py_None) {
    s = PyObject_Str(attr_o);
    if (s) {
      attr = PyString_CopyAsString(s);
      if (attr == 0)
        return NULL;
    } else return NULL;
  }

  if (igraphmodule_PyList_to_matrix_t(matrix, &m)) {
    if (attr != 0)
      free(attr);
    PyErr_SetString(PyExc_TypeError,
                    "Error while converting adjacency matrix");
    return NULL;
  }

  if (igraph_weighted_adjacency(&g, &m, mode, attr ? attr : "weight",
        PyObject_IsTrue(loops))) {
    igraphmodule_handle_igraph_error();
    if (attr != 0)
      free(attr);
    igraph_matrix_destroy(&m);
    return NULL;
  }

  if (attr != 0)
    free(attr);
  igraph_matrix_destroy(&m);

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/**********************************************************************
 * Advanced structural properties of graphs                           *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Calculates the articulation points of a graph.
 * \return the list of articulation points in a PyObject
 * \sa igraph_articulation_points
 */
PyObject *igraphmodule_Graph_articulation_points(igraphmodule_GraphObject *self) {
  igraph_vector_t res;
  PyObject *o;
  if (igraph_vector_init(&res, 0)) {
	igraphmodule_handle_igraph_error();
	return NULL;
  }

  if (igraph_articulation_points(&self->g, &res)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&res);
    return NULL;
  }

  igraph_vector_sort(&res);
  o = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&res);
  return o;
}

/** \ingroup python_interface_graph
 * \brief Calculates the nominal assortativity coefficient
 * \sa igraph_assortativity_nominal
 */
PyObject *igraphmodule_Graph_assortativity_nominal(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = { "types", "directed", NULL };
  PyObject *types_o = Py_None, *directed = Py_True;
  igraph_real_t res;
  int ret;
  igraph_vector_t *types = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &types_o, &directed))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(types_o, self, &types, ATTRIBUTE_TYPE_VERTEX))
    return NULL;

  ret = igraph_assortativity_nominal(&self->g, types, &res, PyObject_IsTrue(directed));

  if (types) {
    igraph_vector_destroy(types); free(types);
  }

  if (ret) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  return Py_BuildValue("d", (double)(res));
}

/** \ingroup python_interface_graph
 * \brief Calculates the assortativity coefficient
 * \sa igraph_assortativity
 */
PyObject *igraphmodule_Graph_assortativity(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = { "types1", "types2", "directed", NULL };
  PyObject *types1_o = Py_None, *types2_o = Py_None, *directed = Py_True;
  igraph_real_t res;
  int ret;
  igraph_vector_t *types1 = 0, *types2 = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OO", kwlist, &types1_o, &types2_o, &directed))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(types1_o, self, &types1, ATTRIBUTE_TYPE_VERTEX))
    return NULL;
  if (igraphmodule_attrib_to_vector_t(types2_o, self, &types2, ATTRIBUTE_TYPE_VERTEX)) {
    if (types1) { igraph_vector_destroy(types1); free(types1); }
    return NULL;
  }

  ret = igraph_assortativity(&self->g, types1, types2, &res, PyObject_IsTrue(directed));

  if (types1) { igraph_vector_destroy(types1); free(types1); }
  if (types2) { igraph_vector_destroy(types2); free(types2); }

  if (ret) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  return Py_BuildValue("d", (double)(res));
}

/** \ingroup python_interface_graph
 * \brief Calculates the assortativity coefficient for degrees
 * \sa igraph_assortativity_degree
 */
PyObject *igraphmodule_Graph_assortativity_degree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = { "directed", NULL };
  PyObject *directed = Py_True;
  igraph_real_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &directed))
    return NULL;

  if (igraph_assortativity_degree(&self->g, &res, PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  return Py_BuildValue("d", (double)(res));
}

/** \ingroup python_interface_graph
 * \brief Calculates Kleinberg's authority scores of the vertices in the graph
 * \sa igraph_authority_score
 */
PyObject *igraphmodule_Graph_authority_score(
  igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] =
    { "weights", "scale", "arpack_options", "return_eigenvalue", NULL };
  PyObject *scale_o = Py_True, *weights_o = Py_None;
  PyObject *arpack_options_o = igraphmodule_arpack_options_default;
  igraphmodule_ARPACKOptionsObject *arpack_options;
  PyObject *return_eigenvalue = Py_False;
  PyObject *res_o;
  igraph_real_t value;
  igraph_vector_t res, *weights = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO!O", kwlist, &weights_o,
	                               &scale_o, &igraphmodule_ARPACKOptionsType,
                                   &arpack_options_o, &return_eigenvalue))
    return NULL;

  if (igraph_vector_init(&res, 0)) return igraphmodule_handle_igraph_error();

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  arpack_options = (igraphmodule_ARPACKOptionsObject*)arpack_options_o;
  if (igraph_authority_score(&self->g, &res, &value, PyObject_IsTrue(scale_o),
      weights, igraphmodule_ARPACKOptions_get(arpack_options))) {
    igraphmodule_handle_igraph_error();
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraph_vector_destroy(&res);
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  res_o = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT); 
  igraph_vector_destroy(&res);
  if (res_o == NULL) return igraphmodule_handle_igraph_error();

  if (PyObject_IsTrue(return_eigenvalue)) {
    PyObject *ev_o = PyFloat_FromDouble((double)value);
    if (ev_o == NULL) {
      Py_DECREF(res_o);
      return igraphmodule_handle_igraph_error();
    }
    return Py_BuildValue("NN", res_o, ev_o);
  }

  return res_o;
}


/** \ingroup python_interface_graph
 * \brief Calculates the average path length in a graph.
 * \return the average path length as a PyObject
 * \sa igraph_average_path_length
 */
PyObject *igraphmodule_Graph_average_path_length(igraphmodule_GraphObject *
                                                 self, PyObject * args,
                                                 PyObject * kwds)
{
  char *kwlist[] = { "directed", "unconn", NULL };
  PyObject *directed = Py_True, *unconn = Py_True;
  igraph_real_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!", kwlist,
                                   &PyBool_Type, &directed,
                                   &PyBool_Type, &unconn))
    return NULL;

  if (igraph_average_path_length(&self->g, &res, (directed == Py_True),
                                 (unconn == Py_True))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  return PyFloat_FromDouble(res);
}

/** \ingroup python_interface_graph
 * \brief Calculates the betweennesses of some vertices in a graph.
 * \return the betweennesses as a list (or a single float)
 * \sa igraph_betweenness
 */
PyObject *igraphmodule_Graph_betweenness(igraphmodule_GraphObject * self,
                                         PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "vertices", "directed", "cutoff", "weights",
    "nobigint", NULL };
  PyObject *directed = Py_True;
  PyObject *vobj = Py_None, *list;
  PyObject *cutoff = Py_None;
  PyObject *weights_o = Py_None;
  PyObject *nobigint = Py_True;
  igraph_vector_t res, *weights = 0;
  igraph_bool_t return_single = 0;
  igraph_vs_t vs;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOOO", kwlist,
                                   &vobj, &directed, &cutoff, &weights_o,
                                   &nobigint)) {
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&res, 0)) {
    igraph_vs_destroy(&vs);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return igraphmodule_handle_igraph_error();
  }

  if (cutoff == Py_None) {
    if (igraph_betweenness(&self->g, &res, vs, PyObject_IsTrue(directed),
          weights, PyObject_IsTrue(nobigint))) {
      igraph_vs_destroy(&vs);
      igraph_vector_destroy(&res);
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else if (PyNumber_Check(cutoff)) {
    PyObject *cutoff_num = PyNumber_Int(cutoff);
    if (cutoff_num == NULL) {
      igraph_vs_destroy(&vs);
      igraph_vector_destroy(&res);
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      return NULL;
    }
    if (igraph_betweenness_estimate(&self->g, &res, vs, PyObject_IsTrue(directed),
        (igraph_integer_t)PyInt_AsLong(cutoff_num), weights,
        PyObject_IsTrue(nobigint))) {
      igraph_vs_destroy(&vs);
      igraph_vector_destroy(&res);
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      Py_DECREF(cutoff_num);
      igraphmodule_handle_igraph_error();
      return NULL;
    }
    Py_DECREF(cutoff_num);
  } else {
    PyErr_SetString(PyExc_TypeError, "cutoff value must be None or integer");
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&res);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
  else
    list = PyFloat_FromDouble(VECTOR(res)[0]);

  igraph_vector_destroy(&res);
  igraph_vs_destroy(&vs);
  if (weights) { igraph_vector_destroy(weights); free(weights); }

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the bibliographic coupling of some vertices in a graph.
 * \return the bibliographic coupling values in a matrix
 * \sa igraph_bibcoupling
 */
PyObject *igraphmodule_Graph_bibcoupling(igraphmodule_GraphObject * self,
                                         PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "vertices", NULL };
  PyObject *vobj = NULL, *list;
  igraph_matrix_t res;
  igraph_vs_t vs;
  igraph_bool_t return_single = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &vobj))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_matrix_init(&res, 1, igraph_vcount(&self->g))) {
    igraph_vs_destroy(&vs);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_bibcoupling(&self->g, &res, vs)) {
    igraph_vs_destroy(&vs);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  /* TODO: Return a single list instead of a matrix if only one vertex was given */
  list = igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);

  igraph_matrix_destroy(&res);
  igraph_vs_destroy(&vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the biconnected components of a graph.
 * \return the list of spanning trees of biconnected components in a PyObject
 * \sa igraph_biconnected_components
 */
PyObject *igraphmodule_Graph_biconnected_components(igraphmodule_GraphObject *self,
	PyObject *args, PyObject *kwds) {
  igraph_vector_ptr_t components;
  igraph_vector_t points;
  igraph_bool_t return_articulation_points;
  igraph_integer_t no;
  PyObject *result, *aps=Py_False;

  static char* kwlist[] = {"return_articulation_points", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &aps)) return NULL;
  return_articulation_points = PyObject_IsTrue(aps);

  if (igraph_vector_ptr_init(&components, 0)) {
	igraphmodule_handle_igraph_error();
	return NULL;
  }
  if (return_articulation_points) {
	if (igraph_vector_init(&points, 0)) {
	  igraphmodule_handle_igraph_error();
	  igraph_vector_ptr_destroy(&components);
	  return NULL;
	}
  }

  if (igraph_biconnected_components(&self->g, &no, &components, 0, 0, return_articulation_points ? &points : 0)) {
    igraphmodule_handle_igraph_error();
	igraph_vector_ptr_destroy(&components);
    if (return_articulation_points) igraph_vector_destroy(&points);
    return NULL;
  }

  result = igraphmodule_vector_ptr_t_to_PyList(&components, IGRAPHMODULE_TYPE_INT);
  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&components, igraph_vector_destroy);
  igraph_vector_ptr_destroy_all(&components);

  if (return_articulation_points) {
	PyObject *result2;
	igraph_vector_sort(&points);
	result2 = igraphmodule_vector_t_to_PyList(&points, IGRAPHMODULE_TYPE_INT);
    igraph_vector_destroy(&points);
	return Py_BuildValue("NN", result, result2); /* references stolen */
  }
  
  return result;
}

/** \ingroup python_interface_graph
 * \brief Returns the one-mode projections of a bipartite graph
 * \return the two projections as new igraph objects
 * \sa igraph_bipartite_projection
 */
PyObject *igraphmodule_Graph_bipartite_projection(igraphmodule_GraphObject * self,
		PyObject* args, PyObject* kwds) {
  PyObject *types_o = Py_None, *multiplicity_o = Py_True, *mul1 = 0, *mul2 = 0;
  igraphmodule_GraphObject *result1, *result2;
  igraph_vector_bool_t* types = 0;
  igraph_vector_t multiplicities[2];
  igraph_t g1, g2;
  long int probe1 = -1;

  static char* kwlist[] = {"types", "multiplicity", "probe1", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|Ol", kwlist, &types_o, &multiplicity_o, &probe1))
    return NULL;

  if (igraphmodule_attrib_to_vector_bool_t(types_o, self, &types, ATTRIBUTE_TYPE_VERTEX))
	return NULL;

  if (PyObject_IsTrue(multiplicity_o)) {
    if (igraph_vector_init(&multiplicities[0], 0)) {
      if (types) { igraph_vector_bool_destroy(types); free(types); }
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (igraph_vector_init(&multiplicities[1], 0)) {
      igraph_vector_destroy(&multiplicities[0]);
      if (types) { igraph_vector_bool_destroy(types); free(types); }
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (igraph_bipartite_projection(&self->g, types, &g1, &g2,
        &multiplicities[0], &multiplicities[1], (igraph_integer_t) probe1)) {
      igraph_vector_destroy(&multiplicities[0]);
      igraph_vector_destroy(&multiplicities[1]);
      if (types) { igraph_vector_bool_destroy(types); free(types); }
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else {
    if (igraph_bipartite_projection(&self->g, types, &g1, &g2, 0, 0,
          (igraph_integer_t) probe1)) {
      if (types) { igraph_vector_bool_destroy(types); free(types); }
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  if (types) { igraph_vector_bool_destroy(types); free(types); }

  CREATE_GRAPH(result1, g1);
  CREATE_GRAPH(result2, g2);

  if (PyObject_IsTrue(multiplicity_o)) {
    mul1 = igraphmodule_vector_t_to_PyList(&multiplicities[0], IGRAPHMODULE_TYPE_INT);
    igraph_vector_destroy(&multiplicities[0]);
    if (mul1 == NULL) {
      igraph_vector_destroy(&multiplicities[1]);
      return NULL;
    }
    mul2 = igraphmodule_vector_t_to_PyList(&multiplicities[1], IGRAPHMODULE_TYPE_INT);
    igraph_vector_destroy(&multiplicities[1]);
    if (mul2 == NULL)
      return NULL;
    return Py_BuildValue("NNNN", result1, result2, mul1, mul2);
  } else {
    return Py_BuildValue("NN", result1, result2);
  }
}

/** \ingroup python_interface_graph
 * \brief Returns the sizes of the two one-mode projections of a bipartite graph
 * \return the two one-mode projections as new igraph objects
 * \sa igraph_bipartite_projection_size
 */
PyObject *igraphmodule_Graph_bipartite_projection_size(igraphmodule_GraphObject * self,
		PyObject* args, PyObject* kwds) {
  PyObject *types_o = Py_None;
  igraph_vector_bool_t* types = 0;
  igraph_integer_t vcount1, vcount2, ecount1, ecount2;

  static char* kwlist[] = {"types", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &types_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_bool_t(types_o, self, &types, ATTRIBUTE_TYPE_VERTEX))
	return NULL;

  if (igraph_bipartite_projection_size(&self->g, types,
			  &vcount1, &ecount1, &vcount2, &ecount2)) {
    if (types) { igraph_vector_bool_destroy(types); free(types); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (types) { igraph_vector_bool_destroy(types); free(types); }

  return Py_BuildValue("llll", (long)vcount1, (long)ecount1, (long)vcount2, (long)ecount2);
}

/** \ingroup python_interface_graph
 * \brief Calculates the closeness centrality of some vertices in a graph.
 * \return the closeness centralities as a list (or a single float)
 * \sa igraph_betweenness
 */
PyObject *igraphmodule_Graph_closeness(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "vertices", "mode", "cutoff", "weights", NULL };
  PyObject *vobj = Py_None, *list = NULL, *cutoff = Py_None,
           *mode_o = Py_None, *weights_o = Py_None;
  igraph_vector_t res, *weights = 0;
  igraph_neimode_t mode = IGRAPH_ALL;
  int return_single = 0;
  igraph_vs_t vs;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO", kwlist, &vobj,
      &mode_o, &cutoff, &weights_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;
  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&res, 0)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) {
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&res);
    return NULL;
  }

  if (cutoff == Py_None) {
    if (igraph_closeness(&self->g, &res, vs, mode, weights)) {
      igraph_vs_destroy(&vs);
      igraph_vector_destroy(&res);
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else if (PyNumber_Check(cutoff)) {
    PyObject *cutoff_num = PyNumber_Int(cutoff);
    if (cutoff_num == NULL) {
      igraph_vs_destroy(&vs); igraph_vector_destroy(&res);
      return NULL;
    }
    if (igraph_closeness_estimate(&self->g, &res, vs, mode,
        (igraph_integer_t)PyInt_AsLong(cutoff_num), weights)) {
      igraph_vs_destroy(&vs);
      igraph_vector_destroy(&res);
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      igraphmodule_handle_igraph_error();
      Py_DECREF(cutoff_num);
      return NULL;
    }
    Py_DECREF(cutoff_num);
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
  else
    list = PyFloat_FromDouble(VECTOR(res)[0]);

  igraph_vector_destroy(&res);
  igraph_vs_destroy(&vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the (weakly or strongly) connected components in a graph.
 * \return a list containing the cluster ID for every vertex in the graph
 * \sa igraph_clusters
 */
PyObject *igraphmodule_Graph_clusters(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "mode", NULL };
  igraph_connectedness_t mode = IGRAPH_STRONG;
  igraph_vector_t res1, res2;
  igraph_integer_t no;
  PyObject *list, *mode_o = Py_None;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_connectedness_t(mode_o, &mode)) return NULL;

  if (mode != IGRAPH_STRONG && mode != IGRAPH_WEAK) {
    PyErr_SetString(PyExc_ValueError, "mode must be either STRONG or WEAK");
    return NULL;
  }

  igraph_vector_init(&res1, igraph_vcount(&self->g));
  igraph_vector_init(&res2, 10);

  if (igraph_clusters(&self->g, &res1, &res2, &no, mode)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&res1);
    igraph_vector_destroy(&res2);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&res1, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&res1);
  igraph_vector_destroy(&res2);
  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates Burt's constraint scores for a given graph
 * \sa igraph_constraint
 */
PyObject *igraphmodule_Graph_constraint(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "vertices", "weights", NULL };
  PyObject *vids_obj = Py_None, *weight_obj = Py_None, *list;
  igraph_vector_t result, weights;
  igraph_vs_t vids;
  igraph_bool_t return_single = 0;

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "|OO", kwlist, &vids_obj, &weight_obj))
    return NULL;

  if (igraph_vector_init(&result, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_PyObject_to_attribute_values(weight_obj, &weights,
                                                self, ATTRHASH_IDX_EDGE,
                                                1.0)) {
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (igraphmodule_PyObject_to_vs_t(vids_obj, &vids, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    igraph_vector_destroy(&weights);
    return NULL;
  }

  if (igraph_constraint(&self->g, &result, vids, &weights)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vids);
    igraph_vector_destroy(&result);
    igraph_vector_destroy(&weights);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vs_destroy(&vids);
  igraph_vector_destroy(&result);
  igraph_vector_destroy(&weights);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the cocitation scores of some vertices in a graph.
 * \return the cocitation scores in a matrix
 * \sa igraph_cocitation
 */
PyObject *igraphmodule_Graph_cocitation(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "vertices", NULL };
  PyObject *vobj = NULL, *list = NULL;
  igraph_matrix_t res;
  int return_single = 0;
  igraph_vs_t vs;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &vobj))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_matrix_init(&res, 1, igraph_vcount(&self->g))) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_cocitation(&self->g, &res, vs)) {
    igraph_matrix_destroy(&res);
    igraph_vs_destroy(&vs);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  /* TODO: Return a single list instead of a matrix if only one vertex was given */
  list = igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);

  igraph_matrix_destroy(&res);
  igraph_vs_destroy(&vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Replaces multiple vertices with a single one.
 * \return None.
 * \sa igraph_contract_vertices
 */
PyObject *igraphmodule_Graph_contract_vertices(igraphmodule_GraphObject * self,
	                                           PyObject * args, PyObject * kwds) {
  static char* kwlist[] = {"mapping", "combine_attrs", NULL };
  PyObject *mapping_o, *combination_o = Py_None;
  igraph_vector_t mapping;
  igraph_attribute_combination_t combination;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &mapping_o,
		&combination_o))
	return NULL;

  if (igraphmodule_PyObject_to_attribute_combination_t(
		combination_o, &combination))
	return NULL;

  if (igraphmodule_PyObject_to_vector_t(mapping_o, &mapping, 1, 0)) {
    igraph_attribute_combination_destroy(&combination);
    return NULL;
  }

  if (igraph_contract_vertices(&self->g, &mapping, &combination)) {
    igraph_attribute_combination_destroy(&combination);
    igraph_vector_destroy(&mapping);
    return NULL;
  }

  igraph_attribute_combination_destroy(&combination);
  igraph_vector_destroy(&mapping);

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Decomposes a graph into components.
 * \return a list of graph objects, each containing a copy of a component in the original graph.
 * \sa igraph_components
 */
PyObject *igraphmodule_Graph_decompose(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "mode", "maxcompno", "minelements", NULL };
  igraph_connectedness_t mode = IGRAPH_STRONG;
  PyObject *list;
  igraphmodule_GraphObject *o;
  long maxcompno = -1, minelements = -1, n, i;
  igraph_vector_ptr_t components;
  igraph_t *g;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lll", kwlist, &mode,
                                   &maxcompno, &minelements))
    return NULL;

  if (mode != IGRAPH_STRONG && mode != IGRAPH_WEAK) {
    PyErr_SetString(PyExc_ValueError, "mode must be either STRONG or WEAK");
    return NULL;
  }

  igraph_vector_ptr_init(&components, 3);
  if (igraph_decompose(&self->g, &components, mode, maxcompno, minelements)) {
    igraph_vector_ptr_destroy(&components);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  /* We have to create a Python igraph object for every graph returned */
  n = igraph_vector_ptr_size(&components);
  list = PyList_New(n);
  for (i = 0; i < n; i++) {
    g = (igraph_t *) VECTOR(components)[i];
    CREATE_GRAPH(o, *g);
    PyList_SET_ITEM(list, i, (PyObject *) o);
    /* reference has been transferred by PyList_SET_ITEM, no need to DECREF
     *
     * we mustn't call igraph_destroy here, because it would free the vertices
     * and the edges as well, but we need them in o->g. So just call free */
    free(g);
  }

  igraph_vector_ptr_destroy(&components);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the eccentricities of some vertices in a graph.
 * \return the eccentricities as a list (or a single float)
 * \sa igraph_eccentricity
 */
PyObject *igraphmodule_Graph_eccentricity(igraphmodule_GraphObject* self,
                                          PyObject* args, PyObject* kwds) {
  static char *kwlist[] = { "vertices", "mode", NULL };
  PyObject *vobj = Py_None, *list = NULL, *mode_o = Py_None;
  igraph_vector_t res;
  igraph_neimode_t mode = IGRAPH_OUT;
  int return_single = 0;
  igraph_vs_t vs;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &vobj, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&res, 0)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_eccentricity(&self->g, &res, vs, mode)) {
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&res);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
  else
    list = PyFloat_FromDouble(VECTOR(res)[0]);

  igraph_vector_destroy(&res);
  igraph_vs_destroy(&vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the edge betweennesses in the graph
 * \return a list containing the edge betweenness for every edge
 * \sa igraph_edge_betweenness
 */
PyObject *igraphmodule_Graph_edge_betweenness(igraphmodule_GraphObject * self,
                                              PyObject * args,
                                              PyObject * kwds)
{
  static char *kwlist[] = { "directed", "cutoff", "weights", NULL };
  igraph_vector_t res, *weights = 0;
  PyObject *list, *directed = Py_True, *cutoff = Py_None;
  PyObject *weights_o = Py_None;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist,
                                   &directed, &cutoff, &weights_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  igraph_vector_init(&res, igraph_ecount(&self->g));

  if (cutoff == Py_None) {
    if (igraph_edge_betweenness(&self->g, &res, PyObject_IsTrue(directed), weights)) {
      igraphmodule_handle_igraph_error();
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      igraph_vector_destroy(&res);
      return NULL;
    }
  } else if (PyNumber_Check(cutoff)) {
    PyObject *cutoff_num = PyNumber_Int(cutoff);
    if (!cutoff_num) {
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      igraph_vector_destroy(&res); return NULL;
    }
    if (igraph_edge_betweenness_estimate(&self->g, &res, PyObject_IsTrue(directed),
        (igraph_integer_t)PyInt_AsLong(cutoff_num), weights)) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(&res);
      if (weights) { igraph_vector_destroy(weights); free(weights); }
      Py_DECREF(cutoff_num);
      return NULL;
    }
    Py_DECREF(cutoff_num);
  } else {
    PyErr_SetString(PyExc_TypeError, "cutoff value must be None or integer");
    igraph_vector_destroy(&res);
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&res);
  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the edge connectivity of the graph
 * \return the edge connectivity
 * \sa igraph_edge_connectivity, igraph_st_edge_connectivity
 */
PyObject *igraphmodule_Graph_edge_connectivity(igraphmodule_GraphObject *self,
        PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "source", "target", "checks", NULL };
  PyObject *checks = Py_True;
  long int source = -1, target = -1, result;
  igraph_integer_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|llO", kwlist,
      &source, &target, &checks))
    return NULL;

  if (source < 0 && target < 0) {
    if (igraph_edge_connectivity(&self->g, &res, PyObject_IsTrue(checks))) {
      igraphmodule_handle_igraph_error();
	  return NULL;
    }
  } else if (source >= 0 && target >= 0) {
    if (igraph_st_edge_connectivity(&self->g, &res, (igraph_integer_t) source,
          (igraph_integer_t) target)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else {
	PyErr_SetString(PyExc_ValueError, "if source or target is given, the other one must also be specified");
	return NULL;
  }

  result = res;

  return Py_BuildValue("l", result);
}

/** \ingroup python_interface_graph
 * \brief Calculates the eigenvector centralities of the vertices in the graph
 * \sa igraph_eigenvector_centrality
 */
PyObject *igraphmodule_Graph_eigenvector_centrality(
  igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] =
    { "directed", "scale", "weights", "arpack_options", "return_eigenvalue",
      NULL };
  PyObject *directed_o = Py_True;
  PyObject *scale_o = Py_True;
  PyObject *weights_o = Py_None;
  PyObject *arpack_options_o = igraphmodule_arpack_options_default;
  igraphmodule_ARPACKOptionsObject *arpack_options;
  PyObject *return_eigenvalue = Py_False;
  PyObject *res_o;
  igraph_real_t value;
  igraph_vector_t *weights=0, res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO!O", kwlist,
                                   &directed_o, &scale_o, &weights_o,
                                   &igraphmodule_ARPACKOptionsType,
                                   &arpack_options, &return_eigenvalue))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  if (igraph_vector_init(&res, 0)) {
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return igraphmodule_handle_igraph_error();
  }

  arpack_options = (igraphmodule_ARPACKOptionsObject*)arpack_options_o;
  if (igraph_eigenvector_centrality(&self->g, &res, &value,
      PyObject_IsTrue(directed_o), PyObject_IsTrue(scale_o),
      weights, igraphmodule_ARPACKOptions_get(arpack_options))) {
    igraphmodule_handle_igraph_error();
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraph_vector_destroy(&res);
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }
  
  res_o = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT); 
  igraph_vector_destroy(&res);
  if (res_o == NULL) return igraphmodule_handle_igraph_error();

  if (PyObject_IsTrue(return_eigenvalue)) {
    PyObject *ev_o = PyFloat_FromDouble((double)value);
    if (ev_o == NULL) {
      Py_DECREF(res_o);
      return igraphmodule_handle_igraph_error();
    }
    return Py_BuildValue("NN", res_o, ev_o);
  }

  return res_o;
}

/** \ingroup python_interface_graph
 * \brief Calculates a feedback arc set for a graph
 * \return a list containing the indices in the chosen feedback arc set
 * \sa igraph_feedback_arc_set
 */
PyObject *igraphmodule_Graph_feedback_arc_set(
    igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "weights", "method", NULL };
  igraph_vector_t* weights = 0;
  igraph_vector_t result;
  igraph_fas_algorithm_t algo = IGRAPH_FAS_APPROX_EADES;
  PyObject *weights_o = Py_None, *result_o = NULL, *algo_o = NULL;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &weights_o, &algo_o))
    return NULL;

  if (igraphmodule_PyObject_to_fas_algorithm_t(algo_o, &algo))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE))
    return NULL;

  if (igraph_vector_init(&result, 0)) {
    if (weights) { igraph_vector_destroy(weights); free(weights); }
  }

  if (igraph_feedback_arc_set(&self->g, &result, weights, algo)) {
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  result_o = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  return result_o;
}


/** \ingroup python_interface_graph
 * \brief Calculates the shortest paths from/to a given node in the graph
 * \return a list containing shortest paths from/to the given node
 * \sa igraph_get_shortest_paths
 */
PyObject *igraphmodule_Graph_get_shortest_paths(igraphmodule_GraphObject *
                                                self, PyObject * args,
                                                PyObject * kwds)
{
  static char *kwlist[] = { "v", "to", "weights", "mode", "output", NULL };
  igraph_vector_t *res, *weights=0;
  igraph_neimode_t mode = IGRAPH_OUT;
  long int i, j;
  igraph_integer_t from, no_of_target_nodes;
  igraph_vs_t to;
  PyObject *list, *item, *mode_o=Py_None, *weights_o=Py_None,
           *output_o=Py_None, *from_o = Py_None, *to_o=Py_None;
  igraph_vector_ptr_t *ptrvec=0;
  igraph_bool_t use_edges = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OOOO!", kwlist, &from_o,
        &to_o, &weights_o, &mode_o, &PyString_Type, &output_o))
    return NULL;

  if (output_o == 0 || output_o == Py_None ||
      PyString_IsEqualToASCIIString(output_o, "vpath")) {
    use_edges = 0;
  } else if (PyString_IsEqualToASCIIString(output_o, "epath")) {
    use_edges = 1;
  } else {
    PyErr_SetString(PyExc_ValueError, "output argument must be \"vpath\" or \"epath\"");
    return NULL;
  }

  if (igraphmodule_PyObject_to_vid(from_o, &from, &self->g))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;
  
  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
      ATTRIBUTE_TYPE_EDGE)) return NULL;

  if (igraphmodule_PyObject_to_vs_t(to_o, &to, &self->g, 0, 0)) return NULL;

  if (igraph_vs_size(&self->g, &to, &no_of_target_nodes)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  ptrvec = (igraph_vector_ptr_t *) calloc(1, sizeof(igraph_vector_ptr_t));
  if (!ptrvec) {
    PyErr_SetString(PyExc_MemoryError, "");
	if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  if (igraph_vector_ptr_init(ptrvec, no_of_target_nodes)) {
    PyErr_SetString(PyExc_MemoryError, "");
    free(ptrvec);
	if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  res = (igraph_vector_t *) calloc(no_of_target_nodes, sizeof(igraph_vector_t));
  if (!res) {
    PyErr_SetString(PyExc_MemoryError, "");
    igraph_vector_ptr_destroy(ptrvec); free(ptrvec);
	if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  for (i = 0; i < no_of_target_nodes; i++) {
    VECTOR(*ptrvec)[i] = &res[i];
    igraph_vector_init(&res[i], 0);
  }

  if (igraph_get_shortest_paths_dijkstra(&self->g, use_edges ? 0 : ptrvec,
        use_edges ? ptrvec : 0, from, to, weights, mode)) {
    igraphmodule_handle_igraph_error();
    for (j = 0; j < no_of_target_nodes; j++) igraph_vector_destroy(&res[j]);
    free(res);
	if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraph_vector_ptr_destroy(ptrvec);
    free(ptrvec);
    return NULL;
  }

  igraph_vector_ptr_destroy(ptrvec);
  free(ptrvec);
  if (weights) { igraph_vector_destroy(weights); free(weights); }

  list = PyList_New(no_of_target_nodes);
  if (!list) {
    for (j = 0; j < no_of_target_nodes; j++) igraph_vector_destroy(&res[j]);
    free(res);
    return NULL;
  }

  for (i = 0; i < no_of_target_nodes; i++) {
    item = igraphmodule_vector_t_to_PyList(&res[i], IGRAPHMODULE_TYPE_INT);
    if (!item || PyList_SetItem(list, i, item)) {
      if (item) {
        Py_DECREF(item);
      }
      Py_DECREF(list);
      for (j = 0; j < no_of_target_nodes; j++) igraph_vector_destroy(&res[j]);
      free(res);
      return NULL;
    }
  }

  for (j = 0; j < no_of_target_nodes; j++) igraph_vector_destroy(&res[j]);
  free(res);
  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates all of the shortest paths from/to a given node in the graph
 * \return a list containing shortest paths from/to the given node
 * \sa igraph_get_shortest_paths
 */
PyObject *igraphmodule_Graph_get_all_shortest_paths(igraphmodule_GraphObject *
                                                    self, PyObject * args,
                                                    PyObject * kwds)
{
  static char *kwlist[] = { "v", "to", "weights", "mode", NULL };
  igraph_vector_ptr_t res;
  igraph_vector_t *weights = 0;
  igraph_neimode_t mode = IGRAPH_OUT;
  long int i, j;
  igraph_integer_t from;
  igraph_vs_t to;
  PyObject *list, *item, *from_o, *mode_o=Py_None, *to_o=Py_None, *weights_o=Py_None;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OOO", kwlist, &from_o,
        &to_o, &weights_o, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (igraphmodule_PyObject_to_vid(from_o, &from, &self->g))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(to_o, &to, &self->g, 0, 0))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
      ATTRIBUTE_TYPE_EDGE)) {
    igraph_vs_destroy(&to);
    return NULL;
  }

  if (igraph_vector_ptr_init(&res, 1)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&to);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  if (igraph_get_all_shortest_paths_dijkstra(&self->g, &res,
        NULL, from, to, weights, mode)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_ptr_destroy(&res);
    igraph_vs_destroy(&to);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  igraph_vs_destroy(&to);
  if (weights) { igraph_vector_destroy(weights); free(weights); }

  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&res, igraph_vector_destroy);

  j = igraph_vector_ptr_size(&res);
  list = PyList_New(j);
  if (!list) {
    igraph_vector_ptr_destroy_all(&res);
    return NULL;
  }

  for (i = 0; i < j; i++) {
    item =
      igraphmodule_vector_t_to_PyList((igraph_vector_t *)
                                      igraph_vector_ptr_e(&res, i),
                    IGRAPHMODULE_TYPE_INT);
    if (!item) {
      Py_DECREF(list);
      igraph_vector_ptr_destroy_all(&res);
      return NULL;
    }
    if (PyList_SetItem(list, i, item)) {
      Py_DECREF(list);
      Py_DECREF(item);
      igraph_vector_ptr_destroy_all(&res);
      return NULL;
    }
  }

  igraph_vector_ptr_destroy_all(&res);
  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates Kleinberg's hub scores of the vertices in the graph
 * \sa igraph_hub_score
 */
PyObject *igraphmodule_Graph_hub_score(
  igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] =
    { "weights", "scale", "arpack_options", "return_eigenvalue", NULL };
  PyObject *scale_o = Py_True, *weights_o = Py_None;
  PyObject *arpack_options_o = igraphmodule_arpack_options_default;
  igraphmodule_ARPACKOptionsObject *arpack_options;
  PyObject *return_eigenvalue = Py_False;
  PyObject *res_o;
  igraph_real_t value;
  igraph_vector_t res, *weights = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO!O", kwlist, &weights_o,
                                   &scale_o, &igraphmodule_ARPACKOptionsType,
                                   &arpack_options, &return_eigenvalue))
    return NULL;

  if (igraph_vector_init(&res, 0)) return igraphmodule_handle_igraph_error();

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  arpack_options = (igraphmodule_ARPACKOptionsObject*)arpack_options_o;
  if (igraph_hub_score(&self->g, &res, &value, PyObject_IsTrue(scale_o),
      weights, igraphmodule_ARPACKOptions_get(arpack_options))) {
    igraphmodule_handle_igraph_error();
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraph_vector_destroy(&res);
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  res_o = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT); 
  igraph_vector_destroy(&res);
  if (res_o == NULL) return igraphmodule_handle_igraph_error();

  if (PyObject_IsTrue(return_eigenvalue)) {
    PyObject *ev_o = PyFloat_FromDouble((double)value);
    if (ev_o == NULL) {
      Py_DECREF(res_o);
      return igraphmodule_handle_igraph_error();
    }
    return Py_BuildValue("NN", res_o, ev_o);
  }

  return res_o;
}


/** \ingroup python_interface_graph
 * \brief Returns the line graph of the graph
 * \return the line graph as a new igraph object
 * \sa igraph_linegraph
 */
PyObject *igraphmodule_Graph_linegraph(igraphmodule_GraphObject * self) {
  igraph_t lg;
  igraphmodule_GraphObject *result;

  if (igraph_linegraph(&self->g, &lg)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  CREATE_GRAPH(result, lg);

  return (PyObject *) result;
}

/**
 * \ingroup python_interface_graph
 * \brief Returns the k-neighborhood of some vertices in the
 *   graph.
 * \sa igraph_neighborhood
 */
PyObject *igraphmodule_Graph_neighborhood(igraphmodule_GraphObject *self,
    PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "vertices", "order", "mode", NULL };
  PyObject *vobj = Py_None;
  PyObject *mode_o = 0;
  PyObject *result;
  long int order = 1;
  igraph_neimode_t mode = IGRAPH_ALL;
  igraph_bool_t return_single = 0;
  igraph_vs_t vs;
  igraph_vector_ptr_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OlO", kwlist,
        &vobj, &order, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_vector_ptr_init(&res, 0)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_neighborhood(&self->g, &res, vs, (igraph_integer_t) order, mode)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  igraph_vs_destroy(&vs);

  if (!return_single)
    result = igraphmodule_vector_ptr_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  else
    result = igraphmodule_vector_t_to_PyList((igraph_vector_t*)VECTOR(res)[0],
        IGRAPHMODULE_TYPE_INT);

  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&res, igraph_vector_destroy);
  igraph_vector_ptr_destroy_all(&res);

  return result;
}

/**
 * \ingroup python_interface_graph
 * \brief Returns the size of the k-neighborhood of some vertices in the
 *   graph.
 * \sa igraph_neighborhood_size
 */
PyObject *igraphmodule_Graph_neighborhood_size(igraphmodule_GraphObject *self,
    PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "vertices", "order", "mode", NULL };
  PyObject *vobj = Py_None;
  PyObject *mode_o = 0;
  PyObject *result;
  long int order = 1;
  igraph_neimode_t mode = IGRAPH_ALL;
  igraph_bool_t return_single = 0;
  igraph_vs_t vs;
  igraph_vector_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OlO", kwlist,
        &vobj, &order, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_vector_init(&res, 0)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_neighborhood_size(&self->g, &res, vs, (igraph_integer_t) order, mode)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  igraph_vs_destroy(&vs);

  if (!return_single)
    result = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  else
    result = PyInt_FromLong((long)VECTOR(res)[0]);

  igraph_vector_destroy(&res);

  return result;
}

/** \ingroup python_interface_graph
 * \brief Calculates the Google personalized PageRank value of some vertices in the graph.
 * \return the personalized PageRank values
 * \sa igraph_personalized_pagerank
 */
PyObject *igraphmodule_Graph_personalized_pagerank(igraphmodule_GraphObject *self,
                                      PyObject *args, PyObject *kwds)
{
  static char *kwlist[] =
    { "vertices", "directed", "damping", "reset", "reset_vertices", "weights",
      "arpack_options", NULL };
  PyObject *directed = Py_True;
  PyObject *vobj = Py_None, *wobj = Py_None, *robj = Py_None, *rvsobj = Py_None;
  PyObject *list;
  PyObject *arpack_options_o = igraphmodule_arpack_options_default;
  igraphmodule_ARPACKOptionsObject *arpack_options;
  double damping = 0.85;
  igraph_vector_t res;
  igraph_vector_t *reset = 0;
  igraph_vector_t weights;
  igraph_bool_t return_single = 0;
  igraph_vs_t vs, reset_vs;
  int retval;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOdOOOO!", kwlist, &vobj,
                                   &directed, &damping, &robj, &rvsobj, &wobj,
                                   &igraphmodule_ARPACKOptionsType,
                                   &arpack_options_o))
    return NULL;

  if (robj != Py_None && rvsobj != Py_None) {
    PyErr_SetString(PyExc_ValueError, "only reset or reset_vs can be defined, not both");
    return NULL;
  }

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (robj != Py_None) {
    if (igraphmodule_attrib_to_vector_t(robj, self, &reset, ATTRIBUTE_TYPE_VERTEX)) {
      igraph_vs_destroy(&vs);
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else if (rvsobj != Py_None) {
    if (igraphmodule_PyObject_to_vs_t(rvsobj, &reset_vs, &self->g, 0, 0)) {
      igraph_vs_destroy(&vs);
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  if (igraphmodule_PyObject_to_attribute_values(wobj, &weights,
                                                self, ATTRHASH_IDX_EDGE,
                                                1.0)) {
    igraph_vs_destroy(&vs);
    if (rvsobj != Py_None) igraph_vs_destroy(&reset_vs);
    if (reset) { igraph_vector_destroy(reset); free(reset); }
    return NULL;
  }

  if (igraph_vector_init(&res, 0)) {
    igraph_vs_destroy(&vs);
    if (rvsobj != Py_None) igraph_vs_destroy(&reset_vs);
    if (reset) { igraph_vector_destroy(reset); free(reset); }
    igraph_vector_destroy(&weights);
    return igraphmodule_handle_igraph_error();
  }

  arpack_options = (igraphmodule_ARPACKOptionsObject*)arpack_options_o;
  if (rvsobj != Py_None)
    retval = igraph_personalized_pagerank_vs(&self->g, &res, 0, vs, PyObject_IsTrue(directed),
        damping, reset_vs, &weights, igraphmodule_ARPACKOptions_get(arpack_options));
  else
    retval = igraph_personalized_pagerank(&self->g, &res, 0, vs, PyObject_IsTrue(directed),
        damping, reset, &weights, igraphmodule_ARPACKOptions_get(arpack_options));

  if (retval) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    if (rvsobj != Py_None) igraph_vs_destroy(&reset_vs);
    if (reset) { igraph_vector_destroy(reset); free(reset); }
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&res);
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
  else
    list = PyFloat_FromDouble(VECTOR(res)[0]);

  igraph_vector_destroy(&res);
  igraph_vs_destroy(&vs);
  if (rvsobj != Py_None) igraph_vs_destroy(&reset_vs);
  igraph_vector_destroy(&weights);
  if (reset) { igraph_vector_destroy(reset); free(reset); }

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the path length histogram of the graph
 * \sa igraph_path_length_hist
 */
PyObject *igraphmodule_Graph_path_length_hist(igraphmodule_GraphObject *self,
                                              PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "directed", NULL };
  PyObject *directed = Py_True, *result;
  igraph_real_t unconn;
  igraph_vector_t res;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &directed))
	return NULL;

  if (igraph_vector_init(&res, 0))
	return igraphmodule_handle_igraph_error();

  if (igraph_path_length_hist(&self->g, &res, &unconn, PyObject_IsTrue(directed))) {
	igraph_vector_destroy(&res);
	return igraphmodule_handle_igraph_error();
  }
  
  result=igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&res);
  return Py_BuildValue("Nd", result, (double)unconn);
}

/** \ingroup python_interface_graph
 * \brief Permutes the vertices of the graph 
 * \return the new graph as a new igraph object
 * \sa igraph_permute_vertices
 */
PyObject *igraphmodule_Graph_permute_vertices(igraphmodule_GraphObject *self,
                                              PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "permutation", NULL };
  igraph_t pg;
  igraph_vector_t perm;
  igraphmodule_GraphObject *result;
  PyObject *list;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, &PyList_Type, &list))
    return NULL;

  if (igraphmodule_PyObject_to_vector_t(list, &perm, 1, 0)) return NULL;

  if (igraph_permute_vertices(&self->g, &pg, &perm)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&perm);
    return NULL;
  }

  igraph_vector_destroy(&perm);

  CREATE_GRAPH(result, pg);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Rewires a graph while preserving degree distribution
 * \return the rewired graph
 * \sa igraph_rewire
 */
PyObject *igraphmodule_Graph_rewire(igraphmodule_GraphObject * self,
                                    PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "n", "mode", NULL };
  long int n = 1000;
  PyObject *mode_o = Py_None;
  igraph_rewiring_t mode = IGRAPH_REWIRING_SIMPLE;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lO", kwlist, &n, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_rewiring_t(mode_o, &mode))
    return NULL;

  if (igraph_rewire(&self->g, (igraph_integer_t) n, mode)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Rewires the edges of a graph wth constant probability
 * \return the rewired graph
 * \sa igraph_rewire_edges
 */
PyObject *igraphmodule_Graph_rewire_edges(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "prob", "loops", "multiple", NULL };
  double prob;
  PyObject *loops_o = Py_False, *multiple_o = Py_False;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|OO", kwlist,
        &prob, &loops_o, &multiple_o))
    return NULL;

  if (igraph_rewire_edges(&self->g, prob, PyObject_IsTrue(loops_o),
        PyObject_IsTrue(multiple_o))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Calculates shortest paths in a graph.
 * \return the shortest path lengths for the given vertices
 * \sa igraph_shortest_paths, igraph_shortest_paths_dijkstra,
 *     igraph_shortest_paths_bellman_ford, igraph_shortest_paths_johnson
 */
PyObject *igraphmodule_Graph_shortest_paths(igraphmodule_GraphObject * self,
                                            PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "source", "target", "weights", "mode", NULL };
  PyObject *from_o = NULL, *to_o = NULL, *mode_o = NULL, *weights_o = Py_None;
  PyObject *list = NULL;
  igraph_matrix_t res;
  igraph_vector_t *weights=0;
  igraph_neimode_t mode = IGRAPH_OUT;
  int return_single_from = 0, return_single_to = 0, e = 0;
  igraph_vs_t from_vs, to_vs;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO", kwlist,
        &from_o, &to_o, &weights_o, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return 0;
  if (igraphmodule_PyObject_to_vs_t(from_o, &from_vs, &self->g, &return_single_from, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_PyObject_to_vs_t(to_o, &to_vs, &self->g, &return_single_to, 0)) {
    igraph_vs_destroy(&from_vs);
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
      ATTRIBUTE_TYPE_EDGE)) {
    igraph_vs_destroy(&from_vs);
    igraph_vs_destroy(&to_vs);
    return NULL;
  }

  if (igraph_matrix_init(&res, 1, igraph_vcount(&self->g))) {
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraph_vs_destroy(&from_vs);
    igraph_vs_destroy(&to_vs);
    return igraphmodule_handle_igraph_error();
  }

  /* Select the most suitable algorithm */
  if (weights) {
    if (igraph_vector_min(weights) > 0) {
      /* Only positive weights, use Dijkstra's algorithm */
      e = igraph_shortest_paths_dijkstra(&self->g, &res, from_vs, to_vs, weights, mode);
    } else {
      /* There are negative weights. For a small number of sources, use Bellman-Ford.
       * Otherwise, use Johnson's algorithm */
      igraph_integer_t vs_size;
      e = igraph_vs_size(&self->g, &from_vs, &vs_size);
      if (!e) {
        if (vs_size <= 100 || mode != IGRAPH_OUT) {
          e = igraph_shortest_paths_bellman_ford(&self->g, &res, from_vs, to_vs, weights, mode);
        } else {
          e = igraph_shortest_paths_johnson(&self->g, &res, from_vs, to_vs, weights);
        }
      }
    }
  } else {
    /* No weights, use a simple BFS */
    e = igraph_shortest_paths(&self->g, &res, from_vs, to_vs, mode);
  }

  if (e) {
    if (weights) igraph_vector_destroy(weights);
    igraph_matrix_destroy(&res);
    igraph_vs_destroy(&from_vs);
    igraph_vs_destroy(&to_vs);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (weights) {
    igraph_vector_destroy(weights);
    list = igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
  } else {
    list = igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }

  igraph_matrix_destroy(&res);
  igraph_vs_destroy(&from_vs);
  igraph_vs_destroy(&to_vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the Jaccard similarities of some vertices in a graph.
 * \return the similarity scores in a matrix
 * \sa igraph_similarity_jaccard
 */
PyObject *igraphmodule_Graph_similarity_jaccard(igraphmodule_GraphObject * self,
  PyObject * args, PyObject * kwds) {
  static char *kwlist[] = { "vertices", "pairs", "mode", "loops", NULL };
  PyObject *vertices_o = Py_None, *pairs_o = Py_None;
  PyObject *list = NULL, *loops = Py_True, *mode_o = Py_None;
  igraph_neimode_t mode = IGRAPH_ALL;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO", kwlist, &vertices_o,
        &pairs_o, &mode_o, &loops))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (vertices_o != Py_None && pairs_o != Py_None) {
    PyErr_SetString(PyExc_ValueError, "at most one of `vertices` and `pairs` "
        "must be given");
    return NULL;
  }

  if (pairs_o == Py_None) {
    /* Case #1: vertices, returning matrix */
    igraph_matrix_t res;
    igraph_vs_t vs;
    int return_single = 0;

    if (igraphmodule_PyObject_to_vs_t(vertices_o, &vs, &self->g, &return_single, 0))
      return NULL;

    if (igraph_matrix_init(&res, 0, 0)) {
      igraph_vs_destroy(&vs);
      return igraphmodule_handle_igraph_error();
    }

    if (igraph_similarity_jaccard(&self->g, &res, vs, mode, PyObject_IsTrue(loops))) {
      igraph_matrix_destroy(&res);
      igraph_vs_destroy(&vs);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    igraph_vs_destroy(&vs);

    list = igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
    igraph_matrix_destroy(&res);
  } else {
    /* Case #2: vertex pairs or edges, returning list */
    igraph_vector_t edges;
    igraph_vector_t res;

    if (igraphmodule_PyObject_to_vector_t(pairs_o, &edges, 1, 1))
      return NULL;

    if (igraph_vector_init(&res, igraph_vector_size(&edges) / 2)) {
      igraph_vector_destroy(&edges);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (igraph_similarity_jaccard_pairs(&self->g, &res, &edges, mode,
          PyObject_IsTrue(loops))) {
      igraph_vector_destroy(&res);
      igraph_vector_destroy(&edges);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    igraph_vector_destroy(&edges);

    list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
    igraph_vector_destroy(&res);
  }

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the Dice similarities of some vertices in a graph.
 * \return the similarity scores in a matrix
 * \sa igraph_similarity_dice
 */
PyObject *igraphmodule_Graph_similarity_dice(igraphmodule_GraphObject * self,
  PyObject * args, PyObject * kwds) {
  static char *kwlist[] = { "vertices", "pairs", "mode", "loops", NULL };
  PyObject *vertices_o = Py_None, *pairs_o = Py_None;
  PyObject *list = NULL, *loops = Py_True, *mode_o = Py_None;
  igraph_neimode_t mode = IGRAPH_ALL;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOO", kwlist, &vertices_o,
        &pairs_o, &mode_o, &loops))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (vertices_o != Py_None && pairs_o != Py_None) {
    PyErr_SetString(PyExc_ValueError, "at most one of `vertices` and `pairs` "
        "must be given");
    return NULL;
  }

  if (pairs_o == Py_None) {
    /* Case #1: vertices, returning matrix */
    igraph_matrix_t res;
    igraph_vs_t vs;
    int return_single = 0;

    if (igraphmodule_PyObject_to_vs_t(vertices_o, &vs, &self->g, &return_single, 0))
      return NULL;

    if (igraph_matrix_init(&res, 0, 0)) {
      igraph_vs_destroy(&vs);
      return igraphmodule_handle_igraph_error();
    }

    if (igraph_similarity_dice(&self->g, &res, vs, mode, PyObject_IsTrue(loops))) {
      igraph_matrix_destroy(&res);
      igraph_vs_destroy(&vs);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    igraph_vs_destroy(&vs);

    list = igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
    igraph_matrix_destroy(&res);
  } else {
    /* Case #2: vertex pairs or edges, returning list */
    igraph_vector_t edges;
    igraph_vector_t res;

    if (igraphmodule_PyObject_to_vector_t(pairs_o, &edges, 1, 1))
      return NULL;

    if (igraph_vector_init(&res, igraph_vector_size(&edges) / 2)) {
      igraph_vector_destroy(&edges);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (igraph_similarity_dice_pairs(&self->g, &res, &edges, mode,
          PyObject_IsTrue(loops))) {
      igraph_vector_destroy(&res);
      igraph_vector_destroy(&edges);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    igraph_vector_destroy(&edges);

    list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);
    igraph_vector_destroy(&res);
  }

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the inverse log-weighted similarities of some vertices in
 * a graph.
 * \return the similarity scores in a matrix
 * \sa igraph_similarity_inverse_log_weighted
 */
PyObject *igraphmodule_Graph_similarity_inverse_log_weighted(
  igraphmodule_GraphObject * self, PyObject * args, PyObject * kwds) {
  static char *kwlist[] = { "vertices", "mode", NULL };
  PyObject *vobj = NULL, *list = NULL, *mode_o = Py_None;
  igraph_matrix_t res;
  igraph_neimode_t mode = IGRAPH_ALL;
  int return_single = 0;
  igraph_vs_t vs;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &vobj, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;
  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) return NULL; 

  if (igraph_matrix_init(&res, 0, 0)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_similarity_inverse_log_weighted(&self->g,&res,vs,mode)) {
    igraph_matrix_destroy(&res);
    igraph_vs_destroy(&vs);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  list = igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_FLOAT);

  igraph_matrix_destroy(&res);
  igraph_vs_destroy(&vs);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates a spanning tree for a graph
 * \return the spanning tree or a list of edges participating in the spanning tree
 * \sa igraph_minimum_spanning_tree_unweighted
 * \sa igraph_minimum_spanning_tree_unweighted
 * \sa igraph_minimum_spanning_tree_prim
 */
PyObject *igraphmodule_Graph_spanning_tree(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "weights", NULL };
  igraph_vector_t* ws = 0;
  igraph_vector_t res;
  PyObject *weights_o = Py_None, *result = NULL;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &weights_o))
    return NULL;

  if (igraph_vector_init(&res, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &ws, ATTRIBUTE_TYPE_EDGE)) {
    igraph_vector_destroy(&res);
    return NULL;
  }

  if (igraph_minimum_spanning_tree(&self->g, &res, ws)) {
	if (ws != 0) { igraph_vector_destroy(ws); free(ws); }
    igraph_vector_destroy(&res);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (ws != 0) { igraph_vector_destroy(ws); free(ws); }
  result = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&res);
  return result;
}

/** \ingroup python_interface_graph
 * \brief Simplifies a graph by removing loops and/or multiple edges
 * \return the simplified graph.
 * \sa igraph_simplify
 */
PyObject *igraphmodule_Graph_simplify(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "multiple", "loops", "combine_edges", NULL };
  PyObject *multiple = Py_True, *loops = Py_True, *comb_o = Py_None;
  igraph_attribute_combination_t comb;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist,
                                   &multiple, &loops, &comb_o))
    return NULL;

  if (igraphmodule_PyObject_to_attribute_combination_t(comb_o, &comb))
	return NULL;

  if (igraph_simplify(&self->g, PyObject_IsTrue(multiple),
                      PyObject_IsTrue(loops), &comb)) {
	igraph_attribute_combination_destroy(&comb);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  igraph_attribute_combination_destroy(&comb);
  Py_INCREF(self);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Calculates the vertex indices within the same component as a given vertex
 * \return the vertex indices in a list
 * \sa igraph_subcomponent
 */
PyObject *igraphmodule_Graph_subcomponent(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "v", "mode", NULL };
  igraph_vector_t res;
  igraph_neimode_t mode = IGRAPH_ALL;
  igraph_integer_t from;
  PyObject *list = NULL, *mode_o = Py_None, *from_o = Py_None;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &from_o, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (igraphmodule_PyObject_to_vid(from_o, &from, &self->g))
    return NULL;

  igraph_vector_init(&res, 0);
  if (igraph_subcomponent(&self->g, &res, from, mode)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&res);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&res);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Returns an induced subgraph of the graph based on the given vertices
 * \return the subgraph as a new igraph object
 * \sa igraph_induced_subgraph
 */
PyObject *igraphmodule_Graph_induced_subgraph(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "vertices", "implementation", NULL };
  igraph_vs_t vs;
  igraph_t sg;
  igraphmodule_GraphObject *result;
  PyObject *list, *impl_o = Py_None;
  igraph_subgraph_implementation_t impl = IGRAPH_SUBGRAPH_AUTO;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &list, &impl_o))
    return NULL;

  if (igraphmodule_PyObject_to_subgraph_implementation_t(impl_o, &impl))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, 0, 0))
    return NULL;

  if (igraph_induced_subgraph(&self->g, &sg, vs, impl)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    return NULL;
  }

  igraph_vs_destroy(&vs);

  CREATE_GRAPH(result, sg);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Returns a subgraph of the graph based on the given edges
 * \return the subgraph as a new igraph object
 * \sa igraph_subgraph_edges
 */
PyObject *igraphmodule_Graph_subgraph_edges(igraphmodule_GraphObject * self,
                                            PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "edges", "delete_vertices", NULL };
  igraph_es_t es;
  igraph_t sg;
  igraphmodule_GraphObject *result;
  PyObject *list, *delete_vertices = Py_True;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &list, &delete_vertices))
    return NULL;

  if (igraphmodule_PyObject_to_es_t(list, &es, &self->g, 0))
    return NULL;

  if (igraph_subgraph_edges(&self->g, &sg, es, PyObject_IsTrue(delete_vertices))) {
    igraphmodule_handle_igraph_error();
    igraph_es_destroy(&es);
    return NULL;
  }

  CREATE_GRAPH(result, sg);

  igraph_es_destroy(&es);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Calculates the graph transitivity (a.k.a. clustering coefficient)
 * \return the clustering coefficient
 * \sa igraph_transitivity_undirected
 */
PyObject *igraphmodule_Graph_transitivity_undirected(igraphmodule_GraphObject
                                                     * self, PyObject * args,
                                                     PyObject * kwds)
{
  static char *kwlist[] = { "mode", NULL };
  igraph_real_t res;
  PyObject *r, *mode_o = Py_None;
  igraph_transitivity_mode_t mode = IGRAPH_TRANSITIVITY_NAN;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_transitivity_mode_t(mode_o, &mode))
    return NULL;


  if (igraph_transitivity_undirected(&self->g, &res, mode)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  r = Py_BuildValue("d", (double)(res));
  return r;
}

/** \ingroup python_interface_graph
 * \brief Calculates the average of vertex transitivities over the graph
 * \sa igraph_transitivity_avglocal_undirected
 */
PyObject *igraphmodule_Graph_transitivity_avglocal_undirected(igraphmodule_GraphObject
                                                              * self, PyObject * args,
                                                              PyObject * kwds)
{
  static char *kwlist[] = { "mode", NULL };
  igraph_real_t res;
  PyObject *r, *mode_o = Py_None;
  igraph_transitivity_mode_t mode = IGRAPH_TRANSITIVITY_NAN;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_transitivity_mode_t(mode_o, &mode))
    return NULL;

  if (igraph_transitivity_avglocal_undirected(&self->g, &res, mode)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  r = Py_BuildValue("d", (double)(res));
  return r;
}

/** \ingroup python_interface_graph
 * \brief Calculates the local transitivity of given vertices
 * \return the transitivities in a list
 * \sa igraph_transitivity_local_undirected
 */
PyObject
  *igraphmodule_Graph_transitivity_local_undirected(igraphmodule_GraphObject *
                                                    self, PyObject * args,
                                                    PyObject * kwds)
{
  static char *kwlist[] = { "vertices", "mode", "weights", NULL };
  PyObject *vobj = NULL, *mode_o = Py_None, *list = NULL;
  PyObject *weights_o = Py_None;
  igraph_vector_t result;
  igraph_vector_t *weights = 0;
  igraph_bool_t return_single = 0;
  igraph_vs_t vs;
  igraph_transitivity_mode_t mode = IGRAPH_TRANSITIVITY_NAN;
  int retval;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist, &vobj, &mode_o, &weights_o))
    return NULL;

  if (igraphmodule_PyObject_to_transitivity_mode_t(mode_o, &mode))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(vobj, &vs, &self->g, &return_single, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&result, 0)) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) {
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (weights == 0) {
    retval = igraph_transitivity_local_undirected(&self->g, &result, vs, mode);
  } else {
    retval = igraph_transitivity_barrat(&self->g, &result, vs, weights, mode);
  }

  igraph_vs_destroy(&vs);
  if (weights) {
    igraph_vector_destroy(weights); free(weights);
  }

  if (retval) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  if (!return_single)
    list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_FLOAT);
  else
    list = PyFloat_FromDouble(VECTOR(result)[0]);

  igraph_vector_destroy(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates a possible topological sorting
 * \return a possible topological sorting as a list
 * \sa igraph_topological_sorting
 */
PyObject *igraphmodule_Graph_topological_sorting(igraphmodule_GraphObject *
                                                 self, PyObject * args,
                                                 PyObject * kwds)
{
  static char *kwlist[] = { "mode", "warnings", NULL };
  PyObject *list, *mode_o=Py_None;
  PyObject *warnings_o=Py_True;
  igraph_neimode_t mode = IGRAPH_OUT;
  igraph_vector_t result;
  igraph_warning_handler_t* old_handler = 0;
  int retval;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &mode_o, &warnings_o))
    return NULL;
  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;

  if (igraph_vector_init(&result, 0)) 
    return igraphmodule_handle_igraph_error();

  if (!PyObject_IsTrue(warnings_o)) {
    /* Turn off the warnings temporarily */
    old_handler = igraph_set_warning_handler(igraph_warning_handler_ignore);
  }

  retval = igraph_topological_sorting(&self->g, &result, mode);

  if (!PyObject_IsTrue(warnings_o)) {
    /* Restore the warning handler */
    igraph_set_warning_handler(old_handler);
  }

  if (retval) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Calculates the vertex connectivity of the graph
 * \return the vertex connectivity
 * \sa igraph_vertex_connectivity, igraph_st_vertex_connectivity
 */
PyObject *igraphmodule_Graph_vertex_connectivity(igraphmodule_GraphObject *self,
        PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "source", "target", "checks", "neighbors", NULL };
  PyObject *checks = Py_True, *neis = Py_None;
  long int source = -1, target = -1, result;
  igraph_integer_t res;
  igraph_vconn_nei_t neighbors = IGRAPH_VCONN_NEI_ERROR;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|llOO", kwlist,
      &source, &target, &checks, &neis))
    return NULL;

  if (source < 0 && target < 0) {
    if (igraph_vertex_connectivity(&self->g, &res, PyObject_IsTrue(checks))) {
      igraphmodule_handle_igraph_error();
	  return NULL;
    }
  } else if (source >= 0 && target >= 0) {
    if (igraphmodule_PyObject_to_vconn_nei_t(neis, &neighbors))
      return NULL;
    if (igraph_st_vertex_connectivity(&self->g, &res,
          (igraph_integer_t) source, (igraph_integer_t) target, neighbors)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else {
	PyErr_SetString(PyExc_ValueError, "if source or target is given, the other one must also be specified");
	return NULL;
  }

  if (!IGRAPH_FINITE(res)) return Py_BuildValue("d", (double)res);
  
  result = (long)res;
  return Py_BuildValue("l", result);
}

/**********************************************************************
 * Bipartite graphs                                                   *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Checks whether a graph is bipartite
 * \return a boolean or a (boolean, list of booleans) pair
 * \sa igraph_is_bipartite
 */
PyObject *igraphmodule_Graph_is_bipartite(igraphmodule_GraphObject *self,
                                          PyObject *args, PyObject *kwds) {
  PyObject *types_o, *return_types_o = Py_False;
  igraph_vector_bool_t types;
  igraph_bool_t return_types = 0, result;

  static char *kwlist[] = { "return_types", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &return_types_o))
    return NULL;
  return_types = PyObject_IsTrue(return_types_o);

  if (return_types) {
    if (igraph_vector_bool_init(&types, igraph_vcount(&self->g))) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (igraph_is_bipartite(&self->g, &result, &types)) {
      igraph_vector_bool_destroy(&types);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (result) {
      types_o = igraphmodule_vector_bool_t_to_PyList(&types);
      if (!types_o) {
        igraph_vector_bool_destroy(&types);
        return NULL;
      }
      igraph_vector_bool_destroy(&types);
      // reference to types_o will be stolen by Py_BuildValue
      return Py_BuildValue("ON", Py_True, types_o);
    } else {
      igraph_vector_bool_destroy(&types);
      return Py_BuildValue("OO", Py_False, Py_None);
    }
  } else {
    if (igraph_is_bipartite(&self->g, &result, 0)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    if (result)
      Py_RETURN_TRUE;
    else
      Py_RETURN_FALSE;
  }
}

/**********************************************************************
 * Motifs, dyad and triad census                                      *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Calculates the dyad census of the graph
 * \return the dyad census as a 3-tuple
 * \sa igraph_dyad_census
 */
PyObject *igraphmodule_Graph_dyad_census(igraphmodule_GraphObject *self) {
  igraph_integer_t mut, asym, nul;
  PyObject *list;

  if (igraph_dyad_census(&self->g, &mut, &asym, &nul)) {
    return igraphmodule_handle_igraph_error();
  }

  list = Py_BuildValue("lll", (long)mut, (long)asym, (long)nul);
  return list;
}

typedef struct {
  PyObject* func;
  PyObject* graph;
} igraphmodule_i_Graph_motifs_randesu_callback_data_t;

igraph_bool_t igraphmodule_i_Graph_motifs_randesu_callback(const igraph_t *graph,
    igraph_vector_t *vids, int isoclass, void* extra) {
  igraphmodule_i_Graph_motifs_randesu_callback_data_t* data =
    (igraphmodule_i_Graph_motifs_randesu_callback_data_t*)extra;
  PyObject* vector;
  PyObject* result;
  igraph_bool_t retval;

  vector = igraphmodule_vector_t_to_PyList(vids, IGRAPHMODULE_TYPE_INT);
  if (vector == NULL) {
    /* Error in conversion, return 1 */
    return 1;
  }

  result = PyObject_CallFunction(data->func, "OOi", data->graph, vector, isoclass);
  Py_DECREF(vector);

  if (result == NULL) {
    /* Error in callback, return 1 */
    return 1;
  }

  retval = PyObject_IsTrue(result);
  Py_DECREF(result);

  return retval;
}

/** \ingroup python_interface_graph
 * \brief Counts the motifs of the graph sorted by isomorphism classes 
 * \return the number of motifs found for each isomorphism class
 * \sa igraph_motifs_randesu
 */
PyObject *igraphmodule_Graph_motifs_randesu(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  igraph_vector_t result, cut_prob;
  long int size=3;
  PyObject* cut_prob_list=Py_None;
  PyObject* callback=Py_None;
  PyObject *list;
  static char* kwlist[] = {"size", "cut_prob", "callback", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lOO", kwlist, &size,
        &cut_prob_list, &callback))
    return NULL;

  if (cut_prob_list == Py_None) {
    if (igraph_vector_init(&cut_prob, size)) {
      return igraphmodule_handle_igraph_error();
    }
    igraph_vector_fill(&cut_prob, 0);
  } else {
    if (igraphmodule_PyObject_float_to_vector_t(cut_prob_list, &cut_prob))
      return NULL;
  }

  if (callback == Py_None) {
    if (igraph_vector_init(&result, 1)) {
      igraph_vector_destroy(&cut_prob);
      return igraphmodule_handle_igraph_error();
    }
    if (igraph_motifs_randesu(&self->g, &result, (igraph_integer_t) size, &cut_prob)) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(&result);
      igraph_vector_destroy(&cut_prob);
      return NULL;
    }
    igraph_vector_destroy(&cut_prob);

    list = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
    igraph_vector_destroy(&result);

    return list;
  } else if (PyCallable_Check(callback)) {
    igraphmodule_i_Graph_motifs_randesu_callback_data_t data;
    data.graph = (PyObject*)self;
    data.func = callback;
    if (igraph_motifs_randesu_callback(&self->g, (igraph_integer_t) size, &cut_prob,
          igraphmodule_i_Graph_motifs_randesu_callback, &data)) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(&cut_prob);
      return NULL;
    }
    igraph_vector_destroy(&cut_prob);
    /* Don't let exceptions from the callback function go unnoticed */
    if (PyErr_Occurred())
      return NULL;
    Py_RETURN_NONE;
  } else {
    PyErr_SetString(PyExc_TypeError, "callback must be callable or None");
    return NULL;
  }
}

/** \ingroup python_interface_graph
 * \brief Counts the total number of motifs of the graph
 * \return the total number of motifs
 * \sa igraph_motifs_randesu
 */
PyObject *igraphmodule_Graph_motifs_randesu_no(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  igraph_vector_t cut_prob;
  igraph_integer_t result;
  long int size=3;
  PyObject* cut_prob_list=Py_None;
  static char* kwlist[] = {"size", "cut_prob", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lO", kwlist, &size, &cut_prob_list))
    return NULL;

  if (cut_prob_list == Py_None) {
    if (igraph_vector_init(&cut_prob, size)) {
      return igraphmodule_handle_igraph_error();
    }
    igraph_vector_fill(&cut_prob, 0);
  } else {
    if (igraphmodule_PyObject_float_to_vector_t(cut_prob_list, &cut_prob)) {
      return NULL;
    }
  }
  if (igraph_motifs_randesu_no(&self->g, &result, (igraph_integer_t) size, &cut_prob)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&cut_prob);
    return NULL;
  }
  igraph_vector_destroy(&cut_prob);

  return PyInt_FromLong((long)result);
}

/** \ingroup python_interface_graph
 * \brief Estimates the total number of motifs of the graph
 * \return the estimated total number of motifs
 * \sa igraph_motifs_randesu_estimate
 */
PyObject *igraphmodule_Graph_motifs_randesu_estimate(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  igraph_vector_t cut_prob;
  igraph_integer_t result;
  long size=3;
  PyObject* cut_prob_list=Py_None;
  PyObject *sample=Py_None;
  static char* kwlist[] = {"size", "cut_prob", "sample", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lOO", kwlist,
      &size, &cut_prob_list, &sample))
    return NULL;

  if (sample == Py_None) {
    PyErr_SetString(PyExc_TypeError, "sample size must be given");
    return NULL;
  }

  if (cut_prob_list == Py_None) {
    if (igraph_vector_init(&cut_prob, size)) {
      return igraphmodule_handle_igraph_error();
    }
    igraph_vector_fill(&cut_prob, 0);
  } else {
    if (igraphmodule_PyObject_float_to_vector_t(cut_prob_list, &cut_prob)) {
      return NULL;
    }
  }

  if (PyInt_Check(sample)) {
    /* samples chosen randomly */
    long int ns = PyInt_AsLong(sample);
    if (igraph_motifs_randesu_estimate(&self->g, &result, (igraph_integer_t) size,
          &cut_prob, (igraph_integer_t) ns, 0)) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(&cut_prob);
      return NULL;
    }
  } else {
    /* samples given in advance */
    igraph_vector_t samp;
    if (igraphmodule_PyObject_to_vector_t(sample, &samp, 1, 0)) {
      igraph_vector_destroy(&cut_prob);
      return NULL;
    }
    if (igraph_motifs_randesu_estimate(&self->g, &result, (igraph_integer_t) size,
          &cut_prob, 0, &samp)) {
      igraphmodule_handle_igraph_error();
      igraph_vector_destroy(&cut_prob);
      return NULL;
    }
  }
  igraph_vector_destroy(&cut_prob);

  return PyInt_FromLong((long)result);
}

/** \ingroup python_interface_graph
 * \brief Calculates the triad census of the graph
 * \return the triad census as a list
 * \sa igraph_triad_census
 */
PyObject *igraphmodule_Graph_triad_census(igraphmodule_GraphObject *self) {
  igraph_vector_t result;
  PyObject *list;

  if (igraph_vector_init(&result, 16)) {
    return igraphmodule_handle_igraph_error();
  }
  if (igraph_triad_census(&self->g, &result)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    return NULL;
  }

  list = igraphmodule_vector_t_to_PyTuple(&result);
  igraph_vector_destroy(&result);

  return list;
}

/**********************************************************************
 * Graph layout algorithms                                            *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Places the vertices of a graph uniformly on a circle.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_circle
 */
PyObject *igraphmodule_Graph_layout_circle(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  igraph_matrix_t m;
  int ret;
  long dim = 2;
  PyObject *result;
  static char *kwlist[] = { "dim", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|l", kwlist, &dim))
    return NULL;

  if (dim != 2 && dim != 3) {
    PyErr_SetString(PyExc_ValueError, "number of dimensions must be either 2 or 3");
    return NULL;
  }

  if (igraph_matrix_init(&m, 1, 1)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (dim == 2)
    ret = igraph_layout_circle(&self->g, &m);
  else
    ret = igraph_layout_sphere(&self->g, &m);

  if (ret) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);

  igraph_matrix_destroy(&m);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices of a graph randomly.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_random
 */
PyObject *igraphmodule_Graph_layout_random(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  igraph_matrix_t m;
  int ret;
  long dim = 2;
  PyObject *result;
  static char *kwlist[] = { "dim", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|l", kwlist, &dim))
    return NULL;

  if (dim != 2 && dim != 3) {
    PyErr_SetString(PyExc_ValueError, "number of dimensions must be either 2 or 3");
    return NULL;
  }

  if (igraph_matrix_init(&m, 1, 1)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (dim == 2)
    ret = igraph_layout_random(&self->g, &m);
  else
    ret = igraph_layout_random_3d(&self->g, &m);

  if (ret) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices on a grid
 * \sa igraph_layout_grid, igraph_layout_grid_3d
 */
PyObject *igraphmodule_Graph_layout_grid(igraphmodule_GraphObject* self,
		PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "width", "height", "dim", NULL };

  igraph_matrix_t m;
  PyObject *result;
  long int width = 0, height = 0, dim = 2;
  int ret;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lll", kwlist,
        &width, &height, &dim))
    return NULL;

  if (dim == 2 && height > 0) {
    PyErr_SetString(PyExc_ValueError, "height must not be given if dim=2");
    return NULL;
  }

  if (dim != 2 && dim != 3) {
    PyErr_SetString(PyExc_ValueError, "number of dimensions must be either 2 or 3");
    return NULL;
  }

  if (igraph_matrix_init(&m, 1, 1)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (dim == 2)
    ret = igraph_layout_grid(&self->g, &m, width);
  else
    ret = igraph_layout_grid_3d(&self->g, &m, width, height);

  if (ret != IGRAPH_SUCCESS) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&m);
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices in a star-like layout
 * \sa igraph_layout_star
 */
PyObject *igraphmodule_Graph_layout_star(igraphmodule_GraphObject* self,
		PyObject *args, PyObject *kwds) {
  static char *kwlist[] =
    { "center", "order", NULL };

  igraph_matrix_t m;
  PyObject *result, *order_o = Py_None, *center_o = Py_None;
  igraph_integer_t center = 0;
  igraph_vector_t* order = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist,
        &center_o, &order_o))
    return NULL;

  if (igraph_matrix_init(&m, 1, 1)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_PyObject_to_vid(center_o, &center, &self->g))
    return NULL;

  if (order_o != Py_None) {
    order = (igraph_vector_t*)calloc(1, sizeof(igraph_vector_t));
    if (!order) {
      igraph_matrix_destroy(&m);
      PyErr_NoMemory();
      return NULL;
    }
    if (igraphmodule_PyObject_to_vector_t(order_o, order, 1, 0)) {
      igraph_matrix_destroy(&m);
      free(order);
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  if (igraph_layout_star(&self->g, &m, center, order)) {
    if (order) {
      igraph_vector_destroy(order); free(order);
    }
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices on a plane according to the Kamada-Kawai algorithm.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_kamada_kawai
 */
PyObject *igraphmodule_Graph_layout_kamada_kawai(igraphmodule_GraphObject *
                                                 self, PyObject * args,
                                                 PyObject * kwds)
{
  static char *kwlist[] =
    { "maxiter", "sigma", "initemp", "coolexp", "kkconst", "seed", "minx", "maxx",
      "miny", "maxy", "minz", "maxz", "dim", NULL };
  igraph_matrix_t m;
  igraph_bool_t use_seed=0;
  int ret;
  long int niter = 1000, dim = 2;
  double sigma, initemp, coolexp, kkconst;
  PyObject *result, *seed_o=Py_None;
  PyObject *minx_o=Py_None, *maxx_o=Py_None;
  PyObject *miny_o=Py_None, *maxy_o=Py_None;
  PyObject *minz_o=Py_None, *maxz_o=Py_None;
  igraph_vector_t *minx=0, *maxx=0;
  igraph_vector_t *miny=0, *maxy=0;
  igraph_vector_t *minz=0, *maxz=0;

#define DESTROY_VECTORS { \
  if (minx)    { igraph_vector_destroy(minx); free(minx); } \
  if (maxx)    { igraph_vector_destroy(maxx); free(maxx); } \
  if (miny)    { igraph_vector_destroy(miny); free(miny); } \
  if (maxy)    { igraph_vector_destroy(maxy); free(maxy); } \
  if (minz)    { igraph_vector_destroy(minz); free(minz); } \
  if (maxz)    { igraph_vector_destroy(maxz); free(maxz); } \
}

  sigma = igraph_vcount(&self->g);
  kkconst = sigma * sigma;
  sigma = sigma / 4.0;
  initemp = 10.0;
  coolexp = 0.99;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lddddOOOOOOOl", kwlist,
                                   &niter, &sigma, &initemp, &coolexp,
                                   &kkconst, &seed_o,
                                   &minx_o, &maxx_o,
                                   &miny_o, &maxy_o,
                                   &minz_o, &maxz_o, &dim))
    return NULL;

  if (dim != 2 && dim != 3) {
    PyErr_SetString(PyExc_ValueError, "number of dimensions must be either 2 or 3");
    return NULL;
  }

  if (seed_o == 0 || seed_o == Py_None) {
    if (igraph_matrix_init(&m, 1, 1)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else {
    use_seed=1;
	if (igraphmodule_PyList_to_matrix_t(seed_o, &m)) return NULL;
  }

  /* Convert minimum and maximum x-y-z values */
  if (igraphmodule_attrib_to_vector_t(minx_o, self, &minx, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_t(maxx_o, self, &maxx, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_t(miny_o, self, &miny, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_t(maxy_o, self, &maxy, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (dim > 2) {
    if (igraphmodule_attrib_to_vector_t(minz_o, self, &minz, ATTRIBUTE_TYPE_EDGE)) {
      igraph_matrix_destroy(&m);
      DESTROY_VECTORS;
      igraphmodule_handle_igraph_error();
      return NULL;
    }
    if (igraphmodule_attrib_to_vector_t(maxz_o, self, &maxz, ATTRIBUTE_TYPE_EDGE)) {
      igraph_matrix_destroy(&m);
      DESTROY_VECTORS;
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }
  if (dim == 2)
    ret = igraph_layout_kamada_kawai
      (&self->g, &m, (igraph_integer_t) niter, sigma, initemp, coolexp, kkconst, use_seed,
       /*bounds*/ minx, maxx, miny, maxy);
  else
    ret = igraph_layout_kamada_kawai_3d
      (&self->g, &m, (igraph_integer_t) niter, sigma, initemp, coolexp, kkconst, use_seed, 0,
       /*bounds*/ minx, maxx, miny, maxy, minz, maxz);

  DESTROY_VECTORS;

#undef DESTROY_VECTORS

  if (ret) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices on a plane according to the DrL algorithm.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_drl
 */
PyObject* igraphmodule_Graph_layout_drl(igraphmodule_GraphObject *self,
          PyObject *args, PyObject *kwds)
{
  static char *kwlist[] =
    { "weights", "seed", "fixed", "options", "dim", NULL };
  igraph_matrix_t m;
  igraph_bool_t use_seed=0;
  igraph_vector_t *weights=0;
  igraph_vector_bool_t *fixed=0;
  igraph_layout_drl_options_t options;
  PyObject *result;
  PyObject *wobj=Py_None, *fixed_o=Py_None, *seed_o=Py_None, *options_o=Py_None;
  long dim = 2;
  int retval;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOOOl", kwlist,
                                   &wobj, &seed_o, &fixed_o, &options_o, &dim))
	  return NULL;

  if (dim != 2 && dim != 3) {
    PyErr_SetString(PyExc_ValueError, "number of dimensions must be either 2 or 3");
    return NULL;
  }

  if (igraphmodule_PyObject_to_drl_options_t(options_o, &options))
    return NULL;

  if (igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_DEFAULT)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (fixed_o != 0 && fixed_o != Py_None) {
    fixed = (igraph_vector_bool_t*)malloc(sizeof(igraph_vector_bool_t));
    if (!fixed) {
	  PyErr_NoMemory();
	  return NULL;
	}
	if (igraphmodule_PyObject_to_vector_bool_t(fixed_o, fixed)) {
	  free(fixed);
	  return NULL;
	}
  }

  if (seed_o == 0 || seed_o == Py_None) {
    if (igraph_matrix_init(&m, 1, 1)) {
      igraphmodule_handle_igraph_error();
	  if (fixed) { igraph_vector_bool_destroy(fixed); free(fixed); }
      return NULL;
    }
  } else {
    if (igraphmodule_PyList_to_matrix_t(seed_o, &m)) {
	  if (fixed) { igraph_vector_bool_destroy(fixed); free(fixed); }
	  return NULL;
	}
	use_seed=1;
  }

  /* Convert the weight parameter to a vector */
  if (igraphmodule_attrib_to_vector_t(wobj, self, &weights, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
	if (fixed) { igraph_vector_bool_destroy(fixed); free(fixed); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (dim == 2) {
	retval = igraph_layout_drl(&self->g, &m, use_seed, &options, weights, fixed);
  } else {
	retval = igraph_layout_drl_3d(&self->g, &m, use_seed, &options, weights, fixed);
  }

  if (retval) {
    igraph_matrix_destroy(&m);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
	if (fixed) { igraph_vector_bool_destroy(fixed); free(fixed); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }
  if (fixed) { igraph_vector_bool_destroy(fixed); free(fixed); }
  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices on a plane according to the Fruchterman-Reingold algorithm.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_fruchterman_reingold
 */
PyObject
  *igraphmodule_Graph_layout_fruchterman_reingold(igraphmodule_GraphObject *
                                                  self, PyObject * args,
                                                  PyObject * kwds)
{
  static char *kwlist[] =
    { "weights", "maxiter", "maxdelta", "area", "coolexp", "repulserad",
      "seed", "minx", "maxx", "miny", "maxy", "minz", "maxz", "dim", NULL };
  igraph_matrix_t m;
  igraph_bool_t use_seed=0;
  igraph_vector_t *weights=0;
  igraph_vector_t *minx=0, *maxx=0;
  igraph_vector_t *miny=0, *maxy=0;
  igraph_vector_t *minz=0, *maxz=0;
  int ret;
  long int niter = 500, dim = 2;
  double maxdelta, area, coolexp, repulserad;
  PyObject *result;
  PyObject *wobj=Py_None, *seed_o=Py_None;
  PyObject *minx_o=Py_None, *maxx_o=Py_None;
  PyObject *miny_o=Py_None, *maxy_o=Py_None;
  PyObject *minz_o=Py_None, *maxz_o=Py_None;

#define DESTROY_VECTORS { \
  if (weights) { igraph_vector_destroy(weights); free(weights); } \
  if (minx)    { igraph_vector_destroy(minx); free(minx); } \
  if (maxx)    { igraph_vector_destroy(maxx); free(maxx); } \
  if (miny)    { igraph_vector_destroy(miny); free(miny); } \
  if (maxy)    { igraph_vector_destroy(maxy); free(maxy); } \
  if (minz)    { igraph_vector_destroy(minz); free(minz); } \
  if (maxz)    { igraph_vector_destroy(maxz); free(maxz); } \
}

  maxdelta = igraph_vcount(&self->g);
  coolexp = 1.5;
  repulserad = -1;
  area = -1;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OlddddOOOOOOOl", kwlist, &wobj,
                                   &niter, &maxdelta, &area, &coolexp,
                                   &repulserad, &seed_o, &minx_o, &maxx_o,
                                   &miny_o, &maxy_o, &minz_o, &maxz_o,
                                   &dim))
    return NULL;

  if (area <= 0) {
    area = maxdelta * maxdelta;
    if (dim == 3)
      area *= maxdelta;
  }
  if (repulserad <= 0)
    repulserad = area * maxdelta;

  if (dim != 2 && dim != 3) {
    PyErr_SetString(PyExc_ValueError, "number of dimensions must be either 2 or 3");
    return NULL;
  }

  if (seed_o == 0 || seed_o == Py_None) {
    if (igraph_matrix_init(&m, 1, 1)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else {
    if (igraphmodule_PyList_to_matrix_t(seed_o, &m)) return NULL;
	use_seed=1;
  }


  /* Convert the weight parameter to a vector */
  if (igraphmodule_attrib_to_vector_t(wobj, self, &weights, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  /* Convert minimum and maximum x-y-z values */
  if (igraphmodule_attrib_to_vector_t(minx_o, self, &minx, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_t(maxx_o, self, &maxx, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_t(miny_o, self, &miny, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_t(maxy_o, self, &maxy, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    DESTROY_VECTORS;
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (dim > 2) {
    if (igraphmodule_attrib_to_vector_t(minz_o, self, &minz, ATTRIBUTE_TYPE_EDGE)) {
      igraph_matrix_destroy(&m);
      DESTROY_VECTORS;
      igraphmodule_handle_igraph_error();
      return NULL;
    }
    if (igraphmodule_attrib_to_vector_t(maxz_o, self, &maxz, ATTRIBUTE_TYPE_EDGE)) {
      igraph_matrix_destroy(&m);
      DESTROY_VECTORS;
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  if (dim == 2)
    ret = igraph_layout_fruchterman_reingold
        (&self->g, &m, (igraph_integer_t) niter, maxdelta, area, coolexp, repulserad,
         use_seed, weights, minx, maxx, miny, maxy);
  else
    ret = igraph_layout_fruchterman_reingold_3d
      (&self->g, &m, (igraph_integer_t) niter, maxdelta, area, coolexp, repulserad,
       use_seed, weights, minx, maxx, miny, maxy, minz, maxz);

  DESTROY_VECTORS;

  if (ret) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

#undef DESTROY_VECTORS

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices on a plane according to the layout algorithm in
 * graphopt 0.4.1
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_graphopt
 */
PyObject *igraphmodule_Graph_layout_graphopt(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  static char *kwlist[] =
    { "niter", "node_charge", "node_mass", "spring_length", "spring_constant",
      "max_sa_movement", "seed", NULL };
  igraph_matrix_t m;
  long int niter = 500;
  double node_charge = 0.001, node_mass = 30;
  long spring_length = 0;
  double spring_constant = 1, max_sa_movement = 5;
  PyObject *result, *seed_o = Py_None;
  igraph_bool_t use_seed=0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lddlddO", kwlist,
                                   &niter, &node_charge, &node_mass,
                                   &spring_length, &spring_constant,
                                   &max_sa_movement, &seed_o))
    return NULL;

  if (seed_o == 0 || seed_o == Py_None) {
    if (igraph_matrix_init(&m, 1, 1)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else {
    use_seed=1;
    if (igraphmodule_PyList_to_matrix_t(seed_o, &m))
			return NULL;
  }

  if (igraph_layout_graphopt(&self->g, &m, (igraph_integer_t) niter,
        node_charge, node_mass, spring_length, spring_constant,
        max_sa_movement, use_seed)) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices on a plane according to the Fruchterman-Reingold grid layout algorithm.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_grid_fruchterman_reingold
 */
PyObject
  *igraphmodule_Graph_layout_grid_fruchterman_reingold
  (igraphmodule_GraphObject * self, PyObject * args, PyObject * kwds)
{
  static char *kwlist[] =
    { "weights", "maxiter", "maxdelta", "area", "coolexp", "repulserad", "cellsize",
	  "seed", NULL };
  igraph_matrix_t m;
  long int niter = 500;
  double maxdelta, area, coolexp, repulserad, cellsize;
  PyObject *result, *seed_o = Py_None, *wobj = Py_None;
  igraph_bool_t use_seed=0;
  igraph_vector_t *weights;

  maxdelta = igraph_vcount(&self->g);
  area = maxdelta * maxdelta;
  coolexp = 1.5;
  repulserad = area * igraph_vcount(&self->g); 
  cellsize = sqrt(sqrt(area));

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OldddddO", kwlist, &wobj,
                                   &niter, &maxdelta, &area, &coolexp,
                                   &repulserad, &cellsize, &seed_o))
    return NULL;

  /* Interpret the seed matrix */
  if (seed_o == 0 || seed_o == Py_None) {
    if (igraph_matrix_init(&m, 1, 1)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  } else {
    use_seed=1;
	if (igraphmodule_PyList_to_matrix_t(seed_o, &m)) return NULL;
  }

  /* Convert the weight parameter to a vector */
  if (igraphmodule_attrib_to_vector_t(wobj, self, &weights, ATTRIBUTE_TYPE_EDGE)) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_layout_grid_fruchterman_reingold
      (&self->g, &m, (igraph_integer_t) niter, maxdelta, area, coolexp, repulserad,
       cellsize, use_seed, weights)) {
    igraph_matrix_destroy(&m);
    if (weights) {
      igraph_vector_destroy(weights); free(weights);
    }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (weights) {
    igraph_vector_destroy(weights); free(weights);
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices of a graph according to the Large Graph Layout
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_lgl
 */
PyObject *igraphmodule_Graph_layout_lgl(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  static char *kwlist[] =
    { "maxiter", "maxdelta", "area", "coolexp", "repulserad", "cellsize", "root",
    NULL };
  igraph_matrix_t m;
  PyObject *result, *root_o = Py_None;
  long int maxiter = 150;
  igraph_integer_t proot = -1;
  double maxdelta, area, coolexp, repulserad, cellsize;

  maxdelta = igraph_vcount(&self->g);
  area = -1;
  coolexp = 1.5;
  repulserad = -1;
  cellsize = -1;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ldddddO", kwlist,
                                   &maxiter, &maxdelta, &area, &coolexp,
                                   &repulserad, &cellsize, &root_o))
    return NULL;

  if (area <= 0)
    area = igraph_vcount(&self->g)*igraph_vcount(&self->g);
  if (repulserad <= 0)
    repulserad = area*igraph_vcount(&self->g);
  if (cellsize <= 0)
    cellsize = sqrt(sqrt(area));

  if (igraphmodule_PyObject_to_vid(root_o, &proot, &self->g))
    return NULL;

  if (igraph_matrix_init(&m, 1, 1)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_layout_lgl(&self->g, &m, (igraph_integer_t) maxiter, maxdelta,
                        area, coolexp, repulserad, cellsize, proot)) {
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices of a graph using multidimensional scaling
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_mds
 */
PyObject *igraphmodule_Graph_layout_mds(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  static char *kwlist[] =
    { "dist", "dim", "arpack_options", NULL };
  igraph_matrix_t m;
  igraph_matrix_t *dist = 0;
  long int dim = 2;
  igraphmodule_ARPACKOptionsObject *arpack_options;
  PyObject *dist_o = Py_None;
  PyObject *arpack_options_o = igraphmodule_arpack_options_default;
  PyObject *result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OlO!", kwlist, &dist_o,
                                   &dim, &igraphmodule_ARPACKOptionsType,
                                   &arpack_options_o))
    return NULL;

  if (dist_o != Py_None) {
    dist = (igraph_matrix_t*)malloc(sizeof(igraph_matrix_t));
    if (!dist) {
      PyErr_NoMemory();
      return NULL;
    }
	if (igraphmodule_PyList_to_matrix_t(dist_o, dist)) {
      free(dist);
      return NULL;
    }
  }

  if (igraph_matrix_init(&m, 1, 1)) {
    if (dist) {
      igraph_matrix_destroy(dist); free(dist);
    }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  arpack_options = (igraphmodule_ARPACKOptionsObject*)arpack_options_o;
  if (igraph_layout_mds(&self->g, &m, dist, dim,
                        igraphmodule_ARPACKOptions_get(arpack_options))) {
    if (dist) {
      igraph_matrix_destroy(dist); free(dist);
    }
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (dist) {
    igraph_matrix_destroy(dist); free(dist);
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices of a graph according to the Reingold-Tilford
 * tree layout algorithm
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_reingold_tilford
 */
PyObject *igraphmodule_Graph_layout_reingold_tilford(igraphmodule_GraphObject
                                                     * self, PyObject * args,
                                                     PyObject * kwds)
{
  static char *kwlist[] = { "mode", "root", "rootlevel", NULL };
  igraph_matrix_t m;
  igraph_vector_t roots, *roots_p = 0;
  igraph_vector_t rootlevels, *rootlevels_p = 0;
  PyObject *roots_o=Py_None, *rootlevels_o=Py_None, *mode_o=Py_None;
  igraph_neimode_t mode = IGRAPH_OUT;
  PyObject *result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist,
    &mode_o, &roots_o, &rootlevels_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode))
    return NULL;

  if (roots_o != Py_None) {
    roots_p = &roots;
    if (igraphmodule_PyObject_to_vector_t(roots_o, roots_p, 1, 0)) return 0;
  }
  if (rootlevels_o != Py_None) {
    rootlevels_p = &rootlevels;
    if (igraphmodule_PyObject_to_vector_t(rootlevels_o, rootlevels_p, 1, 0)) {
      if (roots_p) igraph_vector_destroy(roots_p);
      return 0;
    }
  }

  if (igraph_matrix_init(&m, 1, 1)) {
    if (roots_p) igraph_vector_destroy(roots_p);
    if (rootlevels_p) igraph_vector_destroy(rootlevels_p);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_layout_reingold_tilford(&self->g, &m, mode, roots_p, rootlevels_p)) {
    igraph_matrix_destroy(&m);
    if (roots_p) igraph_vector_destroy(roots_p);
    if (rootlevels_p) igraph_vector_destroy(rootlevels_p);
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (roots_p) igraph_vector_destroy(roots_p);
  if (rootlevels_p) igraph_vector_destroy(rootlevels_p);

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices of a graph according to the Reingold-Tilford
 * tree layout algorithm in a polar coordinate system
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_reingold_tilford
 */
PyObject *igraphmodule_Graph_layout_reingold_tilford_circular(
  igraphmodule_GraphObject * self, PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "mode", "root", "rootlevel", NULL };
  igraph_matrix_t m;
  igraph_vector_t roots, *roots_p = 0;
  igraph_vector_t rootlevels, *rootlevels_p = 0;
  PyObject *roots_o=Py_None, *rootlevels_o=Py_None, *mode_o=Py_None;
  igraph_neimode_t mode = IGRAPH_OUT;
  PyObject *result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist,
    &mode_o, &roots_o, &rootlevels_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;

  if (roots_o != Py_None) {
    roots_p = &roots;
    if (igraphmodule_PyObject_to_vector_t(roots_o, roots_p, 1, 0)) return 0;
  }
  if (rootlevels_o != Py_None) {
    rootlevels_p = &rootlevels;
    if (igraphmodule_PyObject_to_vector_t(rootlevels_o, rootlevels_p, 1, 0)) {
      if (roots_p) igraph_vector_destroy(roots_p);
      return 0;
    }
  }

  if (igraph_matrix_init(&m, 1, 1)) {
    if (roots_p) igraph_vector_destroy(roots_p);
    if (rootlevels_p) igraph_vector_destroy(rootlevels_p);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_layout_reingold_tilford_circular(&self->g, &m, mode, roots_p,
      rootlevels_p)) {
    igraph_matrix_destroy(&m);
    if (roots_p) igraph_vector_destroy(roots_p);
    if (rootlevels_p) igraph_vector_destroy(rootlevels_p);
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (roots_p) igraph_vector_destroy(roots_p);
  if (rootlevels_p) igraph_vector_destroy(rootlevels_p);

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);
  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Places the vertices of a graph according to the Sugiyama layout.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_sugiyama
 */
PyObject *igraphmodule_Graph_layout_sugiyama(
  igraphmodule_GraphObject * self, PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "layers", "weights", "hgap", "vgap", "maxiter",
    "return_extended_graph", NULL };
  igraph_matrix_t m;
  igraph_t extd_graph;
  igraph_vector_t extd_to_orig_eids;
  igraph_vector_t *weights = 0, *layers = 0;
  double hgap = 1, vgap = 1;
  long int maxiter = 100;
  PyObject *layers_o = Py_None, *weights_o = Py_None, *extd_to_orig_eids_o = Py_None;
  PyObject *return_extended_graph = Py_False;
  PyObject *result;
  igraphmodule_GraphObject *graph_o;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOddlO", kwlist,
    &layers_o, &weights_o, &hgap, &vgap, &maxiter, &return_extended_graph))
    return NULL;

  if (igraph_vector_init(&extd_to_orig_eids, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_matrix_init(&m, 1, 1)) {
    igraph_vector_destroy(&extd_to_orig_eids);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(layers_o, self, &layers,
	  ATTRIBUTE_TYPE_VERTEX)) {
    igraph_vector_destroy(&extd_to_orig_eids);
    igraph_matrix_destroy(&m);
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (layers != 0) { igraph_vector_destroy(layers); free(layers); }
    igraph_vector_destroy(&extd_to_orig_eids);
    igraph_matrix_destroy(&m);
    return NULL;
  }

  if (igraph_layout_sugiyama(&self->g, &m,
        (PyObject_IsTrue(return_extended_graph) ? &extd_graph : 0),
        (PyObject_IsTrue(return_extended_graph) ? &extd_to_orig_eids : 0),
        layers, hgap, vgap, maxiter, weights)) {
    if (layers != 0) { igraph_vector_destroy(layers); free(layers); }
    if (weights != 0) { igraph_vector_destroy(weights); free(weights); }
    igraph_vector_destroy(&extd_to_orig_eids);
    igraph_matrix_destroy(&m);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (layers != 0) { igraph_vector_destroy(layers); free(layers); }
  if (weights != 0) { igraph_vector_destroy(weights); free(weights); }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  igraph_matrix_destroy(&m);

  if (PyObject_IsTrue(return_extended_graph)) {
    CREATE_GRAPH(graph_o, extd_graph);
    extd_to_orig_eids_o = igraphmodule_vector_t_to_PyList(&extd_to_orig_eids,
        IGRAPHMODULE_TYPE_INT);
    result = Py_BuildValue("NNN", result, graph_o, extd_to_orig_eids_o);
  }

  igraph_vector_destroy(&extd_to_orig_eids);
  return (PyObject *) result;
}

/**********************************************************************
 * Conversion between various graph representations                   *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Returns the adjacency matrix of a graph.
 * \return the adjacency matrix as a Python list of lists
 * \sa igraph_get_adjacency
 */
PyObject *igraphmodule_Graph_get_adjacency(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "type", "eids", NULL };
  igraph_get_adjacency_t t = IGRAPH_GET_ADJACENCY_BOTH;
  igraph_matrix_t m;
  PyObject *result, *eids = Py_False;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iO", kwlist, &t, &eids))
    return NULL;

  if (t != IGRAPH_GET_ADJACENCY_UPPER && t != IGRAPH_GET_ADJACENCY_LOWER &&
      t != IGRAPH_GET_ADJACENCY_BOTH) {
    PyErr_SetString(PyExc_ValueError,
                    "type must be either GET_ADJACENCY_LOWER or GET_ADJACENCY_UPPER or GET_ADJACENCY_BOTH");
    return NULL;
  }

  if (igraph_matrix_init
      (&m, igraph_vcount(&self->g), igraph_vcount(&self->g))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_get_adjacency(&self->g, &m, t, PyObject_IsTrue(eids))) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&m);
    return NULL;
  }

  result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_INT);
  igraph_matrix_destroy(&m);
  return result;
}

/** \ingroup python_interface_graph
 * \brief Returns the incidence matrix of a bipartite graph.
 * \return the incidence matrix as a Python list of lists
 * \sa igraph_get_incidence
 */
PyObject *igraphmodule_Graph_get_incidence(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "types", NULL };
  igraph_matrix_t matrix;
  igraph_vector_t row_ids, col_ids;
  igraph_vector_bool_t *types;
  PyObject *matrix_o, *row_ids_o, *col_ids_o, *types_o;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &types_o))
    return NULL;

  if (igraph_vector_init(&row_ids, 0))
	return NULL;

  if (igraph_vector_init(&col_ids, 0)) {
    igraph_vector_destroy(&row_ids);
	return NULL;
  }

  if (igraphmodule_attrib_to_vector_bool_t(types_o, self, &types, ATTRIBUTE_TYPE_VERTEX)) {
    igraph_vector_destroy(&row_ids);
    igraph_vector_destroy(&col_ids);
	return NULL;
  }

  if (igraph_matrix_init(&matrix, 1, 1)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&row_ids);
    igraph_vector_destroy(&col_ids);
    if (types) { igraph_vector_bool_destroy(types); free(types); }
    return NULL;
  }

  if (igraph_get_incidence(&self->g, types, &matrix, &row_ids, &col_ids)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&row_ids);
    igraph_vector_destroy(&col_ids);
    if (types) { igraph_vector_bool_destroy(types); free(types); }
    igraph_matrix_destroy(&matrix);
    return NULL;
  }
  
  if (types) { igraph_vector_bool_destroy(types); free(types); }

  matrix_o = igraphmodule_matrix_t_to_PyList(&matrix, IGRAPHMODULE_TYPE_INT);
  igraph_matrix_destroy(&matrix);

  row_ids_o = igraphmodule_vector_t_to_PyList(&row_ids, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&row_ids);
  col_ids_o = igraphmodule_vector_t_to_PyList(&col_ids, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&col_ids);

  return Py_BuildValue("NNN", matrix_o, row_ids_o, col_ids_o);
}

/** \ingroup python_interface_graph
 * \brief Returns the Laplacian matrix of a graph.
 * \return the Laplacian matrix as a Python list of lists
 * \sa igraph_laplacian
 */
PyObject *igraphmodule_Graph_laplacian(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "weights", "normalized", NULL };
  igraph_matrix_t m;
  PyObject *result;
  PyObject *weights_o = Py_None;
  PyObject *normalized = Py_False;
  igraph_vector_t *weights = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist,
        &weights_o, &normalized))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  if (igraph_matrix_init
      (&m, igraph_vcount(&self->g), igraph_vcount(&self->g))) {
    igraphmodule_handle_igraph_error();
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }

  if (igraph_laplacian(&self->g, &m, /*sparseres=*/ 0, 
		       PyObject_IsTrue(normalized), weights)) {
    igraphmodule_handle_igraph_error();
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    igraph_matrix_destroy(&m);
    return NULL;
  }

  if (PyObject_IsTrue(normalized) || weights) {
    result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
  }
  else {
    result = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_INT);
  }

  if (weights) { igraph_vector_destroy(weights); free(weights); }
  igraph_matrix_destroy(&m);

  return result;
}

/** \ingroup python_interface_graph
 * \brief Returns the list of edges in a graph.
 * \return the list of edges, every edge is represented by a pair
 * \sa igraph_get_edgelist
 */
PyObject *igraphmodule_Graph_get_edgelist(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  igraph_vector_t edgelist;
  PyObject *result;

  igraph_vector_init(&edgelist, igraph_ecount(&self->g));
  if (igraph_get_edgelist(&self->g, &edgelist, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&edgelist);
    return NULL;
  }

  result = igraphmodule_vector_t_to_PyList_pairs(&edgelist);
  igraph_vector_destroy(&edgelist);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \function igraphmodule_Graph_to_undirected
 * \brief Converts a directed graph to an undirected one.
 * \return The undirected graph.
 * \sa igraph_to_undirected
 */
PyObject *igraphmodule_Graph_to_undirected(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  PyObject *mode_o = Py_None, *comb_o = Py_None;
  igraph_to_undirected_t mode = IGRAPH_TO_UNDIRECTED_COLLAPSE;
  igraph_attribute_combination_t comb;
  static char *kwlist[] = { "mode", "combine_edges", NULL };

  if (kwds && PyDict_GetItemString(kwds, "mode") && PyDict_GetItemString(kwds, "collapse"))
    PyErr_SetString(PyExc_ValueError, "cannot specify mode=... and collapse=... at the same time");

  if (kwds && PyDict_GetItemString(kwds, "collapse")) {
    mode_o = PyDict_GetItemString(kwds, "collapse");
    PY_IGRAPH_DEPRECATED("The collapse=... keyword argument of Graph.to_undirected() is deprecated, please use mode=... instead. This warning will be removed from igraph 0.7.");
    PyDict_DelItemString(kwds, "collapse");
  }

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &mode_o, &comb_o))
    return NULL;

  if (igraphmodule_PyObject_to_to_undirected_t(mode_o, &mode))
    return NULL;

  if (igraphmodule_PyObject_to_attribute_combination_t(comb_o, &comb))
	return NULL;

  if (igraph_to_undirected(&self->g, mode, &comb)) {
    igraph_attribute_combination_destroy(&comb);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  igraph_attribute_combination_destroy(&comb);

  Py_RETURN_NONE;
}


/** \ingroup python_interface_graph
 * \function igraphmodule_Graph_to_directed
 * \brief Converts an undirected graph to a directed one.
 * \return The directed graph.
 * \sa igraph_to_directed
 */
PyObject *igraphmodule_Graph_to_directed(igraphmodule_GraphObject * self,
                                         PyObject * args, PyObject * kwds)
{
  PyObject *mutual = Py_True;
  igraph_to_directed_t mode = IGRAPH_TO_DIRECTED_MUTUAL;
  static char *kwlist[] = { "mutual", NULL };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &mutual))
    return NULL;
  mode =
    (PyObject_IsTrue(mutual) ? IGRAPH_TO_DIRECTED_MUTUAL :
     IGRAPH_TO_DIRECTED_ARBITRARY);
  if (igraph_to_directed(&self->g, mode)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  Py_RETURN_NONE;
}

/**********************************************************************
 * Reading/writing foreing graph formats                              *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Reads a DIMACS file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_dimacs
 */
PyObject *igraphmodule_Graph_Read_DIMACS(PyTypeObject * type,
                                         PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraphmodule_filehandle_t fobj;
  igraph_integer_t source = 0, target = 0;
  igraph_vector_t capacity;
  igraph_t g;
  PyObject *fname = NULL, *directed = Py_False, *capacity_obj;

  static char *kwlist[] = { "f", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "O|O", kwlist, &fname, &directed))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_vector_init(&capacity, 0)) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  if (igraph_read_graph_dimacs(&g, igraphmodule_filehandle_get(&fobj),
                               0, 0, &source, &target,
                               &capacity, PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&capacity);
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  igraphmodule_filehandle_destroy(&fobj);

  capacity_obj = igraphmodule_vector_t_to_PyList(&capacity, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&capacity);
  if (!capacity_obj)
    return NULL;

  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return Py_BuildValue("NiiN", (PyObject *) self, (long)source,
                       (long)target, capacity_obj);
}

/** \ingroup python_interface_graph
 * \brief Reads an UCINET DL file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_dl
 */
PyObject *igraphmodule_Graph_Read_DL(PyTypeObject * type,
                                     PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  igraph_t g;
  igraphmodule_filehandle_t fobj;
  PyObject *fname = NULL, *directed = Py_True;

  static char *kwlist[] = { "f", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "O|O", kwlist, &fname, &directed))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_dl(&g, igraphmodule_filehandle_get(&fobj),
        PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject*)self;
}

/** \ingroup python_interface_graph
 * \brief Reads an edge list from a file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_edgelist
 */
PyObject *igraphmodule_Graph_Read_Edgelist(PyTypeObject * type,
                                           PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  PyObject *directed = Py_True, *fname = NULL;
  igraphmodule_filehandle_t fobj;
  igraph_t g;

  static char *kwlist[] = { "f", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist,
                                   &fname, &directed))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_edgelist(&g, igraphmodule_filehandle_get(&fobj),
        0, PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Reads an edge list from an NCOL file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_ncol
 */
PyObject *igraphmodule_Graph_Read_Ncol(PyTypeObject * type, PyObject * args,
                                       PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  PyObject *names = Py_True, *weights = Py_None, *directed = Py_True;
  PyObject *fname = NULL;
  igraphmodule_filehandle_t fobj;
  igraph_add_weights_t add_weights = IGRAPH_ADD_WEIGHTS_IF_PRESENT;
  igraph_t g;

  static char *kwlist[] = { "f", "names", "weights", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OOO", kwlist,
                                   &fname, &names, &weights, &directed))
    return NULL;

  if (igraphmodule_PyObject_to_add_weights_t(weights, &add_weights))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_ncol(&g, igraphmodule_filehandle_get(&fobj), 0,
      PyObject_IsTrue(names), add_weights,
      PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Reads an edge list from an LGL file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_lgl
 */
PyObject *igraphmodule_Graph_Read_Lgl(PyTypeObject * type, PyObject * args,
                                      PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  PyObject *names = Py_True, *weights = Py_None, *directed = Py_True;
  PyObject *fname = NULL;
  igraphmodule_filehandle_t fobj;
  igraph_add_weights_t add_weights = IGRAPH_ADD_WEIGHTS_IF_PRESENT;
  igraph_t g;

  static char *kwlist[] = { "f", "names", "weights", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OOO", kwlist,
                                   &fname, &names, &weights, &directed))
    return NULL;

  if (igraphmodule_PyObject_to_add_weights_t(weights, &add_weights))
    return NULL;

  if (kwds && PyDict_Check(kwds) && \
      PyDict_GetItemString(kwds, "directed") == NULL) {
    if (PyErr_Occurred())
      return NULL;
    PY_IGRAPH_WARN("Graph.Read_Lgl creates directed networks by default from igraph 0.6. To get rid of this warning, specify directed=... explicitly. This warning will be removed from igraph 0.7.");
  }

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_lgl(&g, igraphmodule_filehandle_get(&fobj),
        PyObject_IsTrue(names), add_weights,
        PyObject_IsTrue(directed))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Reads an edge list from a Pajek file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_pajek
 */
PyObject *igraphmodule_Graph_Read_Pajek(PyTypeObject * type, PyObject * args,
                                        PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  PyObject *fname = NULL;
  igraphmodule_filehandle_t fobj;
  igraph_t g;

  static char *kwlist[] = { "f", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &fname))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_pajek(&g, igraphmodule_filehandle_get(&fobj))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  
  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Reads a GML file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_gml
 */
PyObject *igraphmodule_Graph_Read_GML(PyTypeObject * type,
                                      PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  PyObject *fname = NULL;
  igraphmodule_filehandle_t fobj;
  igraph_t g;

  static char *kwlist[] = { "f", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &fname))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_gml(&g, igraphmodule_filehandle_get(&fobj))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);

  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Reads a GraphDB file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_graphdb
 */
PyObject *igraphmodule_Graph_Read_GraphDB(PyTypeObject * type,
                                          PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  PyObject *fname = NULL, *directed_o = Py_False;
  igraph_t g;
  igraphmodule_filehandle_t fobj;

  static char *kwlist[] = { "f", "directed", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &fname, &directed_o))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_graphdb(&g, igraphmodule_filehandle_get(&fobj),
        PyObject_IsTrue(directed_o))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  
  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);
  
  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Reads a GraphML file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_graphml
 */
PyObject *igraphmodule_Graph_Read_GraphML(PyTypeObject * type,
                                          PyObject * args, PyObject * kwds)
{
  igraphmodule_GraphObject *self;
  PyObject *fname = NULL;
  long int index = 0;
  igraph_t g;
  igraphmodule_filehandle_t fobj;

  static char *kwlist[] = { "f", "index", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|l", kwlist, &fname, &index))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "r"))
    return NULL;

  if (igraph_read_graph_graphml(&g, igraphmodule_filehandle_get(&fobj),
        (igraph_integer_t) index)) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  
  igraphmodule_filehandle_destroy(&fobj);
  CREATE_GRAPH_FROM_TYPE(self, g, type);
  
  return (PyObject *) self;
}

/** \ingroup python_interface_graph
 * \brief Writes the graph as a DIMACS file
 * \return none
 * \sa igraph_write_graph_dimacs
 */
PyObject *igraphmodule_Graph_write_dimacs(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  long source = 0, target = 0;
  PyObject *capacity_obj = Py_None, *fname = NULL;
  igraphmodule_filehandle_t fobj;
  igraph_vector_t* capacity = 0;

  static char *kwlist[] = {
    "f", "source", "target", "capacity", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oll|O", kwlist, &fname,
                                   &source, &target, &capacity_obj))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (capacity_obj == Py_None) {
    capacity_obj = PyString_FromString("capacity");
  } else {
    Py_INCREF(capacity_obj);
  }

  if (igraphmodule_attrib_to_vector_t(capacity_obj, self, &capacity, ATTRIBUTE_TYPE_EDGE)) {
    igraphmodule_filehandle_destroy(&fobj);
    Py_DECREF(capacity_obj);
    return NULL;
  }

  Py_DECREF(capacity_obj);

  if (igraph_write_graph_dimacs(&self->g, igraphmodule_filehandle_get(&fobj),
        source, target, capacity)) {
    igraphmodule_handle_igraph_error();
    if (capacity) {
      igraph_vector_destroy(capacity); free(capacity);
    }
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }

  if (capacity) {
    igraph_vector_destroy(capacity); free(capacity);
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}


/** \ingroup python_interface_graph
 * \brief Writes the graph as a DOT (GraphViz) file
 * \return none
 * \sa igraph_write_graph_dot
 */
PyObject *igraphmodule_Graph_write_dot(igraphmodule_GraphObject * self,
  PyObject * args, PyObject * kwds) {
  PyObject *fname = NULL;
  igraphmodule_filehandle_t fobj;
  static char *kwlist[] = { "f", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &fname))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (igraph_write_graph_dot(&self->g, igraphmodule_filehandle_get(&fobj))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Writes the edge list to a file
 * \return none
 * \sa igraph_write_graph_edgelist
 */
PyObject *igraphmodule_Graph_write_edgelist(igraphmodule_GraphObject * self,
                                            PyObject * args, PyObject * kwds)
{
  PyObject *fname = NULL;
  igraphmodule_filehandle_t fobj;
  static char *kwlist[] = { "f", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &fname))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (igraph_write_graph_edgelist(&self->g, igraphmodule_filehandle_get(&fobj))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}


/** \ingroup python_interface_graph
 * \brief Writes the graph as a GML file
 * \return none
 * \sa igraph_write_graph_gml
 */
PyObject *igraphmodule_Graph_write_gml(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  PyObject *ids = Py_None, *fname = NULL;
  PyObject *creator = Py_None;
  igraph_vector_t idvec, *idvecptr=0;
  char *creator_str=0;
  igraphmodule_filehandle_t fobj;

  static char *kwlist[] = {
    "f", "creator", "ids", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OO", kwlist, &fname, &creator, &ids))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (PyList_Check(ids)) {
    idvecptr = &idvec;
    if (igraphmodule_PyObject_to_vector_t(ids, idvecptr, 0, 0)) {
      igraphmodule_filehandle_destroy(&fobj);
      return NULL;
    }
  }

  if (creator != Py_None) {
    PyObject* o = PyObject_Str(creator);
    if (o == 0) {
      if (idvecptr)
        igraph_vector_destroy(idvecptr);
      igraphmodule_filehandle_destroy(&fobj);
    }

    creator_str = PyString_CopyAsString(o);
    Py_DECREF(o);

    if (creator_str == 0) {
      if (idvecptr)
        igraph_vector_destroy(idvecptr);
      igraphmodule_filehandle_destroy(&fobj);
    }
  }

  if (igraph_write_graph_gml(&self->g, igraphmodule_filehandle_get(&fobj),
        idvecptr, creator_str)) {
    if (idvecptr) { igraph_vector_destroy(idvecptr); }
    if (creator_str)
      free(creator_str);
    igraphmodule_filehandle_destroy(&fobj);
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  if (idvecptr) { igraph_vector_destroy(idvecptr); }
  if (creator_str)
    free(creator_str);
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Writes the edge list to a file in .ncol format
 * \return none
 * \sa igraph_write_graph_ncol
 */
PyObject *igraphmodule_Graph_write_ncol(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  PyObject *fname = NULL;
  char *names = "name";
  char *weights = "weight";
  igraphmodule_filehandle_t fobj;

  static char *kwlist[] = { "f", "names", "weights", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|zz", kwlist,
                                   &fname, &names, &weights))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (igraph_write_graph_ncol(&self->g, igraphmodule_filehandle_get(&fobj),
        names, weights)) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Writes the edge list to a file in .lgl format
 * \return none
 * \sa igraph_write_graph_lgl
 */
PyObject *igraphmodule_Graph_write_lgl(igraphmodule_GraphObject * self,
                                       PyObject * args, PyObject * kwds)
{
  PyObject *fname = NULL;
  char *names = "name";
  char *weights = "weight";
  PyObject *isolates = Py_True;
  igraphmodule_filehandle_t fobj;

  static char *kwlist[] = { "f", "names", "weights", "isolates", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|zzO", kwlist,
                                   &fname, &names, &weights, &isolates))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (igraph_write_graph_lgl(&self->g, igraphmodule_filehandle_get(&fobj),
        names, weights, PyObject_IsTrue(isolates))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Writes the graph as a Pajek .net file
 * \return none
 * \sa igraph_write_graph_pajek
 */
PyObject *igraphmodule_Graph_write_pajek(igraphmodule_GraphObject * self,
  PyObject * args, PyObject * kwds) {
  PyObject *fname = NULL;
  static char *kwlist[] = { "f", NULL };
  igraphmodule_filehandle_t fobj;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &fname))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (igraph_write_graph_pajek(&self->g, igraphmodule_filehandle_get(&fobj))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Writes the graph to a GraphML file
 * \return none
 * \sa igraph_write_graph_graphml
 */
PyObject *igraphmodule_Graph_write_graphml(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  PyObject *fname = NULL;
  static char *kwlist[] = { "f", NULL };
  igraphmodule_filehandle_t fobj;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &fname))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (igraph_write_graph_graphml(&self->g, igraphmodule_filehandle_get(&fobj))) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}

/** \ingroup python_interface_graph
 * \brief Writes the edge list to a file in LEDA native format
 * \return none
 * \sa igraph_write_graph_leda
 */
PyObject *igraphmodule_Graph_write_leda(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  PyObject *fname = NULL;
  char *vertex_attr_name = "name";
  char *edge_attr_name = "weight";
  igraphmodule_filehandle_t fobj;

  static char *kwlist[] = { "f", "names", "weights", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|zz", kwlist,
                                   &fname, &vertex_attr_name,
								   &edge_attr_name))
    return NULL;

  if (igraphmodule_filehandle_init(&fobj, fname, "w"))
    return NULL;

  if (igraph_write_graph_leda(&self->g, igraphmodule_filehandle_get(&fobj),
        vertex_attr_name, edge_attr_name)) {
    igraphmodule_handle_igraph_error();
    igraphmodule_filehandle_destroy(&fobj);
    return NULL;
  }
  igraphmodule_filehandle_destroy(&fobj);

  Py_RETURN_NONE;
}

/**********************************************************************
 * Routines related to graph isomorphism                              *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Calculates the isomorphy class of a graph or its subgraph
 * \sa igraph_isoclass, igraph_isoclass_subgraph
 */
PyObject *igraphmodule_Graph_isoclass(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  Py_ssize_t n;
  igraph_integer_t isoclass = 0;
  PyObject *vids = 0;
  char *kwlist[] = { "vertices", NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "|O!", kwlist, &PyList_Type, &vids))
    return NULL;

  n = vids ? PyList_Size(vids) : igraph_vcount(&self->g);
  if (n < 3 || n > 4) {
    PyErr_SetString(PyExc_ValueError,
                    "Graph or subgraph must have 3 or 4 vertices.");
    return NULL;
  }

  if (vids) {
    igraph_vector_t vidsvec;
    if (igraphmodule_PyObject_to_vector_t(vids, &vidsvec, 1, 0)) {
      PyErr_SetString(PyExc_ValueError,
                      "Error while converting PyList to igraph_vector_t");
      return NULL;
    }
    if (igraph_isoclass_subgraph(&self->g, &vidsvec, &isoclass)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }
  else {
    if (igraph_isoclass(&self->g, &isoclass)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  return PyInt_FromLong((long)isoclass);
}

/** \ingroup python_interface_graph
 * \brief Determines whether the graph is isomorphic to another graph.
 *
 * \sa igraph_isomorphic
 */
PyObject *igraphmodule_Graph_isomorphic(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  igraph_bool_t result = 0;
  PyObject *o = Py_None;
  igraphmodule_GraphObject *other;
  static char *kwlist[] = { "other", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!", kwlist,
      &igraphmodule_GraphType, &o))
    return NULL;
  if (o == Py_None) other = self; else other = (igraphmodule_GraphObject *) o;

  if (igraph_isomorphic(&self->g, &other->g, &result)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (result) Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Determines whether the graph is isomorphic to another graph,
 *   using the BLISS isomorphism algorithm
 *
 * The actual code is almost the same as igraphmodule_Graph_isomorphic_vf2.
 * Be sure to correct bugs in both interfaces if applicable!
 *
 * \sa igraph_isomorphic_bliss
 */
PyObject *igraphmodule_Graph_isomorphic_bliss(igraphmodule_GraphObject * self,
                                              PyObject * args, PyObject * kwds)
{
  igraph_bool_t result = 0;
  PyObject *o=Py_None, *return1=Py_False, *return2=Py_False;
  PyObject *sho1=Py_None, *sho2=Py_None;
  igraphmodule_GraphObject *other;
  igraph_vector_t mapping_12, mapping_21, *map12=0, *map21=0;
  igraph_bliss_sh_t sh1=IGRAPH_BLISS_FM, sh2=IGRAPH_BLISS_FM;

  static char *kwlist[] = { "other", "return_mapping_12",
    "return_mapping_21", "sh1", "sh2", NULL };
  /* TODO: convert igraph_bliss_info_t when needed */
  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "|O!OOOO", kwlist, &igraphmodule_GraphType, &o,
       &return1, &return2, &sho1, &sho2))
    return NULL;
  if (igraphmodule_PyObject_to_bliss_sh_t(sho1, &sh1)) return NULL;
  if (igraphmodule_PyObject_to_bliss_sh_t(sho2, &sh2)) return NULL;
  if (o == Py_None) other = self; else other = (igraphmodule_GraphObject *) o;

  if (PyObject_IsTrue(return1)) {
	igraph_vector_init(&mapping_12, 0);
	map12 = &mapping_12;
  }
  if (PyObject_IsTrue(return2)) {
	igraph_vector_init(&mapping_21, 0);
	map21 = &mapping_21;
  }

  if (igraph_isomorphic_bliss(&self->g, &other->g, &result, map12, map21,
      sh1, sh2, 0, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (!map12 && !map21) {
    if (result) Py_RETURN_TRUE;
    Py_RETURN_FALSE;
  } else {
	PyObject *iso, *m1, *m2;
	iso = result ? Py_True : Py_False;
	Py_INCREF(iso);
	if (map12) {
	  m1 = igraphmodule_vector_t_to_PyList(map12, IGRAPHMODULE_TYPE_INT);
	  igraph_vector_destroy(map12);
	  if (!m1) {
	    Py_DECREF(iso);
		if (map21) igraph_vector_destroy(map21);
		return NULL;
	  }
	} else { m1 = Py_None; Py_INCREF(m1); }
	if (map21) {
	  m2 = igraphmodule_vector_t_to_PyList(map21, IGRAPHMODULE_TYPE_INT);
	  igraph_vector_destroy(map21);
	  if (!m2) {
	    Py_DECREF(iso); Py_DECREF(m1);
		return NULL;
	  }
	} else { m2 = Py_None; Py_INCREF(m2); }
	return Py_BuildValue("NNN", iso, m1, m2);
  }
}


typedef struct {
  PyObject* node_compat_fn;
  PyObject* edge_compat_fn;
  PyObject* callback_fn;
  PyObject* graph1;
  PyObject* graph2;
} igraphmodule_i_Graph_isomorphic_vf2_callback_data_t;

igraph_bool_t igraphmodule_i_Graph_isomorphic_vf2_callback_fn(
    const igraph_vector_t *map12, const igraph_vector_t *map21,
    void* extra) {
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t* data =
    (igraphmodule_i_Graph_isomorphic_vf2_callback_data_t*)extra;
  igraph_bool_t retval;
  PyObject *map12_o, *map21_o;
  PyObject *result;

  map12_o = igraphmodule_vector_t_to_PyList(map12, IGRAPHMODULE_TYPE_INT);
  if (map12_o == NULL) {
    /* Error in conversion, return 0 to stop the search */
    PyErr_WriteUnraisable(data->callback_fn);
    return 0;
  }

  map21_o = igraphmodule_vector_t_to_PyList(map21, IGRAPHMODULE_TYPE_INT);
  if (map21_o == NULL) {
    /* Error in conversion, return 0 to stop the search */
    PyErr_WriteUnraisable(data->callback_fn);
    Py_DECREF(map21_o);
    return 0;
  }

  result = PyObject_CallFunction(data->callback_fn, "OOOO", data->graph1, data->graph2,
      map12_o, map21_o);
  Py_DECREF(map12_o);
  Py_DECREF(map21_o);

  if (result == NULL) {
    /* Error in callback, return 0 */
    PyErr_WriteUnraisable(data->callback_fn);
    return 0;
  }

  retval = PyObject_IsTrue(result);
  Py_DECREF(result);

  return retval;
}

igraph_bool_t igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn(
    const igraph_t *graph1, const igraph_t *graph2,
    const igraph_integer_t cand1, const igraph_integer_t cand2,
    void* extra) {
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t* data =
    (igraphmodule_i_Graph_isomorphic_vf2_callback_data_t*)extra;
  igraph_bool_t retval;
  PyObject *result;

  result = PyObject_CallFunction(data->node_compat_fn, "OOll",
      data->graph1, data->graph2, (long)cand1, (long)cand2);

  if (result == NULL) {
    /* Error in callback, return 0 */
    PyErr_WriteUnraisable(data->node_compat_fn);
    return 0;
  }

  retval = PyObject_IsTrue(result);
  Py_DECREF(result);

  return retval;
}

igraph_bool_t igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn(
    const igraph_t *graph1, const igraph_t *graph2,
    const igraph_integer_t cand1, const igraph_integer_t cand2,
    void* extra) {
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t* data =
    (igraphmodule_i_Graph_isomorphic_vf2_callback_data_t*)extra;
  igraph_bool_t retval;
  PyObject *result;

  result = PyObject_CallFunction(data->edge_compat_fn, "OOll",
      data->graph1, data->graph2, (long)cand1, (long)cand2);

  if (result == NULL) {
    /* Error in callback, return 0 */
    PyErr_WriteUnraisable(data->edge_compat_fn);
    return 0;
  }

  retval = PyObject_IsTrue(result);
  Py_DECREF(result);

  return retval;
}

/** \ingroup python_interface_graph
 * \brief Determines whether the graph is isomorphic to another graph,
 *   using the VF2 isomorphism algorithm
 *
 * The actual code is almost the same as igraphmodule_Graph_subisomorphic.
 * Be sure to correct bugs in both interfaces if applicable!
 *
 * \sa igraph_isomorphic_vf2
 */
PyObject *igraphmodule_Graph_isomorphic_vf2(igraphmodule_GraphObject * self,
                                            PyObject * args, PyObject * kwds)
{
  igraph_bool_t result = 0;
  PyObject *o=Py_None, *return1=Py_False, *return2=Py_False;
  PyObject *color1_o=Py_None, *color2_o=Py_None;
  PyObject *edge_color1_o=Py_None, *edge_color2_o=Py_None;
  PyObject *callback_fn=Py_None;
  PyObject *node_compat_fn=Py_None, *edge_compat_fn=Py_None;
  igraphmodule_GraphObject *other;
  igraph_vector_t mapping_12, mapping_21;
  igraph_vector_t *map12=0, *map21=0;
  igraph_vector_int_t *color1=0, *color2=0;
  igraph_vector_int_t *edge_color1=0, *edge_color2=0;
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t callback_data;
  int retval;

  static char *kwlist[] = {
    "other", "color1", "color2", "edge_color1", "edge_color2",
    "return_mapping_12", "return_mapping_21", "callback",
    "node_compat_fn", "edge_compat_fn", NULL
  };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "|O!OOOOOOOOO", kwlist, &igraphmodule_GraphType, &o,
       &color1_o, &color2_o, &edge_color1, &edge_color2, &return1, &return2,
       &callback_fn, &node_compat_fn, &edge_compat_fn))
    return NULL;

  if (o == Py_None)
    other=self;
  else
    other=(igraphmodule_GraphObject*)o;

  if (callback_fn != Py_None && !PyCallable_Check(callback_fn)) {
    PyErr_SetString(PyExc_TypeError, "callback must be None or callable");
    return NULL;
  }

  if (node_compat_fn != Py_None && !PyCallable_Check(node_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "node_compat_fn must be None or callable");
    return NULL;
  }

  if (edge_compat_fn != Py_None && !PyCallable_Check(edge_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "edge_compat_fn must be None or callable");
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_int_t(color1_o, self, &color1,
	  ATTRIBUTE_TYPE_VERTEX)) return NULL;
  if (igraphmodule_attrib_to_vector_int_t(color2_o, other, &color2,
	  ATTRIBUTE_TYPE_VERTEX)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color1_o, self, &edge_color1,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color2_o, other, &edge_color2,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    return NULL;
  }

  if (PyObject_IsTrue(return1)) {
	igraph_vector_init(&mapping_12, 0);
	map12 = &mapping_12;
  }
  if (PyObject_IsTrue(return2)) {
	igraph_vector_init(&mapping_21, 0);
	map21 = &mapping_21;
  }

  callback_data.graph1 = (PyObject*)self;
  callback_data.graph2 = (PyObject*)other;
  callback_data.callback_fn = callback_fn == Py_None ? 0 : callback_fn;
  callback_data.node_compat_fn = node_compat_fn == Py_None ? 0 : node_compat_fn;
  callback_data.edge_compat_fn = edge_compat_fn == Py_None ? 0 : edge_compat_fn;

  if (callback_data.callback_fn == 0) {
    retval = igraph_isomorphic_vf2(&self->g, &other->g,
        color1, color2, edge_color1, edge_color2, &result, map12, map21,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data);
  } else {
    retval = igraph_isomorphic_function_vf2(&self->g, &other->g,
        color1, color2, edge_color1, edge_color2, map12, map21,
        igraphmodule_i_Graph_isomorphic_vf2_callback_fn,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data);
  }

  if (color1) { igraph_vector_int_destroy(color1); free(color1); }
  if (color2) { igraph_vector_int_destroy(color2); free(color2); }
  if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
  if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }

  if (retval) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (!map12 && !map21) {
    if (result) Py_RETURN_TRUE;
    Py_RETURN_FALSE;
  } else {
	PyObject *m1, *m2;
	if (map12) {
	  m1 = igraphmodule_vector_t_to_PyList(map12, IGRAPHMODULE_TYPE_INT);
	  igraph_vector_destroy(map12);
	  if (!m1) {
		if (map21) igraph_vector_destroy(map21);
		return NULL;
	  }
	} else { m1 = Py_None; Py_INCREF(m1); }
	if (map21) {
	  m2 = igraphmodule_vector_t_to_PyList(map21, IGRAPHMODULE_TYPE_INT);
	  igraph_vector_destroy(map21);
	  if (!m2) {
	    Py_DECREF(m1);
		return NULL;
	  }
	} else { m2 = Py_None; Py_INCREF(m2); }
	return Py_BuildValue("ONN", result ? Py_True : Py_False, m1, m2);
  }
}

/** \ingroup python_interface_graph
 * \brief Counts the number of isomorphisms of two given graphs 
 *
 * The actual code is almost the same as igraphmodule_Graph_count_subisomorphisms.
 * Make sure you correct bugs in both interfaces if applicable!
 *
 * \sa igraph_count_isomorphisms_vf2
 */
PyObject *igraphmodule_Graph_count_isomorphisms_vf2(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  igraph_integer_t result = 0;
  PyObject *o = Py_None;
  PyObject *color1_o=Py_None, *color2_o=Py_None;
  PyObject *edge_color1_o=Py_None, *edge_color2_o=Py_None;
  PyObject *node_compat_fn=Py_None, *edge_compat_fn=Py_None;
  igraphmodule_GraphObject *other;
  igraph_vector_int_t *color1=0, *color2=0;
  igraph_vector_int_t *edge_color1=0, *edge_color2=0;
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t callback_data;

  static char *kwlist[] = { "other", "color1", "color2", "edge_color1", "edge_color2",
    "node_compat_fn", "edge_compat_fn", NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "|O!OOOOOO", kwlist, &igraphmodule_GraphType, &o,
       &color1_o, &color2_o, &edge_color1_o, &edge_color2_o,
       &node_compat_fn, &edge_compat_fn))
    return NULL;

  if (o == Py_None)
    other=self;
  else
    other=(igraphmodule_GraphObject*)o;

  if (node_compat_fn != Py_None && !PyCallable_Check(node_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "node_compat_fn must be None or callable");
    return NULL;
  }

  if (edge_compat_fn != Py_None && !PyCallable_Check(edge_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "edge_compat_fn must be None or callable");
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_int_t(color1_o, self, &color1,
	  ATTRIBUTE_TYPE_VERTEX)) return NULL;
  if (igraphmodule_attrib_to_vector_int_t(color2_o, other, &color2,
	  ATTRIBUTE_TYPE_VERTEX)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color1_o, self, &edge_color1,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color2_o, other, &edge_color2,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    return NULL;
  }

  callback_data.graph1 = (PyObject*)self;
  callback_data.graph2 = (PyObject*)other;
  callback_data.callback_fn = 0;
  callback_data.node_compat_fn = node_compat_fn == Py_None ? 0 : node_compat_fn;
  callback_data.edge_compat_fn = edge_compat_fn == Py_None ? 0 : edge_compat_fn;

  if (igraph_count_isomorphisms_vf2(&self->g, &other->g,
        color1, color2, edge_color1, edge_color2, &result,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (color1) { igraph_vector_int_destroy(color1); free(color1); }
  if (color2) { igraph_vector_int_destroy(color2); free(color2); }
  if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
  if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }

  return Py_BuildValue("l", (long)result);
}

/** \ingroup python_interface_graph
 * \brief Returns all isomorphisms of two given graphs 
 *
 * The actual code is almost the same as igraphmodule_Graph_get_subisomorphisms.
 * Make sure you correct bugs in both interfaces if applicable!
 *
 * \sa igraph_get_isomorphisms_vf2
 */
PyObject *igraphmodule_Graph_get_isomorphisms_vf2(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  igraph_vector_ptr_t result;
  PyObject *o = Py_None;
  PyObject *color1_o = Py_None, *color2_o = Py_None;
  PyObject *edge_color1_o=Py_None, *edge_color2_o=Py_None;
  PyObject *node_compat_fn=Py_None, *edge_compat_fn=Py_None;
  PyObject *res;
  igraphmodule_GraphObject *other;
  igraph_vector_int_t *color1=0, *color2=0;
  igraph_vector_int_t *edge_color1=0, *edge_color2=0;
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t callback_data;

  static char *kwlist[] = { "other", "color1", "color2",
    "edge_color1", "edge_color2", "node_compat_fn", "edge_compat_fn", NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "|O!OOOOOO", kwlist, &igraphmodule_GraphType, &o,
         &color1_o, &color2_o, &edge_color1_o, &edge_color2_o,
         &node_compat_fn, &edge_compat_fn))
    return NULL;

  if (o == Py_None)
    other=self;
  else
    other=(igraphmodule_GraphObject*)o;

  if (node_compat_fn != Py_None && !PyCallable_Check(node_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "node_compat_fn must be None or callable");
    return NULL;
  }

  if (edge_compat_fn != Py_None && !PyCallable_Check(edge_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "edge_compat_fn must be None or callable");
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_int_t(color1_o, self, &color1,
	  ATTRIBUTE_TYPE_VERTEX)) return NULL;
  if (igraphmodule_attrib_to_vector_int_t(color2_o, other, &color2,
	  ATTRIBUTE_TYPE_VERTEX)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color1_o, self, &edge_color1,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color2_o, other, &edge_color2,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    return NULL;
  }

  if (igraph_vector_ptr_init(&result, 0)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }
    return igraphmodule_handle_igraph_error();
  }

  callback_data.graph1 = (PyObject*)self;
  callback_data.graph2 = (PyObject*)other;
  callback_data.callback_fn = 0;
  callback_data.node_compat_fn = node_compat_fn == Py_None ? 0 : node_compat_fn;
  callback_data.edge_compat_fn = edge_compat_fn == Py_None ? 0 : edge_compat_fn;

  if (igraph_get_isomorphisms_vf2(&self->g, &other->g,
        color1, color2, edge_color1, edge_color2, &result,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data)) {
    igraphmodule_handle_igraph_error();
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }
    igraph_vector_ptr_destroy(&result);
    return NULL;
  }

  if (color1) { igraph_vector_int_destroy(color1); free(color1); }
  if (color2) { igraph_vector_int_destroy(color2); free(color2); }
  if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
  if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }

  res = igraphmodule_vector_ptr_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);

  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&result, igraph_vector_destroy);
  igraph_vector_ptr_destroy_all(&result);

  return res;
}

/** \ingroup python_interface_graph
 * \brief Determines whether a subgraph of the graph is isomorphic to another graph
 *
 * The actual code is almost the same as igraphmodule_Graph_isomorphic. Make sure
 * you correct bugs in both interfaces if applicable!
 *
 * \sa igraph_subisomorphic_vf2
 */
PyObject *igraphmodule_Graph_subisomorphic_vf2(igraphmodule_GraphObject * self,
                                        PyObject * args, PyObject * kwds)
{
  igraph_bool_t result = 0;
  PyObject *o, *return1=Py_False, *return2=Py_False;
  PyObject *color1_o=Py_None, *color2_o=Py_None;
  PyObject *edge_color1_o=Py_None, *edge_color2_o=Py_None;
  PyObject *callback_fn=Py_None;
  PyObject *node_compat_fn=Py_None, *edge_compat_fn=Py_None;
  igraphmodule_GraphObject *other;
  igraph_vector_t mapping_12, mapping_21, *map12=0, *map21=0;
  igraph_vector_int_t *color1=0, *color2=0;
  igraph_vector_int_t *edge_color1=0, *edge_color2=0;
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t callback_data;
  int retval;

  static char *kwlist[] = { "other", "color1", "color2", "edge_color1", "edge_color2",
    "return_mapping_12", "return_mapping_21",
    "callback", "node_compat_fn", "edge_compat_fn",
    NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "O!|OOOOOOOOO", kwlist, &igraphmodule_GraphType, &o,
       &color1_o, &color2_o, &edge_color1_o, &edge_color2_o, &return1, &return2,
       &callback_fn, &node_compat_fn, &edge_compat_fn))
    return NULL;

  other=(igraphmodule_GraphObject*)o;

  if (callback_fn != Py_None && !PyCallable_Check(callback_fn)) {
    PyErr_SetString(PyExc_TypeError, "callback must be None or callable");
    return NULL;
  }

  if (node_compat_fn != Py_None && !PyCallable_Check(node_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "node_compat_fn must be None or callable");
    return NULL;
  }

  if (edge_compat_fn != Py_None && !PyCallable_Check(edge_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "edge_compat_fn must be None or callable");
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_int_t(color1_o, self, &color1,
	  ATTRIBUTE_TYPE_VERTEX)) return NULL;
  if (igraphmodule_attrib_to_vector_int_t(color2_o, other, &color2,
	  ATTRIBUTE_TYPE_VERTEX)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color1_o, self, &edge_color1,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color2_o, other, &edge_color2,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    return NULL;
  }

  if (PyObject_IsTrue(return1)) {
	igraph_vector_init(&mapping_12, 0);
	map12 = &mapping_12;
  }
  if (PyObject_IsTrue(return2)) {
	igraph_vector_init(&mapping_21, 0);
	map21 = &mapping_21;
  }

  callback_data.graph1 = (PyObject*)self;
  callback_data.graph2 = (PyObject*)other;
  callback_data.callback_fn = callback_fn == Py_None ? 0 : callback_fn;
  callback_data.node_compat_fn = node_compat_fn == Py_None ? 0 : node_compat_fn;
  callback_data.edge_compat_fn = edge_compat_fn == Py_None ? 0 : edge_compat_fn;

  if (callback_data.callback_fn == 0) {
    retval = igraph_subisomorphic_vf2(&self->g, &other->g,
        color1, color2, edge_color1, edge_color2, &result, map12, map21,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data);
  } else {
    retval = igraph_subisomorphic_function_vf2(&self->g, &other->g,
        color1, color2, edge_color1, edge_color2, map12, map21,
        igraphmodule_i_Graph_isomorphic_vf2_callback_fn,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data);
  }

  if (color1) { igraph_vector_int_destroy(color1); free(color1); }
  if (color2) { igraph_vector_int_destroy(color2); free(color2); }
  if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
  if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }

  if (retval) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (!map12 && !map21) {
    if (result)
      Py_RETURN_TRUE;
    Py_RETURN_FALSE;
  } else {
	PyObject *m1, *m2;
	if (map12) {
	  m1 = igraphmodule_vector_t_to_PyList(map12, IGRAPHMODULE_TYPE_INT);
	  igraph_vector_destroy(map12);
	  if (!m1) {
		if (map21) igraph_vector_destroy(map21);
		return NULL;
	  }
	} else {
      m1 = Py_None; Py_INCREF(m1);
    }
	if (map21) {
	  m2 = igraphmodule_vector_t_to_PyList(map21, IGRAPHMODULE_TYPE_INT);
	  igraph_vector_destroy(map21);
	  if (!m2) {
	    Py_DECREF(m1);
		return NULL;
	  }
	} else {
      m2 = Py_None; Py_INCREF(m2);
    }
	return Py_BuildValue("ONN", result ? Py_True : Py_False, m1, m2);
  }
}

/** \ingroup python_interface_graph
 * \brief Counts the number of subisomorphisms of two given graphs 
 *
 * The actual code is almost the same as igraphmodule_Graph_count_isomorphisms.
 * Make sure you correct bugs in both interfaces if applicable!
 *
 * \sa igraph_count_subisomorphisms_vf2
 */
PyObject *igraphmodule_Graph_count_subisomorphisms_vf2(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  igraph_integer_t result = 0;
  PyObject *o = Py_None;
  PyObject *color1_o = Py_None, *color2_o = Py_None;
  PyObject *edge_color1_o=Py_None, *edge_color2_o=Py_None;
  PyObject *node_compat_fn=Py_None, *edge_compat_fn=Py_None;
  igraph_vector_int_t *color1=0, *color2=0;
  igraph_vector_int_t *edge_color1=0, *edge_color2=0;
  igraphmodule_GraphObject *other;
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t callback_data;

  static char *kwlist[] = { "other", "color1", "color2", "edge_color1",
    "edge_color2", "node_compat_fn", "edge_compat_fn", NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "O!|OOOOOO", kwlist, &igraphmodule_GraphType, &o,
         &color1_o, &color2_o, &edge_color1_o, &edge_color2_o,
         &node_compat_fn, &edge_compat_fn))
    return NULL;

  other=(igraphmodule_GraphObject*)o;

  if (node_compat_fn != Py_None && !PyCallable_Check(node_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "node_compat_fn must be None or callable");
    return NULL;
  }

  if (edge_compat_fn != Py_None && !PyCallable_Check(edge_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "edge_compat_fn must be None or callable");
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_int_t(color1_o, self, &color1,
	  ATTRIBUTE_TYPE_VERTEX)) return NULL;
  if (igraphmodule_attrib_to_vector_int_t(color2_o, other, &color2,
	  ATTRIBUTE_TYPE_VERTEX)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color1_o, self, &edge_color1,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color2_o, other, &edge_color2,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    return NULL;
  }

  callback_data.graph1 = (PyObject*)self;
  callback_data.graph2 = (PyObject*)other;
  callback_data.callback_fn = 0;
  callback_data.node_compat_fn = node_compat_fn == Py_None ? 0 : node_compat_fn;
  callback_data.edge_compat_fn = edge_compat_fn == Py_None ? 0 : edge_compat_fn;

  if (igraph_count_subisomorphisms_vf2(&self->g, &other->g, color1, color2,
        edge_color1, edge_color2, &result,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data)) {
    igraphmodule_handle_igraph_error();
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }
    return NULL;
  }

  if (color1) { igraph_vector_int_destroy(color1); free(color1); }
  if (color2) { igraph_vector_int_destroy(color2); free(color2); }
  if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
  if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }

  return Py_BuildValue("l", (long)result);
}

/** \ingroup python_interface_graph
 * \brief Returns all subisomorphisms of two given graphs 
 *
 * The actual code is almost the same as igraphmodule_Graph_get_isomorphisms.
 * Make sure you correct bugs in both interfaces if applicable!
 *
 * \sa igraph_get_isomorphisms_vf2
 */
PyObject *igraphmodule_Graph_get_subisomorphisms_vf2(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  igraph_vector_ptr_t result;
  PyObject *o;
  PyObject *color1_o=Py_None, *color2_o=Py_None;
  PyObject *edge_color1_o=Py_None, *edge_color2_o=Py_None;
  PyObject *node_compat_fn=Py_None, *edge_compat_fn=Py_None;
  PyObject *res;
  igraphmodule_GraphObject *other;
  igraph_vector_int_t *color1=0, *color2=0;
  igraph_vector_int_t *edge_color1=0, *edge_color2=0;
  igraphmodule_i_Graph_isomorphic_vf2_callback_data_t callback_data;

  static char *kwlist[] = { "other", "color1", "color2", "edge_color1",
    "edge_color2", "node_compat_fn", "edge_compat_fn", NULL };

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "O!|OOOOOO", kwlist, &igraphmodule_GraphType, &o,
       &color1_o, &color2_o, &edge_color1_o, &edge_color2_o,
       &node_compat_fn, &edge_compat_fn))
    return NULL;

  if (igraph_vector_ptr_init(&result, 0)) {
    return igraphmodule_handle_igraph_error();
  }

  other=(igraphmodule_GraphObject*)o;

  if (node_compat_fn != Py_None && !PyCallable_Check(node_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "node_compat_fn must be None or callable");
    return NULL;
  }

  if (edge_compat_fn != Py_None && !PyCallable_Check(edge_compat_fn)) {
    PyErr_SetString(PyExc_TypeError, "edge_compat_fn must be None or callable");
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_int_t(color1_o, self, &color1,
	  ATTRIBUTE_TYPE_VERTEX)) return NULL;
  if (igraphmodule_attrib_to_vector_int_t(color2_o, other, &color2,
	  ATTRIBUTE_TYPE_VERTEX)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color1_o, self, &edge_color1,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    return NULL;
  }
  if (igraphmodule_attrib_to_vector_int_t(edge_color2_o, other, &edge_color2,
	  ATTRIBUTE_TYPE_EDGE)) {
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    return NULL;
  }

  callback_data.graph1 = (PyObject*)self;
  callback_data.graph2 = (PyObject*)other;
  callback_data.callback_fn = 0;
  callback_data.node_compat_fn = node_compat_fn == Py_None ? 0 : node_compat_fn;
  callback_data.edge_compat_fn = edge_compat_fn == Py_None ? 0 : edge_compat_fn;

  if (igraph_get_subisomorphisms_vf2(&self->g, &other->g, color1, color2,
        edge_color1, edge_color2, &result,
        node_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_node_compat_fn,
        edge_compat_fn == Py_None ? 0 : igraphmodule_i_Graph_isomorphic_vf2_edge_compat_fn,
        &callback_data)) {
    igraphmodule_handle_igraph_error();
    if (color1) { igraph_vector_int_destroy(color1); free(color1); }
    if (color2) { igraph_vector_int_destroy(color2); free(color2); }
    if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
    if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }
    igraph_vector_ptr_destroy(&result);
    return NULL;
  }

  if (color1) { igraph_vector_int_destroy(color1); free(color1); }
  if (color2) { igraph_vector_int_destroy(color2); free(color2); }
  if (edge_color1) { igraph_vector_int_destroy(edge_color1); free(edge_color1); }
  if (edge_color2) { igraph_vector_int_destroy(edge_color2); free(edge_color2); }

  res = igraphmodule_vector_ptr_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);

  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&result, igraph_vector_destroy);
  igraph_vector_ptr_destroy_all(&result);

  return res;
}

/**********************************************************************
 * Graph attribute handling                                           *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Returns the number of graph attributes
 */
Py_ssize_t igraphmodule_Graph_attribute_count(igraphmodule_GraphObject * self)
{
  return PyDict_Size(((PyObject **) self->g.attr)[ATTRHASH_IDX_GRAPH]);
}

/** \ingroup python_interface_graph
 * \brief Handles the subscript operator on the graph.
 * 
 * When the subscript is a string, returns the corresponding value of the
 * given attribute in the graph. When the subscript is a tuple of length
 * 2, retrieves the adjacency matrix representation of the graph between
 * some vertices.
 */
PyObject *igraphmodule_Graph_mp_subscript(igraphmodule_GraphObject * self,
                                          PyObject * s)
{
  PyObject *result = 0;

  if (PyTuple_Check(s) && PyTuple_Size(s) >= 2) {
    /* Adjacency matrix representation */
    PyObject *ri = PyTuple_GET_ITEM(s, 0);
    PyObject *ci = PyTuple_GET_ITEM(s, 1);
    PyObject *attr;
    
    if (PyTuple_Size(s) == 2) {
      attr = 0;
    } else if (PyTuple_Size(s) == 3) {
      attr = PyTuple_GET_ITEM(s, 2);
    } else {
      PyErr_SetString(PyExc_TypeError, "adjacency matrix indexing must use at most three arguments");
      return 0;
    }

    return igraphmodule_Graph_adjmatrix_get_index(&self->g, ri, ci, attr);
  }

  /* Ordinary attribute retrieval */
  result = PyDict_GetItem(ATTR_STRUCT_DICT(&self->g)[ATTRHASH_IDX_GRAPH], s);
  if (result) {
    Py_INCREF(result);
    return result;
  }

  /* result is NULL, check whether there was an error */
  if (!PyErr_Occurred())
    PyErr_SetString(PyExc_KeyError, "Attribute does not exist");

  return NULL;
}

/** \ingroup python_interface_graph
 * \brief Handles the subscript assignment operator on the graph.
 * 
 * If k is a string, sets the value of the corresponding attribute of the graph.
 * If k is a tuple of length 2, sets part of the adjacency matrix.
 *
 * \return 0 if everything's ok, -1 in case of error
 */
int igraphmodule_Graph_mp_assign_subscript(igraphmodule_GraphObject * self,
                                     PyObject * k, PyObject * v)
{
  PyObject* dict = ATTR_STRUCT_DICT(&self->g)[ATTRHASH_IDX_GRAPH];

  if (PyTuple_Check(k) && PyTuple_Size(k) >= 2) {
    /* Adjacency matrix representation */
    PyObject *ri, *ci, *attr;
    
    if (v == NULL) {
      PyErr_SetString(PyExc_NotImplementedError, "cannot delete parts "
          "of the adjacency matrix of a graph");
      return -1;
    }

    ri = PyTuple_GET_ITEM(k, 0);
    ci = PyTuple_GET_ITEM(k, 1);

    if (PyTuple_Size(k) == 2) {
      attr = 0;
    } else if (PyTuple_Size(k) == 3) {
      attr = PyTuple_GET_ITEM(k, 2);
    } else {
      PyErr_SetString(PyExc_TypeError, "adjacency matrix indexing must use at most three arguments");
      return 0;
    }

    return igraphmodule_Graph_adjmatrix_set_index(&self->g, ri, ci, attr, v);
  }

  /* Ordinary attribute setting/deletion */
  if (v == NULL)
    return PyDict_DelItem(dict, k);

  if (PyDict_SetItem(dict, k, v) == -1)
    return -1;
  return 0;
}

/** \ingroup python_interface_graph
 * \brief Returns the attribute list of the graph
 */
PyObject *igraphmodule_Graph_attributes(igraphmodule_GraphObject * self)
{
  return PyDict_Keys(ATTR_STRUCT_DICT(&self->g)[ATTRHASH_IDX_GRAPH]);
}

/** \ingroup python_interface_graph
 * \brief Returns the attribute list of the graph's vertices
 */
PyObject *igraphmodule_Graph_vertex_attributes(igraphmodule_GraphObject *
                                               self)
{
  return PyDict_Keys(ATTR_STRUCT_DICT(&self->g)[ATTRHASH_IDX_VERTEX]);
}

/** \ingroup python_interface_graph
 * \brief Returns the attribute list of the graph's edges
 */
PyObject *igraphmodule_Graph_edge_attributes(igraphmodule_GraphObject * self)
{
  return PyDict_Keys(ATTR_STRUCT_DICT(&self->g)[ATTRHASH_IDX_EDGE]);
}

/**********************************************************************
 * Graph operations (union, intersection etc)                         *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Creates the disjoint union of two graphs (operator version)
 */
PyObject *igraphmodule_Graph_disjoint_union(igraphmodule_GraphObject * self,
                                            PyObject * other)
{
  PyObject *it;
  igraphmodule_GraphObject *o, *result;
  igraph_t g;

  /* Did we receive an iterable? */
  it = PyObject_GetIter(other);
  if (it) {
    /* Get all elements, store the graphs in an igraph_vector_ptr */
    igraph_vector_ptr_t gs;
    if (igraph_vector_ptr_init(&gs, 0)) {
      Py_DECREF(it);
      return igraphmodule_handle_igraph_error();
    }
    if (igraph_vector_ptr_push_back(&gs, &self->g)) {
      Py_DECREF(it);
      igraph_vector_ptr_destroy(&gs);
      return igraphmodule_handle_igraph_error();
    }
    if (igraphmodule_append_PyIter_to_vector_ptr_t(it, &gs)) {
      igraph_vector_ptr_destroy(&gs);
      Py_DECREF(it);
      return NULL;
    }
    Py_DECREF(it);

    /* Create disjoint union */
    if (igraph_disjoint_union_many(&g, &gs)) {
      igraph_vector_ptr_destroy(&gs);
      return igraphmodule_handle_igraph_error();
    }

    igraph_vector_ptr_destroy(&gs);
  }
  else {
    PyErr_Clear();
    if (!PyObject_TypeCheck(other, &igraphmodule_GraphType)) {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    }
    o = (igraphmodule_GraphObject *) other;

    if (igraph_disjoint_union(&g, &self->g, &o->g)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  /* this is correct as long as attributes are not copied by the
   * operator. if they are copied, the initialization should not empty
   * the attribute hashes */
  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Creates the union of two graphs (operator version)
 */
PyObject *igraphmodule_Graph_union(igraphmodule_GraphObject * self,
                                   PyObject * other)
{
  PyObject *it;
  igraphmodule_GraphObject *o, *result;
  igraph_t g;

  /* Did we receive an iterable? */
  it = PyObject_GetIter(other);
  if (it) {
    /* Get all elements, store the graphs in an igraph_vector_ptr */
    igraph_vector_ptr_t gs;
    if (igraph_vector_ptr_init(&gs, 0)) {
      Py_DECREF(it);
      return igraphmodule_handle_igraph_error();
    }
    if (igraph_vector_ptr_push_back(&gs, &self->g)) {
      Py_DECREF(it);
      igraph_vector_ptr_destroy(&gs);
      return igraphmodule_handle_igraph_error();
    }
    if (igraphmodule_append_PyIter_to_vector_ptr_t(it, &gs)) {
      Py_DECREF(it);
      igraph_vector_ptr_destroy(&gs);
      return NULL;
    }
    Py_DECREF(it);

    /* Create union */
    if (igraph_union_many(&g, &gs)) {
      igraph_vector_ptr_destroy(&gs);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    igraph_vector_ptr_destroy(&gs);
  }
  else {
    PyErr_Clear();
    if (!PyObject_TypeCheck(other, &igraphmodule_GraphType)) {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    }
    o = (igraphmodule_GraphObject *) other;

    if (igraph_union(&g, &self->g, &o->g)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  /* this is correct as long as attributes are not copied by the
   * operator. if they are copied, the initialization should not empty
   * the attribute hashes */
  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Creates the intersection of two graphs (operator version)
 */
PyObject *igraphmodule_Graph_intersection(igraphmodule_GraphObject * self,
                                          PyObject * other)
{
  PyObject *it;
  igraphmodule_GraphObject *o, *result;
  igraph_t g;

  /* Did we receive an iterable? */
  it = PyObject_GetIter(other);
  if (it) {
    /* Get all elements, store the graphs in an igraph_vector_ptr */
    igraph_vector_ptr_t gs;
    if (igraph_vector_ptr_init(&gs, 0)) {
      Py_DECREF(it);
      return igraphmodule_handle_igraph_error();
    }
    if (igraph_vector_ptr_push_back(&gs, &self->g)) {
      Py_DECREF(it);
      igraph_vector_ptr_destroy(&gs);
      return igraphmodule_handle_igraph_error();
    }
    if (igraphmodule_append_PyIter_to_vector_ptr_t(it, &gs)) {
      Py_DECREF(it);
      igraph_vector_ptr_destroy(&gs);
      return NULL;
    }
    Py_DECREF(it);

    /* Create union */
    if (igraph_intersection_many(&g, &gs)) {
      igraph_vector_ptr_destroy(&gs);
      igraphmodule_handle_igraph_error();
      return NULL;
    }

    igraph_vector_ptr_destroy(&gs);
  }
  else {
    PyErr_Clear();
    if (!PyObject_TypeCheck(other, &igraphmodule_GraphType)) {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    }
    o = (igraphmodule_GraphObject *) other;

    if (igraph_intersection(&g, &self->g, &o->g)) {
      igraphmodule_handle_igraph_error();
      return NULL;
    }
  }

  /* this is correct as long as attributes are not copied by the
   * operator. if they are copied, the initialization should not empty
   * the attribute hashes */
  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Creates the difference of two graphs (operator version)
 */
PyObject *igraphmodule_Graph_difference(igraphmodule_GraphObject * self,
                                        PyObject * other)
{
  igraphmodule_GraphObject *o, *result;
  igraph_t g;

  if (!PyObject_TypeCheck(other, &igraphmodule_GraphType)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  o = (igraphmodule_GraphObject *) other;

  if (igraph_difference(&g, &self->g, &o->g)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  /* this is correct as long as attributes are not copied by the
   * operator. if they are copied, the initialization should not empty
   * the attribute hashes */
  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Creates the complementer of a graph
 */
PyObject *igraphmodule_Graph_complementer(igraphmodule_GraphObject * self,
                                          PyObject * args)
{
  igraphmodule_GraphObject *result;
  PyObject *o = Py_True;
  igraph_t g;

  if (!PyArg_ParseTuple(args, "|O", &o))
    return NULL;
  if (igraph_complementer(&g, &self->g, PyObject_IsTrue(o))) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  /* this is correct as long as attributes are not copied by the
   * operator. if they are copied, the initialization should not empty
   * the attribute hashes */
  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Creates the complementer of a graph (operator version)
 */
PyObject *igraphmodule_Graph_complementer_op(igraphmodule_GraphObject * self)
{
  igraphmodule_GraphObject *result;
  igraph_t g;

  if (igraph_complementer(&g, &self->g, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  /* this is correct as long as attributes are not copied by the
   * operator. if they are copied, the initialization should not empty
   * the attribute hashes */
  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/** \ingroup python_interface_graph
 * \brief Creates the composition of two graphs
 */
PyObject *igraphmodule_Graph_compose(igraphmodule_GraphObject * self,
                                     PyObject * other)
{
  igraphmodule_GraphObject *o, *result;
  igraph_t g;

  if (!PyObject_TypeCheck(other, &igraphmodule_GraphType)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  o = (igraphmodule_GraphObject *) other;

  if (igraph_compose(&g, &self->g, &o->g)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  /* this is correct as long as attributes are not copied by the
   * operator. if they are copied, the initialization should not empty
   * the attribute hashes */
  CREATE_GRAPH(result, g);

  return (PyObject *) result;
}

/**********************************************************************
 * Graph traversal algorithms                                         *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Conducts a breadth first search (BFS) on the graph
 */
PyObject *igraphmodule_Graph_bfs(igraphmodule_GraphObject * self,
                                 PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "vid", "mode", NULL };
  long vid;
  PyObject *l1, *l2, *l3, *result, *mode_o=Py_None;
  igraph_neimode_t mode = IGRAPH_OUT;
  igraph_vector_t vids;
  igraph_vector_t layers;
  igraph_vector_t parents;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|O", kwlist, &vid, &mode_o))
    return NULL;
  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;

  if (igraph_vector_init(&vids, igraph_vcount(&self->g))) 
    return igraphmodule_handle_igraph_error();
  if (igraph_vector_init(&layers, igraph_vcount(&self->g))) {
    igraph_vector_destroy(&vids);
    return igraphmodule_handle_igraph_error();
  }
  if (igraph_vector_init(&parents, igraph_vcount(&self->g))) {
    igraph_vector_destroy(&vids); igraph_vector_destroy(&parents);
    return igraphmodule_handle_igraph_error();
  }
  if (igraph_i_bfs
      (&self->g, (igraph_integer_t) vid, mode, &vids, &layers, &parents)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  l1 = igraphmodule_vector_t_to_PyList(&vids, IGRAPHMODULE_TYPE_INT);
  l2 = igraphmodule_vector_t_to_PyList(&layers, IGRAPHMODULE_TYPE_INT);
  l3 = igraphmodule_vector_t_to_PyList(&parents, IGRAPHMODULE_TYPE_INT);
  if (l1 && l2 && l3) {
    result = Py_BuildValue("NNN", l1, l2, l3);    /* references stolen */
  } else {
    if (l1) { Py_DECREF(l1); }
    if (l2) { Py_DECREF(l2); }
    if (l3) { Py_DECREF(l3); }
    result = NULL;
  }
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&layers);
  igraph_vector_destroy(&parents);
  return result;
}

/** \ingroup python_interface_graph
 * \brief Constructs a breadth first search (BFS) iterator of the graph
 */
PyObject *igraphmodule_Graph_bfsiter(igraphmodule_GraphObject * self,
                                     PyObject * args, PyObject * kwds)
{
  char *kwlist[] = { "vid", "mode", "advanced", NULL };
  PyObject *root, *adv = Py_False, *mode_o = Py_None;
  igraph_neimode_t mode = IGRAPH_OUT;

  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "O|OO", kwlist, &root, &mode_o, &adv))
    return NULL;
  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;
  return igraphmodule_BFSIter_new(self, root, mode, PyObject_IsTrue(adv));
}

/** \ingroup python_interface_graph
 * \brief Unfolds a graph into a tree using BFS
 */
PyObject *igraphmodule_Graph_unfold_tree(igraphmodule_GraphObject * self,
                                 PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "roots", "mode", NULL };
  igraphmodule_GraphObject *result_o;
  PyObject *mapping_o, *mode_o=Py_None, *roots_o=Py_None;
  igraph_neimode_t mode = IGRAPH_OUT;
  igraph_vs_t vs;
  igraph_vector_t mapping, vids;
  igraph_t result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &roots_o, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;
  if (igraphmodule_PyObject_to_vs_t(roots_o, &vs, &self->g, 0, 0)) return NULL;

  if (igraph_vector_init(&mapping, igraph_vcount(&self->g))) {
    igraph_vs_destroy(&vs);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_vector_init(&vids, 0)) {
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&mapping);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_vs_as_vector(&self->g, vs, &vids)) {
    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&vids);
    igraph_vector_destroy(&mapping);
    return igraphmodule_handle_igraph_error();
  }

  igraph_vs_destroy(&vs);

  if (igraph_unfold_tree(&self->g, &result, mode, &vids, &mapping)) {
    igraph_vector_destroy(&vids);
    igraph_vector_destroy(&mapping);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  igraph_vector_destroy(&vids);

  mapping_o = igraphmodule_vector_t_to_PyList(&mapping, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&mapping);

  if (!mapping_o) {
    igraph_destroy(&result);
    return NULL;
  }

  CREATE_GRAPH(result_o, result);

  return Py_BuildValue("NN", result_o, mapping_o);
}

/**********************************************************************
 * Maximum flows                                                      *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Calculates the value of the maximum flow in the graph
 */
PyObject *igraphmodule_Graph_maxflow_value(igraphmodule_GraphObject * self,
                                           PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "source", "target", "capacity", NULL };
  PyObject *capacity_object = Py_None;
  igraph_vector_t capacity_vector;
  igraph_real_t result;
  long int vid1 = -1, vid2 = -1;
  igraph_integer_t v1, v2;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|O", kwlist,
                                   &vid1, &vid2, &capacity_object))
    return NULL;

  v1 = (igraph_integer_t) vid1;
  v2 = (igraph_integer_t) vid2;
  if (igraphmodule_PyObject_to_attribute_values(capacity_object,
                                                &capacity_vector,
                                                self, ATTRHASH_IDX_EDGE, 1.0))
    return igraphmodule_handle_igraph_error();


  if (igraph_maxflow_value(&self->g, &result, v1, v2, &capacity_vector)) {
    igraph_vector_destroy(&capacity_vector);
    return igraphmodule_handle_igraph_error();
  }

  igraph_vector_destroy(&capacity_vector);
  return Py_BuildValue("d", (double)result);
}

/** \ingroup python_interface_graph
 * \brief Calculates the maximum flow of the graph
 */
PyObject *igraphmodule_Graph_maxflow(igraphmodule_GraphObject * self,
                                     PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "source", "target", "capacity", NULL };
  PyObject *capacity_object = Py_None, *flow_o, *cut_o, *partition_o;
  igraph_vector_t capacity_vector;
  igraph_real_t result;
  long int vid1 = -1, vid2 = -1;
  igraph_integer_t v1, v2;
  igraph_vector_t flow, cut, partition;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|O", kwlist,
                                   &vid1, &vid2, &capacity_object))
    return NULL;

  v1 = (igraph_integer_t) vid1;
  v2 = (igraph_integer_t) vid2;
  if (igraphmodule_PyObject_to_attribute_values(capacity_object,
                                                &capacity_vector,
                                                self, ATTRHASH_IDX_EDGE, 1.0))
    return igraphmodule_handle_igraph_error();

  if (igraph_vector_init(&flow, 0)) {
    igraph_vector_destroy(&capacity_vector);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_vector_init(&cut, 0)) {
    igraph_vector_destroy(&capacity_vector);
    igraph_vector_destroy(&flow);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_vector_init(&partition, 0)) {
    igraph_vector_destroy(&capacity_vector);
    igraph_vector_destroy(&flow);
    igraph_vector_destroy(&cut);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_maxflow(&self->g, &result, &flow, &cut, &partition, 0,
        v1, v2, &capacity_vector)) {
    igraph_vector_destroy(&capacity_vector);
    igraph_vector_destroy(&flow);
    igraph_vector_destroy(&cut);
    igraph_vector_destroy(&partition);
    return igraphmodule_handle_igraph_error();
  }

  igraph_vector_destroy(&capacity_vector);

  flow_o = igraphmodule_vector_t_to_PyList(&flow, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&flow);
  
  if (flow_o == NULL) {
    igraph_vector_destroy(&cut);
    igraph_vector_destroy(&partition);
    return NULL;
  }

  cut_o = igraphmodule_vector_t_to_PyList(&cut, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&cut);
  
  if (cut_o == NULL) {
    igraph_vector_destroy(&partition);
    return NULL;
  }

  partition_o = igraphmodule_vector_t_to_PyList(&partition, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&partition);
  
  if (partition_o == NULL)
    return NULL;

  return Py_BuildValue("dOOO", (double)result, flow_o, cut_o, partition_o);
}

/**********************************************************************
 * Minimum cuts (edge separators)                                     *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Calculates all s-t cuts in a graph
 */
PyObject *igraphmodule_Graph_all_st_cuts(igraphmodule_GraphObject * self,
                                         PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "source", "target", NULL };
  igraph_integer_t source, target;
  igraph_vector_ptr_t cuts, partition1s;
  PyObject *source_o, *target_o;
  PyObject *cuts_o, *partition1s_o;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist,
                                   &source_o, &target_o))
    return NULL;

  if (igraphmodule_PyObject_to_vid(source_o, &source, &self->g))
    return NULL;
  if (igraphmodule_PyObject_to_vid(target_o, &target, &self->g))
    return NULL;

  if (igraph_vector_ptr_init(&partition1s, 0)) {
	return igraphmodule_handle_igraph_error();
  }
  if (igraph_vector_ptr_init(&cuts, 0)) {
    igraph_vector_ptr_destroy(&partition1s);
	return igraphmodule_handle_igraph_error();
  }

  if (igraph_all_st_cuts(&self->g, &cuts, &partition1s,
        source, target)) {
    igraph_vector_ptr_destroy(&cuts);
    igraph_vector_ptr_destroy(&partition1s);
    return igraphmodule_handle_igraph_error();
  }

  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&cuts, igraph_vector_destroy);
  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&partition1s, igraph_vector_destroy);

  cuts_o = igraphmodule_vector_ptr_t_to_PyList(&cuts, IGRAPHMODULE_TYPE_INT);
  igraph_vector_ptr_destroy_all(&cuts);
  if (cuts_o == NULL) {
    igraph_vector_ptr_destroy_all(&partition1s);
    return NULL;
  }

  partition1s_o = igraphmodule_vector_ptr_t_to_PyList(&partition1s, IGRAPHMODULE_TYPE_INT);
  igraph_vector_ptr_destroy_all(&partition1s);
  if (partition1s_o == NULL)
    return NULL;

  return Py_BuildValue("NN", cuts_o, partition1s_o);
}

/** \ingroup python_interface_graph
 * \brief Calculates all minimum s-t cuts in a graph
 */
PyObject *igraphmodule_Graph_all_st_mincuts(igraphmodule_GraphObject * self,
                                         PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "source", "target", "capacity", NULL };
  igraph_integer_t source, target;
  igraph_real_t value;
  igraph_vector_ptr_t cuts, partition1s;
  igraph_vector_t capacity_vector;
  PyObject *source_o, *target_o, *capacity_o = Py_None;
  PyObject *cuts_o, *partition1s_o;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOO", kwlist,
                                   &source_o, &target_o, &capacity_o))
    return NULL;

  if (igraphmodule_PyObject_to_vid(source_o, &source, &self->g))
    return NULL;
  if (igraphmodule_PyObject_to_vid(target_o, &target, &self->g))
    return NULL;

  if (igraph_vector_ptr_init(&partition1s, 0)) {
	return igraphmodule_handle_igraph_error();
  }
  if (igraph_vector_ptr_init(&cuts, 0)) {
    igraph_vector_ptr_destroy(&partition1s);
	return igraphmodule_handle_igraph_error();
  }

  if (igraphmodule_PyObject_to_attribute_values(capacity_o,
                                                &capacity_vector,
                                                self, ATTRHASH_IDX_EDGE, 1.0)) {
    igraph_vector_ptr_destroy(&cuts);
    igraph_vector_ptr_destroy(&partition1s);
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_all_st_mincuts(&self->g, &value, &cuts, &partition1s,
        source, target, &capacity_vector)) {
    igraph_vector_ptr_destroy(&cuts);
    igraph_vector_ptr_destroy(&partition1s);
    igraph_vector_destroy(&capacity_vector);
    return igraphmodule_handle_igraph_error();
  }

  igraph_vector_destroy(&capacity_vector);

  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&cuts, igraph_vector_destroy);
  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&partition1s, igraph_vector_destroy);

  cuts_o = igraphmodule_vector_ptr_t_to_PyList(&cuts, IGRAPHMODULE_TYPE_INT);
  igraph_vector_ptr_destroy_all(&cuts);
  if (cuts_o == NULL) {
    igraph_vector_ptr_destroy_all(&partition1s);
    return NULL;
  }

  partition1s_o = igraphmodule_vector_ptr_t_to_PyList(&partition1s, IGRAPHMODULE_TYPE_INT);
  igraph_vector_ptr_destroy_all(&partition1s);
  if (partition1s_o == NULL)
    return NULL;

  return Py_BuildValue("dNN", (double)value, cuts_o, partition1s_o);
}

/** \ingroup python_interface_graph
 * \brief Calculates the value of the minimum cut in the graph
 */
PyObject *igraphmodule_Graph_mincut_value(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "source", "target", "capacity", NULL };
  PyObject *capacity_object = Py_None;
  igraph_vector_t capacity_vector;
  igraph_real_t result, mincut;
  igraph_integer_t v1, v2;
  long vid1 = -1, vid2 = -1;
  long n;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|llO", kwlist,
                                   &vid1, &vid2, &capacity_object))
    return NULL;

  if (igraphmodule_PyObject_to_attribute_values(capacity_object,
                                                &capacity_vector,
                                                self, ATTRHASH_IDX_EDGE, 1.0))
    return igraphmodule_handle_igraph_error();

  v1 = (igraph_integer_t) vid1;
  v2 = (igraph_integer_t) vid2;
  if (v1 == -1 && v2 == -1) {
    if (igraph_mincut_value(&self->g, &result, &capacity_vector)) {
      igraph_vector_destroy(&capacity_vector);
      return igraphmodule_handle_igraph_error();
    }
  }
  else if (v1 == -1) {
    n = igraph_vcount(&self->g);
    result = -1;
    for (v1 = 0; v1 < n; v1++) {
      if (v2 == v1)
        continue;
      if (igraph_st_mincut_value(&self->g, &mincut, v1, v2, &capacity_vector)) {
        igraph_vector_destroy(&capacity_vector);
        return igraphmodule_handle_igraph_error();
      }
      if (result < 0 || result > mincut)
        result = mincut;
    }
    if (result < 0)
      result = 0.0;
  }
  else if (v2 == -1) {
    n = igraph_vcount(&self->g);
    result = -1;
    for (v2 = 0; v2 < n; v2++) {
      if (v2 == v1)
        continue;
      if (igraph_st_mincut_value(&self->g, &mincut, v1, v2, &capacity_vector)) {
        igraph_vector_destroy(&capacity_vector);
        return igraphmodule_handle_igraph_error();
      }
      if (result < 0.0 || result > mincut)
        result = mincut;
    }
    if (result < 0)
      result = 0.0;
  }
  else {
    if (igraph_st_mincut_value(&self->g, &result, v1, v2, &capacity_vector)) {
      igraph_vector_destroy(&capacity_vector);
      return igraphmodule_handle_igraph_error();
    }
  }

  igraph_vector_destroy(&capacity_vector);
  return Py_BuildValue("d", (double)result);
}

/** \ingroup python_interface_graph
 * \brief Calculates a minimum cut in a graph
 */
PyObject *igraphmodule_Graph_mincut(igraphmodule_GraphObject * self,
                                    PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "capacity", NULL };
  PyObject *capacity_object = Py_None, *cut_o, *part_o, *part2_o, *result;
  igraph_vector_t capacity_vector;
  igraph_real_t value;
  igraph_vector_t partition, partition2, cut;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist,
                                   &capacity_object))
    return NULL;

  if (igraphmodule_PyObject_to_attribute_values(capacity_object,
                                                &capacity_vector,
                                                self, ATTRHASH_IDX_EDGE, 1.0))
    return igraphmodule_handle_igraph_error();

  if (igraph_vector_init(&partition, 0)) {
    igraph_vector_destroy(&capacity_vector);
	return igraphmodule_handle_igraph_error();
  }
  if (igraph_vector_init(&partition2, 0)) {
    igraph_vector_destroy(&partition);
    igraph_vector_destroy(&capacity_vector);
	return igraphmodule_handle_igraph_error();
  }
  if (igraph_vector_init(&cut, 0)) {
    igraph_vector_destroy(&partition);
    igraph_vector_destroy(&partition2);
    igraph_vector_destroy(&capacity_vector);
	return igraphmodule_handle_igraph_error();
  }

  if (igraph_mincut(&self->g, &value, &partition, &partition2,
      &cut, &capacity_vector)) {
    igraph_vector_destroy(&cut);
    igraph_vector_destroy(&partition);
    igraph_vector_destroy(&partition2);
    igraph_vector_destroy(&capacity_vector);
    return igraphmodule_handle_igraph_error();
  }

  igraph_vector_destroy(&capacity_vector);

  cut_o=igraphmodule_vector_t_to_PyList(&cut, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&cut);
  if (!cut_o) {
    igraph_vector_destroy(&partition);
    igraph_vector_destroy(&partition2);
	return 0;
  }

  part_o=igraphmodule_vector_t_to_PyList(&partition, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&partition);
  if (!part_o) {
    Py_DECREF(cut_o);
    igraph_vector_destroy(&partition2);
    return 0;
  }

  part2_o=igraphmodule_vector_t_to_PyList(&partition2, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&partition2);
  if (!part2_o) {
    Py_DECREF(part_o);
	Py_DECREF(cut_o);
    return 0;
  }

  result = Py_BuildValue("dNNN", (double)value, cut_o, part_o, part2_o);
  return result;
}

/**********************************************************************
 * Vertex separators                                                  *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Returns all minimal s-t separators of a graph
 */
PyObject *igraphmodule_Graph_all_minimal_st_separators(
    igraphmodule_GraphObject * self) {
  PyObject* result_o;
  igraph_vector_ptr_t result;

  if (igraph_vector_ptr_init(&result, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_all_minimal_st_separators(&self->g, &result)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_ptr_destroy(&result);
    return NULL;
  }

  result_o = igraphmodule_vector_ptr_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&result, igraph_vector_destroy);
  igraph_vector_ptr_destroy_all(&result);

  return result_o;
}

/** \ingroup python_interface_graph
 * \brief Checks whether a given vertex set is a vertex separator
 */
PyObject *igraphmodule_Graph_is_separator(igraphmodule_GraphObject * self,
                                          PyObject * args, PyObject * kwds)
{
  PyObject* list = Py_None;
  igraph_bool_t result;
  igraph_vs_t vs;

  static char *kwlist[] = { "vertices", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &list))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, 0, 0)) {
    return NULL;
  }

  if (igraph_is_separator(&self->g, vs, &result)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    return NULL;
  }

  igraph_vs_destroy(&vs);

  if (result)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Checks whether a given vertex set is a minimal vertex separator
 */
PyObject *igraphmodule_Graph_is_minimal_separator(igraphmodule_GraphObject * self,
                                                  PyObject * args, PyObject * kwds)
{
  PyObject* list = Py_None;
  igraph_bool_t result;
  igraph_vs_t vs;

  static char *kwlist[] = { "vertices", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &list))
    return NULL;

  if (igraphmodule_PyObject_to_vs_t(list, &vs, &self->g, 0, 0)) {
    return NULL;
  }

  if (igraph_is_minimal_separator(&self->g, vs, &result)) {
    igraphmodule_handle_igraph_error();
    igraph_vs_destroy(&vs);
    return NULL;
  }

  igraph_vs_destroy(&vs);

  if (result)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

/** \ingroup python_interface_graph
 * \brief Returns the minimum size separators of the graph
 */
PyObject *igraphmodule_Graph_minimum_size_separators(
    igraphmodule_GraphObject * self) {
  PyObject* result_o;
  igraph_vector_ptr_t result;

  if (igraph_vector_ptr_init(&result, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_minimum_size_separators(&self->g, &result)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_ptr_destroy(&result);
    return NULL;
  }

  result_o = igraphmodule_vector_ptr_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&result, igraph_vector_destroy);
  igraph_vector_ptr_destroy_all(&result);

  return result_o;
}

/**********************************************************************
 * Cohesive blocks                                                    *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Calculates the cohesive block structure of a graph
 */
PyObject *igraphmodule_Graph_cohesive_blocks(igraphmodule_GraphObject *self,
    PyObject *args, PyObject *kwds) {
  PyObject *blocks_o, *cohesion_o, *parents_o, *result_o;
  igraph_vector_ptr_t blocks;
  igraph_vector_t cohesion, parents;

  if (igraph_vector_ptr_init(&blocks, 0)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&cohesion, 0)) {
    igraph_vector_ptr_destroy(&blocks);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_vector_init(&parents, 0)) {
    igraph_vector_ptr_destroy(&blocks);
    igraph_vector_destroy(&cohesion);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_cohesive_blocks(&self->g, &blocks, &cohesion, &parents, 0)) {
    igraph_vector_ptr_destroy(&blocks);
    igraph_vector_destroy(&cohesion);
    igraph_vector_destroy(&parents);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  blocks_o = igraphmodule_vector_ptr_t_to_PyList(&blocks, IGRAPHMODULE_TYPE_INT);
  IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&blocks, igraph_vector_destroy);
  igraph_vector_ptr_destroy_all(&blocks);
  if (blocks_o == NULL) {
    igraph_vector_destroy(&parents);
    igraph_vector_destroy(&cohesion);
    return NULL;
  }

  cohesion_o = igraphmodule_vector_t_to_PyList(&cohesion, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&cohesion);
  if (cohesion_o == NULL) {
    Py_DECREF(blocks_o);
    igraph_vector_destroy(&parents);
    return NULL;
  }

  parents_o = igraphmodule_vector_t_to_PyList(&parents, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&parents); 

  if (parents_o == NULL) {
    Py_DECREF(blocks_o);
    Py_DECREF(cohesion_o);
    return NULL;
  }

  result_o = Py_BuildValue("NNN", blocks_o, cohesion_o, parents_o);
  if (result_o == NULL) {
    Py_DECREF(parents_o);
    Py_DECREF(blocks_o);
    Py_DECREF(cohesion_o);
    return NULL;
  }
  return result_o;
}

/**********************************************************************
 * Cliques and independent sets                                       *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Find all or some cliques in a graph
 */
PyObject *igraphmodule_Graph_cliques(igraphmodule_GraphObject * self,
                                     PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "min", "max", NULL };
  PyObject *list, *item;
  long int min_size = 0, max_size = 0;
  long int i, j, n;
  igraph_vector_ptr_t result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ll", kwlist,
                                   &min_size, &max_size))
    return NULL;

  if (igraph_vector_ptr_init(&result, 0)) {
    PyErr_SetString(PyExc_MemoryError, "");
    return NULL;
  }

  if (igraph_cliques(&self->g, &result, (igraph_integer_t) min_size,
        (igraph_integer_t) max_size)) {
    igraph_vector_ptr_destroy(&result);
    return igraphmodule_handle_igraph_error();
  }

  n = (long)igraph_vector_ptr_size(&result);
  list = PyList_New(n);
  if (!list)
    return NULL;

  for (i = 0; i < n; i++) {
    igraph_vector_t *vec = (igraph_vector_t *) VECTOR(result)[i];
    item = igraphmodule_vector_t_to_PyTuple(vec);
    if (!item) {
      for (j = i; j < n; j++)
        igraph_vector_destroy((igraph_vector_t *) VECTOR(result)[j]);
      igraph_vector_ptr_destroy_all(&result);
      Py_DECREF(list);
      return NULL;
    }
    else {
      PyList_SET_ITEM(list, i, item);
    }
    igraph_vector_destroy(vec);
  }
  igraph_vector_ptr_destroy_all(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Find all largest cliques in a graph
 */
PyObject *igraphmodule_Graph_largest_cliques(igraphmodule_GraphObject * self)
{
  PyObject *list, *item;
  long int i, j, n;
  igraph_vector_ptr_t result;

  if (igraph_vector_ptr_init(&result, 0)) {
    PyErr_SetString(PyExc_MemoryError, "");
    return NULL;
  }

  if (igraph_largest_cliques(&self->g, &result)) {
    igraph_vector_ptr_destroy(&result);
    return igraphmodule_handle_igraph_error();
  }

  n = (long)igraph_vector_ptr_size(&result);
  list = PyList_New(n);
  if (!list)
    return NULL;

  for (i = 0; i < n; i++) {
    igraph_vector_t *vec = (igraph_vector_t *) VECTOR(result)[i];
    item = igraphmodule_vector_t_to_PyTuple(vec);
    if (!item) {
      for (j = i; j < n; j++)
        igraph_vector_destroy((igraph_vector_t *) VECTOR(result)[j]);
      igraph_vector_ptr_destroy_all(&result);
      Py_DECREF(list);
      return NULL;
    }
    else {
      PyList_SET_ITEM(list, i, item);
    }
    igraph_vector_destroy(vec);
  }
  igraph_vector_ptr_destroy_all(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Finds a maximum matching in a bipartite graph
 */
PyObject *igraphmodule_Graph_maximum_bipartite_matching(igraphmodule_GraphObject* self,
    PyObject* args, PyObject* kwds) {
  static char* kwlist[] = { "types", "weights", "eps", NULL };
  PyObject *types_o = Py_None, *weights_o = Py_None, *result_o;
  igraph_vector_bool_t* types = 0;
  igraph_vector_t* weights = 0;
  igraph_vector_long_t result;
  double eps = -1;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|Od", kwlist, &types_o,
        &weights_o, &eps))
    return NULL;

  if (eps < 0)
    eps = DBL_EPSILON * 1000;

  if (igraphmodule_attrib_to_vector_bool_t(types_o, self, &types, ATTRIBUTE_TYPE_VERTEX))
	return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights, ATTRIBUTE_TYPE_EDGE)) {
    if (types != 0) { igraph_vector_bool_destroy(types); free(types); }
	return NULL;
  }

  if (igraph_vector_long_init(&result, 0)) {
    if (types != 0) { igraph_vector_bool_destroy(types); free(types); }
    if (weights != 0) { igraph_vector_destroy(weights); free(weights); }
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_maximum_bipartite_matching(&self->g, types, 0, 0, &result, weights, eps)) {
    if (types != 0) { igraph_vector_bool_destroy(types); free(types); }
    if (weights != 0) { igraph_vector_destroy(weights); free(weights); }
    igraph_vector_long_destroy(&result);
    igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (types != 0) { igraph_vector_bool_destroy(types); free(types); }
  if (weights != 0) { igraph_vector_destroy(weights); free(weights); }

  result_o = igraphmodule_vector_long_t_to_PyList(&result);
  igraph_vector_long_destroy(&result);

  return result_o;
}

/** \ingroup python_interface_graph
 * \brief Find all maximal cliques in a graph
 */
PyObject *igraphmodule_Graph_maximal_cliques(igraphmodule_GraphObject * self,
    PyObject* args, PyObject* kwds) {
  static char* kwlist[] = { "min", "max", NULL };
  PyObject *list, *item;
  long int i = 0, j = 0;
  Py_ssize_t n;
  igraph_vector_ptr_t result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ll", kwlist, &i, &j))
    return NULL;

  if (igraph_vector_ptr_init(&result, 0)) {
    PyErr_SetString(PyExc_MemoryError, "");
    return NULL;
  }

  if (igraph_maximal_cliques(&self->g, &result, (igraph_integer_t) i, (igraph_integer_t) j)) {
    igraph_vector_ptr_destroy(&result);
    return igraphmodule_handle_igraph_error();
  }

  n = (Py_ssize_t)igraph_vector_ptr_size(&result);
  list = PyList_New(n);
  if (!list)
    return NULL;

  for (i = 0; i < n; i++) {
    igraph_vector_t *vec = (igraph_vector_t *) VECTOR(result)[i];
    item = igraphmodule_vector_t_to_PyTuple(vec);
    if (!item) {
      for (j = i; j < n; j++)
        igraph_vector_destroy((igraph_vector_t *) VECTOR(result)[j]);
      igraph_vector_ptr_destroy_all(&result);
      Py_DECREF(list);
      return NULL;
    }
    else {
      PyList_SET_ITEM(list, i, item);
    }
    igraph_vector_destroy(vec);
  }
  igraph_vector_ptr_destroy_all(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Returns the clique number of the graph
 */
PyObject *igraphmodule_Graph_clique_number(igraphmodule_GraphObject * self)
{
  PyObject *result;
  igraph_integer_t i;

  if (igraph_clique_number(&self->g, &i))
    return igraphmodule_handle_igraph_error();

  result = PyInt_FromLong((long)i);
  return result;
}

/** \ingroup python_interface_graph
 * \brief Find all or some independent vertex sets in a graph
 */
PyObject *igraphmodule_Graph_independent_vertex_sets(igraphmodule_GraphObject
                                                     * self, PyObject * args,
                                                     PyObject * kwds)
{
  static char *kwlist[] = { "min", "max", NULL };
  PyObject *list, *item;
  long int min_size = 0, max_size = 0;
  long int i, j, n;
  igraph_vector_ptr_t result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ll", kwlist,
                                   &min_size, &max_size))
    return NULL;

  if (igraph_vector_ptr_init(&result, 0)) {
    PyErr_SetString(PyExc_MemoryError, "");
    return NULL;
  }

  if (igraph_independent_vertex_sets(&self->g, &result,
        (igraph_integer_t) min_size, (igraph_integer_t) max_size)) {
    igraph_vector_ptr_destroy(&result);
    return igraphmodule_handle_igraph_error();
  }

  n = (long)igraph_vector_ptr_size(&result);
  list = PyList_New(n);
  if (!list)
    return NULL;

  for (i = 0; i < n; i++) {
    igraph_vector_t *vec = (igraph_vector_t *) VECTOR(result)[i];
    item = igraphmodule_vector_t_to_PyTuple(vec);
    if (!item) {
      for (j = i; j < n; j++)
        igraph_vector_destroy((igraph_vector_t *) VECTOR(result)[j]);
      igraph_vector_ptr_destroy_all(&result);
      Py_DECREF(list);
      return NULL;
    }
    else {
      PyList_SET_ITEM(list, i, item);
    }
    igraph_vector_destroy(vec);
  }
  igraph_vector_ptr_destroy_all(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Find all largest independent_vertex_sets in a graph
 */
PyObject
  *igraphmodule_Graph_largest_independent_vertex_sets(igraphmodule_GraphObject
                                                      * self)
{
  PyObject *list, *item;
  long int i, j, n;
  igraph_vector_ptr_t result;

  if (igraph_vector_ptr_init(&result, 0)) {
    PyErr_SetString(PyExc_MemoryError, "");
    return NULL;
  }

  if (igraph_largest_independent_vertex_sets(&self->g, &result)) {
    igraph_vector_ptr_destroy(&result);
    return igraphmodule_handle_igraph_error();
  }

  n = (long)igraph_vector_ptr_size(&result);
  list = PyList_New(n);
  if (!list)
    return NULL;

  for (i = 0; i < n; i++) {
    igraph_vector_t *vec = (igraph_vector_t *) VECTOR(result)[i];
    item = igraphmodule_vector_t_to_PyTuple(vec);
    if (!item) {
      for (j = i; j < n; j++)
        igraph_vector_destroy((igraph_vector_t *) VECTOR(result)[j]);
      igraph_vector_ptr_destroy_all(&result);
      Py_DECREF(list);
      return NULL;
    }
    else {
      PyList_SET_ITEM(list, i, item);
    }
    igraph_vector_destroy(vec);
  }
  igraph_vector_ptr_destroy_all(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Find all maximal independent vertex sets in a graph
 */
PyObject
  *igraphmodule_Graph_maximal_independent_vertex_sets(igraphmodule_GraphObject
                                                      * self)
{
  PyObject *list, *item;
  long int i, j, n;
  igraph_vector_ptr_t result;

  if (igraph_vector_ptr_init(&result, 0)) {
    PyErr_SetString(PyExc_MemoryError, "");
    return NULL;
  }

  if (igraph_maximal_independent_vertex_sets(&self->g, &result)) {
    igraph_vector_ptr_destroy(&result);
    return igraphmodule_handle_igraph_error();
  }

  n = (long)igraph_vector_ptr_size(&result);
  list = PyList_New(n);
  if (!list)
    return NULL;

  for (i = 0; i < n; i++) {
    igraph_vector_t *vec = (igraph_vector_t *) VECTOR(result)[i];
    item = igraphmodule_vector_t_to_PyTuple(vec);
    if (!item) {
      for (j = i; j < n; j++)
        igraph_vector_destroy((igraph_vector_t *) VECTOR(result)[j]);
      igraph_vector_ptr_destroy_all(&result);
      Py_DECREF(list);
      return NULL;
    }
    else {
      PyList_SET_ITEM(list, i, item);
    }
    igraph_vector_destroy(vec);
  }
  igraph_vector_ptr_destroy_all(&result);

  return list;
}

/** \ingroup python_interface_graph
 * \brief Returns the independence number of the graph
 */
PyObject *igraphmodule_Graph_independence_number(igraphmodule_GraphObject *
                                                 self)
{
  PyObject *result;
  igraph_integer_t i;

  if (igraph_independence_number(&self->g, &i))
    return igraphmodule_handle_igraph_error();

  result = PyInt_FromLong((long)i);
  return result;
}

/**********************************************************************
 * K-core decomposition                                               *
 **********************************************************************/

/** \ingroup python_interface_graph
 * \brief Returns the corenesses of the graph vertices
 * \return a new PyCObject
 */
PyObject *igraphmodule_Graph_coreness(igraphmodule_GraphObject * self,
                                      PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "mode", NULL };
  igraph_neimode_t mode = IGRAPH_ALL;
  igraph_vector_t result;
  PyObject *o, *mode_o = Py_None;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &mode_o))
    return NULL;

  if (igraphmodule_PyObject_to_neimode_t(mode_o, &mode)) return NULL;

  if (igraph_vector_init(&result, igraph_vcount(&self->g)))
    return igraphmodule_handle_igraph_error();

  if (igraph_coreness(&self->g, &result, mode)) {
    igraph_vector_destroy(&result);
    return igraphmodule_handle_igraph_error();
  }

  o = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  return o;
}

/**********************************************************************
 * Community structure detection and related routines                 *
 **********************************************************************/

/**
 * Modularity calculation
 */
PyObject *igraphmodule_Graph_modularity(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"membership", "weights", 0};
  igraph_vector_t membership;
  igraph_vector_t *weights=0;
  igraph_real_t modularity;
  PyObject *mvec, *wvec=Py_None;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &mvec, &wvec))
    return NULL;

  if (igraphmodule_PyObject_to_vector_t(mvec, &membership, 1, 0)) return NULL;
  if (igraphmodule_attrib_to_vector_t(wvec, self, &weights, ATTRIBUTE_TYPE_EDGE)){
    igraph_vector_destroy(&membership);
    return NULL;
  }
  if (igraph_modularity(&self->g, &membership, &modularity, weights)) {
    igraph_vector_destroy(&membership);
    if (weights) { igraph_vector_destroy(weights); free(weights); }
    return NULL;
  }
  igraph_vector_destroy(&membership);
  if (weights) { igraph_vector_destroy(weights); free(weights); }
  return Py_BuildValue("d", (double)modularity);
}

/**
 * Newman's edge betweenness method
 */
PyObject *igraphmodule_Graph_community_edge_betweenness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "directed", "weights", NULL };
  PyObject *directed = Py_True;
  PyObject *weights_o = Py_None;
  PyObject *res, *qs, *ms;
  igraph_matrix_t merges;
  igraph_vector_t q;
  igraph_vector_t *weights = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &directed, &weights_o))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) return NULL;

  if (igraph_matrix_init(&merges, 0, 0)) {
    if (weights != 0) {
      igraph_vector_destroy(weights); free(weights);
    }
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_vector_init(&q, 0)) {
    igraph_matrix_destroy(&merges);
    if (weights != 0) {
      igraph_vector_destroy(weights); free(weights);
    }
    return igraphmodule_handle_igraph_error();
  }

  if (igraph_community_edge_betweenness(&self->g,
        /* removed_edges = */ 0,
        /* edge_betweenness = */ 0,
        /* merges = */ &merges,
        /* bridges = */ 0,
        /* modularity = */ &q,
        /* membership = */ 0,
        PyObject_IsTrue(directed),
        weights)) {
    igraphmodule_handle_igraph_error();
    if (weights != 0) {
      igraph_vector_destroy(weights); free(weights);
    }
    igraph_matrix_destroy(&merges);
    igraph_vector_destroy(&q);
    return NULL;
  }

  if (weights != 0) {
    igraph_vector_destroy(weights); free(weights);
  }

  qs=igraphmodule_vector_t_to_PyList(&q, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&q);
  if (!qs) {
    igraph_matrix_destroy(&merges);
    return NULL;
  }

  ms=igraphmodule_matrix_t_to_PyList(&merges, IGRAPHMODULE_TYPE_INT);
  igraph_matrix_destroy(&merges);

  if (ms == NULL) {
    Py_DECREF(qs);
    return NULL;
  }

  res=Py_BuildValue("NN", ms, qs); /* steals references */

  return res;
}

/**
 * Newman's leading eigenvector method, precise implementation
 */
PyObject *igraphmodule_Graph_community_leading_eigenvector(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "n", "arpack_options", NULL };
  long int n=-1;
  PyObject *cl, *res, *merges;
  igraph_vector_t members;
  igraph_matrix_t m;
  igraph_real_t q;
  igraphmodule_ARPACKOptionsObject *arpack_options;
  PyObject *arpack_options_o = igraphmodule_arpack_options_default;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lO!", kwlist,
    &n, &igraphmodule_ARPACKOptionsType, &arpack_options_o)) {
    return NULL;
  }

  if (igraph_vector_init(&members, 0))
    return igraphmodule_handle_igraph_error();

  if (igraph_matrix_init(&m, 0, 0))
    return igraphmodule_handle_igraph_error();

  if (n<0)
    n = igraph_vcount(&self->g);
  else
    n -= 1;

  arpack_options = (igraphmodule_ARPACKOptionsObject*)arpack_options_o;
  if (igraph_community_leading_eigenvector(&self->g, &m, &members, (igraph_integer_t) n,
      igraphmodule_ARPACKOptions_get(arpack_options), &q, 0, 0, 0, 0, 0, 0)){
    igraph_matrix_destroy(&m);
    igraph_vector_destroy(&members);
    return igraphmodule_handle_igraph_error();
  }

  cl = igraphmodule_vector_t_to_PyList(&members, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&members);
  if (cl == 0) {
    igraph_matrix_destroy(&m);
    return 0;
  }

  merges = igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_INT);
  igraph_matrix_destroy(&m);
  if (merges == 0)
    return 0;

  res=Py_BuildValue("NNd", cl, merges, (double)q);

  return res;
}

/**
 * Clauset et al's greedy modularity optimization algorithm
 */
PyObject *igraphmodule_Graph_community_fastgreedy(igraphmodule_GraphObject * self,
  PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "weights", NULL };
  PyObject *ms, *qs, *res, *weights = Py_None;
  igraph_matrix_t merges;
  igraph_vector_t q, *ws=0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &weights)) {
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights, self, &ws, ATTRIBUTE_TYPE_EDGE))
    return NULL;

  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&q, 0);
  if (igraph_community_fastgreedy(&self->g, ws, &merges, &q, 0)) {
    if (ws) {
      igraph_vector_destroy(ws); free(ws);
    }
    igraph_vector_destroy(&q);
    igraph_matrix_destroy(&merges);
    return igraphmodule_handle_igraph_error();
  }
  if (ws) {
    igraph_vector_destroy(ws); free(ws);
  }

  qs=igraphmodule_vector_t_to_PyList(&q, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&q);
  if (!qs) {
    igraph_matrix_destroy(&merges);
    return NULL;
  }

  ms=igraphmodule_matrix_t_to_PyList(&merges, IGRAPHMODULE_TYPE_INT);
  igraph_matrix_destroy(&merges);

  if (ms == NULL) {
    Py_DECREF(qs);
    return NULL;
  }

  res=Py_BuildValue("NN", ms, qs); /* steals references */

  return res;
}

/**
 * Infomap community detection algorithm of Martin Rosvall and Carl T. Bergstrom,
 */
PyObject *igraphmodule_Graph_community_infomap(igraphmodule_GraphObject * self,
  PyObject * args, PyObject * kwds)
{
  static char *kwlist[] = { "edge_weights", "vertex_weights", "trials", NULL };
  PyObject *e_weights = Py_None, *v_weights = Py_None;
  unsigned int nb_trials = 10;
  igraph_vector_t *e_ws = 0, *v_ws = 0;
  
  igraph_vector_t membership;
  PyObject *res = Py_False;
  igraph_real_t codelength;
  
  
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOI", kwlist, &e_weights,
        &v_weights, &nb_trials)) {
    return NULL;
  }

  if (igraph_vector_init(&membership, igraph_vcount(&self->g))) {
	igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(e_weights, self, &e_ws, ATTRIBUTE_TYPE_EDGE)) {
    igraph_vector_destroy(&membership);
    return NULL;
  }
  
  if (igraphmodule_attrib_to_vector_t(v_weights, self, &v_ws, ATTRIBUTE_TYPE_VERTEX)){
    igraph_vector_destroy(&membership);
    if (e_ws) {
      igraph_vector_destroy(e_ws);
      free(e_ws);
    }
    return NULL;
  }
  
  if (igraph_community_infomap(/*in */ &self->g, 
                                    /*e_weight=*/ e_ws, /*v_weight=*/ v_ws,
                                    /*nb_trials=*/nb_trials,
                              /*out*/ &membership, &codelength)) {
	igraphmodule_handle_igraph_error();
	igraph_vector_destroy(&membership);
    if (e_ws) {
      igraph_vector_destroy(e_ws);
      free(e_ws);
    }
    if (v_ws) {
      igraph_vector_destroy(v_ws);
      free(v_ws);
    }
    return NULL;
  }
  
  if (e_ws) {
    igraph_vector_destroy(e_ws);
    free(e_ws);
  }
  
  if (v_ws) {
    igraph_vector_destroy(v_ws);
    free(v_ws);
  }
  
  res = igraphmodule_vector_t_to_PyList(&membership, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&membership);
  
  if (!res)
	return NULL;
  
  return Py_BuildValue("Nd", res, (double)codelength);
}


/**
 * The label propagation algorithm of Raghavan et al
 */
PyObject *igraphmodule_Graph_community_label_propagation(
    igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = { "weights", "initial", "fixed", NULL };
  PyObject *weights_o = Py_None, *initial_o = Py_None, *fixed_o = Py_None;
  PyObject *result;
  igraph_vector_t membership, *ws = 0, *initial = 0;
  igraph_vector_bool_t fixed;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist, &weights_o, &initial_o, &fixed_o)) {
    return NULL;
  }

  if (fixed_o != Py_None) {
    if (igraphmodule_PyObject_to_vector_bool_t(fixed_o, &fixed))
      return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &ws, ATTRIBUTE_TYPE_EDGE)) {
    if (fixed_o != Py_None) igraph_vector_bool_destroy(&fixed);
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(initial_o, self, &initial, ATTRIBUTE_TYPE_VERTEX)){
    if (fixed_o != Py_None) igraph_vector_bool_destroy(&fixed);
    if (ws) { igraph_vector_destroy(ws); free(ws); }
    return NULL;
  }

  igraph_vector_init(&membership, igraph_vcount(&self->g));
  if (igraph_community_label_propagation(&self->g, &membership,
        ws, initial, (fixed_o != Py_None ? &fixed : 0), 0)) {
    if (fixed_o != Py_None) igraph_vector_bool_destroy(&fixed);
    if (ws) { igraph_vector_destroy(ws); free(ws); }
    if (initial) { igraph_vector_destroy(initial); free(initial); }
    igraph_vector_destroy(&membership);
    return igraphmodule_handle_igraph_error();
  }

  if (fixed_o != Py_None) igraph_vector_bool_destroy(&fixed);
  if (ws) { igraph_vector_destroy(ws); free(ws); }
  if (initial) { igraph_vector_destroy(initial); free(initial); }

  result=igraphmodule_vector_t_to_PyList(&membership, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&membership);

  return result;
}

/**
 * Multilevel algorithm of Blondel et al
 */
PyObject *igraphmodule_Graph_community_multilevel(igraphmodule_GraphObject *self,
  PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = { "weights", "return_levels", NULL };
  PyObject *return_levels = Py_False;
  PyObject *mss, *qs, *res, *weights = Py_None;
  igraph_matrix_t memberships;
  igraph_vector_t membership, modularity;
  igraph_vector_t *ws;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &weights, &return_levels)) {
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights, self, &ws, ATTRIBUTE_TYPE_EDGE))
    return NULL;

  igraph_matrix_init(&memberships, 0, 0);
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&modularity, 0);

  if (igraph_community_multilevel(&self->g, ws, &membership, &memberships,
        &modularity)) {
    if (ws) { igraph_vector_destroy(ws); free(ws); }
    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&modularity);
    igraph_matrix_destroy(&memberships);
    return igraphmodule_handle_igraph_error();
  }

  if (ws) { igraph_vector_destroy(ws); free(ws); }

  qs=igraphmodule_vector_t_to_PyList(&modularity, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&modularity);
  if (!qs) {
    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&memberships);
    return NULL;
  }

  if (PyObject_IsTrue(return_levels)) {
    mss=igraphmodule_matrix_t_to_PyList(&memberships, IGRAPHMODULE_TYPE_INT);
    if (!mss) {
      res = mss;
    } else {
      res=Py_BuildValue("NN", mss, qs); /* steals references */
    }
  } else {
    res=igraphmodule_vector_t_to_PyList(&membership, IGRAPHMODULE_TYPE_INT);
  }

  igraph_vector_destroy(&membership);
  igraph_matrix_destroy(&memberships);

  return res;
}

/**
 * Optimal modularity by integer programming
 */
PyObject *igraphmodule_Graph_community_optimal_modularity(
	igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = { "verbose", NULL };
  igraph_real_t modularity;
  igraph_vector_t membership;
  PyObject *res, *verbose_o = Py_False;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &verbose_o)) {
    return NULL;
  }

  if (igraph_vector_init(&membership, igraph_vcount(&self->g))) {
	igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraph_community_optimal_modularity(&self->g, &modularity, &membership,
        PyObject_IsTrue(verbose_o))) {
	igraphmodule_handle_igraph_error();
	igraph_vector_destroy(&membership);
    return NULL;
  }

  res = igraphmodule_vector_t_to_PyList(&membership, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&membership);

  if (!res)
	return NULL;

  return Py_BuildValue("Nd", res, (double)modularity);
}

/**
 * Spinglass community detection method of Reichardt & Bornholdt
 */
PyObject *igraphmodule_Graph_community_spinglass(igraphmodule_GraphObject *self,
        PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"weights", "spins", "parupdate",
      "start_temp", "stop_temp", "cool_fact", "update_rule",
      "gamma", "implementation", "lambda_", NULL};
  PyObject *weights_o = Py_None;
  PyObject *parupdate_o = Py_False;
  PyObject *update_rule_o = Py_None;
  PyObject *impl_o = Py_None;
  PyObject *res;

  long int spins = 25;
  double start_temp = 1.0;
  double stop_temp = 0.01;
  double cool_fact = 0.99;
  igraph_spinglass_implementation_t impl = IGRAPH_SPINCOMM_IMP_ORIG;
  igraph_spincomm_update_t update_rule = IGRAPH_SPINCOMM_UPDATE_CONFIG;
  double gamma = 1;
  double lambda = 1;
  igraph_vector_t *weights = 0, membership;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OlOdddOdOd", kwlist,
        &weights_o, &spins, &parupdate_o, &start_temp, &stop_temp,
        &cool_fact, &update_rule_o, &gamma, &impl_o, &lambda))
    return NULL;

  if (igraphmodule_PyObject_to_spincomm_update_t(update_rule_o, &update_rule)) {
    return NULL;
  }

  if (igraphmodule_PyObject_to_spinglass_implementation_t(impl_o, &impl)) {
    return NULL;
  }

  if (igraph_vector_init(&membership, igraph_vcount(&self->g))) {
	igraphmodule_handle_igraph_error();
    return NULL;
  }

  if (igraphmodule_attrib_to_vector_t(weights_o, self, &weights,
	  ATTRIBUTE_TYPE_EDGE)) {
    igraph_vector_destroy(&membership);
    return NULL;
  }

  if (igraph_community_spinglass(&self->g, weights,
              0, 0, &membership, 0, (igraph_integer_t) spins,
              PyObject_IsTrue(parupdate_o),
              start_temp, stop_temp, cool_fact,
              update_rule, gamma, impl, lambda)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&membership);
    if (weights != 0) {
      igraph_vector_destroy(weights);
      free(weights);
    }
    return NULL;
  }

  if (weights != 0) {
    igraph_vector_destroy(weights);
    free(weights);
  }

  res = igraphmodule_vector_t_to_PyList(&membership, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&membership);

  return res;
}

/**
 * Walktrap community detection of Latapy & Pons
 */
PyObject *igraphmodule_Graph_community_walktrap(igraphmodule_GraphObject * self,
  PyObject * args, PyObject * kwds) {
  static char *kwlist[] = { "weights", "steps", NULL };
  PyObject *ms, *qs, *res, *weights = Py_None;
  igraph_matrix_t merges;
  int steps=4;
  igraph_vector_t q, *ws=0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Oi", kwlist, &weights,
      &steps))
    return NULL;

  if (igraphmodule_attrib_to_vector_t(weights, self, &ws, ATTRIBUTE_TYPE_EDGE))
    return NULL;

  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&q, 0);

  if (igraph_community_walktrap(&self->g, ws, steps, &merges, &q, 0)) {
    if (ws) {
      igraph_vector_destroy(ws); free(ws);
    }
    igraph_vector_destroy(&q);
    igraph_matrix_destroy(&merges);
    return igraphmodule_handle_igraph_error();
  }
  if (ws) {
    igraph_vector_destroy(ws); free(ws);
  }

  qs = igraphmodule_vector_t_to_PyList(&q, IGRAPHMODULE_TYPE_FLOAT);
  igraph_vector_destroy(&q);
  if (!qs) {
    igraph_matrix_destroy(&merges);
    return NULL;
  }

  ms = igraphmodule_matrix_t_to_PyList(&merges, IGRAPHMODULE_TYPE_INT);
  igraph_matrix_destroy(&merges);

  if (ms == NULL) {
    Py_DECREF(qs);
    return NULL;
  }

  res=Py_BuildValue("NN", ms, qs); /* steals references */

  return res;
}

/* }}} */

/**********************************************************************
 * Special internal methods that you won't need to mess around with   *
 **********************************************************************/

/** \defgroup python_interface_internal Internal functions
 * \ingroup python_interface */

#ifdef IGRAPH_PYTHON3
PyObject *igraphmodule_Graph___graph_as_capsule__(igraphmodule_GraphObject *
                                                  self, PyObject * args,
                                                  PyObject * kwds)
{
  return PyCapsule_New((void *)&self->g, 0, 0);
}
#else
/** \ingroup python_interface_internal
 * \brief Returns the encapsulated igraph graph as a PyCObject
 * \return a new PyCObject
 */
PyObject *igraphmodule_Graph___graph_as_cobject__(igraphmodule_GraphObject *
                                                  self, PyObject * args,
                                                  PyObject * kwds)
{
  return PyCObject_FromVoidPtr((void *)&self->g, 0);
}
#endif

/** \ingroup python_interface_internal
 * \brief Registers a destructor to be called when the object is destroyed
 * \return the previous destructor (if any)
 * Unimplemented.
 */
PyObject *igraphmodule_Graph___register_destructor__(igraphmodule_GraphObject
                                                     * self, PyObject * args,
                                                     PyObject * kwds)
{
  char *kwlist[] = { "destructor", NULL };
  PyObject *destructor = NULL, *result;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &destructor))
    return NULL;

  if (!PyCallable_Check(destructor)) {
    PyErr_SetString(PyExc_TypeError, "The destructor must be callable!");
    return NULL;
  }

  result = self->destructor;
  self->destructor = destructor;
  Py_INCREF(self->destructor);

  if (!result)
    Py_RETURN_NONE;

  return result;
}

/** \ingroup python_interface
 * \brief Member list of the \c igraph.Graph object type
 */
#define OFF(x) offsetof(igraphmodule_GraphObject, x)

/** \ingroup python_interface
 * \brief Method list of the \c igraph.Graph object type
 */
struct PyMethodDef igraphmodule_Graph_methods[] = {
  ////////////////////////////
  // BASIC IGRAPH INTERFACE //
  ////////////////////////////

  // interface to igraph_vcount
  {"vcount", (PyCFunction) igraphmodule_Graph_vcount,
   METH_NOARGS,
   "vcount()\n\n"
   "Counts the number of vertices.\n"
   "@return: the number of vertices in the graph.\n" "@rtype: integer"},

  // interface to igraph_ecount
  {"ecount", (PyCFunction) igraphmodule_Graph_ecount,
   METH_NOARGS,
   "ecount()\n\n"
   "Counts the number of edges.\n"
   "@return: the number of edges in the graph.\n" "@rtype: integer"},

  // interface to igraph_is_dag
  {"is_dag", (PyCFunction) igraphmodule_Graph_is_dag,
   METH_NOARGS,
   "is_dag()\n\n"
   "Checks whether the graph is a DAG (directed acyclic graph).\n\n"
   "A DAG is a directed graph with no directed cycles.\n\n"
   "@return: C{True} if it is a DAG, C{False} otherwise.\n"
   "@rtype: boolean"},

  // interface to igraph_is_directed
  {"is_directed", (PyCFunction) igraphmodule_Graph_is_directed,
   METH_NOARGS,
   "is_directed()\n\n"
   "Checks whether the graph is directed.\n"
   "@return: C{True} if it is directed, C{False} otherwise.\n"
   "@rtype: boolean"},

  // interface to igraph_is_simple
  {"is_simple", (PyCFunction) igraphmodule_Graph_is_simple,
   METH_NOARGS,
   "is_simple()\n\n"
   "Checks whether the graph is simple (no loop or multiple edges).\n\n"
   "@return: C{True} if it is simple, C{False} otherwise.\n"
   "@rtype: boolean"},

  /* interface to igraph_add_vertices */
  {"add_vertices", (PyCFunction) igraphmodule_Graph_add_vertices,
   METH_VARARGS,
   "add_vertices(n)\n\n"
   "Adds vertices to the graph.\n\n"
   "@param n: the number of vertices to be added\n"},

  /* interface to igraph_delete_vertices */
  {"delete_vertices", (PyCFunction) igraphmodule_Graph_delete_vertices,
   METH_VARARGS,
   "delete_vertices(vs)\n\n"
   "Deletes vertices and all its edges from the graph.\n\n"
   "@param vs: a single vertex ID or the list of vertex IDs\n"
   "  to be deleted.\n"},

  /* interface to igraph_add_edges */
  {"add_edges", (PyCFunction) igraphmodule_Graph_add_edges,
   METH_VARARGS,
   "add_edges(es)\n\n"
   "Adds edges to the graph.\n\n"
   "@param es: the list of edges to be added. Every edge is\n"
   "  represented with a tuple, containing the vertex IDs of the\n"
   "  two endpoints. Vertices are enumerated from zero.\n"},

  /* interface to igraph_delete_edges */
  {"delete_edges", (PyCFunction) igraphmodule_Graph_delete_edges,
   METH_VARARGS | METH_KEYWORDS,
   "delete_edges(es)\n\n"
   "Removes edges from the graph.\n\n"
   "All vertices will be kept, even if they lose all their edges.\n"
   "Nonexistent edges will be silently ignored.\n\n"
   "@param es: the list of edges to be removed. Edges are identifed by\n"
   "  edge IDs. L{EdgeSeq} objects are also accepted here.\n"},

  /* interface to igraph_degree */
  {"degree", (PyCFunction) igraphmodule_Graph_degree,
   METH_VARARGS | METH_KEYWORDS,
   "degree(vertices, mode=ALL, loops=True)\n\n"
   "Returns some vertex degrees from the graph.\n\n"
   "This method accepts a single vertex ID or a list of vertex IDs as a\n"
   "parameter, and returns the degree of the given vertices (in the\n"
   "form of a single integer or a list, depending on the input\n"
   "parameter).\n"
   "\n"
   "@param vertices: a single vertex ID or a list of vertex IDs\n"
   "@param mode: the type of degree to be returned (L{OUT} for\n"
   "  out-degrees, L{IN} IN for in-degrees or L{ALL} for the sum of\n"
   "  them).\n" "@param loops: whether self-loops should be counted.\n"},

  /* interface to igraph_strength */
  {"strength", (PyCFunction) igraphmodule_Graph_strength,
   METH_VARARGS | METH_KEYWORDS,
   "strength(vertices, mode=ALL, loops=True, weights=None)\n\n"
   "Returns the strength (weighted degree) of some vertices from the graph\n\n"
   "This method accepts a single vertex ID or a list of vertex IDs as a\n"
   "parameter, and returns the strength (that is, the sum of the weights\n"
   "of all incident edges) of the given vertices (in the\n"
   "form of a single integer or a list, depending on the input\n"
   "parameter).\n"
   "\n"
   "@param vertices: a single vertex ID or a list of vertex IDs\n"
   "@param mode: the type of degree to be returned (L{OUT} for\n"
   "  out-degrees, L{IN} IN for in-degrees or L{ALL} for the sum of\n"
   "  them).\n"
   "@param loops: whether self-loops should be counted.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
  },

  /* interface to igraph_is_loop */
  {"is_loop", (PyCFunction) igraphmodule_Graph_is_loop,
   METH_VARARGS | METH_KEYWORDS,
   "is_loop(edges=None)\n\n"
   "Checks whether a specific set of edges contain loop edges\n\n"
   "@param edges: edge indices which we want to check. If C{None}, all\n"
   "  edges are checked.\n"
   "@return: a list of booleans, one for every edge given\n"},

  /* interface to igraph_is_multiple */
  {"is_multiple", (PyCFunction) igraphmodule_Graph_is_multiple,
   METH_VARARGS | METH_KEYWORDS,
   "is_multiple(edges=None)\n\n"
   "Checks whether an edge is a multiple edge.\n\n"
   "Also works for a set of edges -- in this case, every edge is checked\n"
   "one by one. Note that if there are multiple edges going between a\n"
   "pair of vertices, there is always one of them that is I{not}\n"
   "reported as multiple (only the others). This allows one to easily\n"
   "detect the edges that have to be deleted in order to make the graph\n"
   "free of multiple edges.\n\n"
   "@param edges: edge indices which we want to check. If C{None}, all\n"
   "  edges are checked.\n"
   "@return: a list of booleans, one for every edge given\n"},

  /* interface to igraph_has_multiple */
  {"has_multiple", (PyCFunction) igraphmodule_Graph_has_multiple,
   METH_NOARGS,
   "has_multiple()\n\n"
   "Checks whether the graph has multiple edges.\n\n"
   "@return: C{True} if the graph has at least one multiple edge,\n"
   "         C{False} otherwise.\n"
   "@rtype: boolean"},

  /* interface to igraph_is_mutual */
  {"is_mutual", (PyCFunction) igraphmodule_Graph_is_mutual,
   METH_VARARGS | METH_KEYWORDS,
   "is_mutual(edges=None)\n\n"
   "Checks whether an edge has an opposite pair.\n\n"
   "Also works for a set of edges -- in this case, every edge is checked\n"
   "one by one. The result will be a list of booleans (or a single boolean\n"
   "if only an edge index is supplied), every boolean corresponding to an\n"
   "edge in the edge set supplied. C{True} is returned for a given edge\n"
   "M{a} --> M{b} if there exists another edge M{b} --> M{a} in the\n"
   "original graph (not the given edge set!). All edges in an undirected\n"
   "graph are mutual. In case there are multiple edges between M{a}\n"
   "and M{b}, it is enough to have at least one edge in either direction\n"
   "to report all edges between them as mutual, so the multiplicity\n"
   "of edges do not matter.\n\n"
   "@param edges: edge indices which we want to check. If C{None}, all\n"
   "  edges are checked.\n"
   "@return: a list of booleans, one for every edge given\n"},

  /* interface to igraph_count_multiple */
  {"count_multiple", (PyCFunction) igraphmodule_Graph_count_multiple,
   METH_VARARGS | METH_KEYWORDS,
   "count_multiple(edges=None)\n\n"
   "Counts the multiplicities of the given edges.\n\n"
   "@param edges: edge indices for which we want to count their\n"
   "  multiplicity. If C{None}, all edges are counted.\n"
   "@return: the multiplicities of the given edges as a list.\n"},

  /* interface to igraph_neighbors */
  {"neighbors", (PyCFunction) igraphmodule_Graph_neighbors,
   METH_VARARGS | METH_KEYWORDS,
   "neighbors(vertex, mode=ALL)\n\n"
   "Returns adjacent vertices to a given vertex.\n\n"
   "@param vertex: a vertex ID\n"
   "@param mode: whether to return only predecessors (L{OUT}),\n"
   "  successors (L{IN}) or both (L{ALL}). Ignored for undirected\n"
   "  graphs."},

  {"successors", (PyCFunction) igraphmodule_Graph_successors,
   METH_VARARGS | METH_KEYWORDS,
   "successors(vertex)\n\n"
   "Returns the successors of a given vertex.\n\n"
   "Equivalent to calling the L{Graph.neighbors} method with type=L{OUT}."},

  {"predecessors", (PyCFunction) igraphmodule_Graph_predecessors,
   METH_VARARGS | METH_KEYWORDS,
   "predecessors(vertex)\n\n"
   "Returns the predecessors of a given vertex.\n\n"
   "Equivalent to calling the L{Graph.neighbors} method with type=L{IN}."},

  /* interface to igraph_get_eid */
  {"get_eid", (PyCFunction) igraphmodule_Graph_get_eid,
   METH_VARARGS | METH_KEYWORDS,
   "get_eid(v1, v2, directed=True, error=True)\n\n"
   "Returns the edge ID of an arbitrary edge between vertices v1 and v2\n\n"
   "@param v1: the first vertex ID\n"
   "@param v2: the second vertex ID\n"
   "@param directed: whether edge directions should be considered in\n"
   "  directed graphs. The default is C{True}. Ignored for undirected\n"
   "  graphs.\n"
   "@param error: if C{True}, an exception will be raised when the\n"
   "  given edge does not exist. If C{False}, -1 will be returned in\n"
   "  that case.\n"
   "@return: the edge ID of an arbitrary edge between vertices v1 and v2\n"},

  /* interface to igraph_get_eids */
  {"get_eids", (PyCFunction) igraphmodule_Graph_get_eids,
   METH_VARARGS | METH_KEYWORDS,
   "get_eids(pairs=None, path=None, directed=True, error=True)\n\n"
   "Returns the edge IDs of some edges between some vertices.\n\n"
   "This method can operate in two different modes, depending on which\n"
   "of the keyword arguments C{pairs} and C{path} are given.\n\n"
   "The method does not consider multiple edges; if there are multiple\n"
   "edges between a pair of vertices, only the ID of one of the edges\n"
   "is returned.\n\n"
   "@param pairs: a list of integer pairs. Each integer pair is considered\n"
   "  as a source-target vertex pair; the corresponding edge is looked up\n"
   "  in the graph and the edge ID is returned for each pair.\n"
   "@param path: a list of vertex IDs. The list is considered as a\n"
   "  continuous path from the first vertex to the last, passing\n"
   "  through the intermediate vertices. The corresponding edge IDs\n"
   "  between the first and the second, the second and the third and\n"
   "  so on are looked up in the graph and the edge IDs are returned.\n"
   "  If both C{path} and C{pairs} are given, the two lists are\n"
   "  concatenated.\n"
   "@param directed: whether edge directions should be considered in\n"
   "  directed graphs. The default is C{True}. Ignored for undirected\n"
   "  graphs.\n"
   "@param error: if C{True}, an exception will be raised if a given\n"
   "  edge does not exist. If C{False}, -1 will be returned in\n"
   "  that case.\n"
   "@return: the edge IDs in a list\n"},

  /* interface to igraph_incident */
  {"incident", (PyCFunction) igraphmodule_Graph_incident,
   METH_VARARGS | METH_KEYWORDS,
   "incident(vertex, mode=OUT)\n\n"
   "Returns the edges a given vertex is incident on.\n\n"
   "@param vertex: a vertex ID\n"
   "@param mode: whether to return only predecessors (L{OUT}),\n"
   "  successors (L{IN}) or both (L{ALL}). Ignored for undirected\n"
   "  graphs."},

  //////////////////////
  // GRAPH GENERATORS //
  //////////////////////

  /* interface to igraph_adjacency */
  {"Adjacency", (PyCFunction) igraphmodule_Graph_Adjacency,
   METH_CLASS | METH_VARARGS | METH_KEYWORDS,
   "Adjacency(matrix, mode=ADJ_DIRECTED)\n\n"
   "Generates a graph from its adjacency matrix.\n\n"
   "@param matrix: the adjacency matrix\n"
   "@param mode: the mode to be used. Possible values are:\n"
   "\n"
   "  - C{ADJ_DIRECTED} - the graph will be directed and a matrix\n"
   "    element gives the number of edges between two vertex.\n"
   "  - C{ADJ_UNDIRECTED} - alias to C{ADJ_MAX} for convenience.\n"
   "  - C{ADJ_MAX}   - undirected graph will be created and the number of\n"
   "    edges between vertex M{i} and M{j} is M{max(A(i,j), A(j,i))}\n"
   "  - C{ADJ_MIN}   - like C{ADJ_MAX}, but with M{min(A(i,j), A(j,i))}\n"
   "  - C{ADJ_PLUS}  - like C{ADJ_MAX}, but with M{A(i,j) + A(j,i)}\n"
   "  - C{ADJ_UPPER} - undirected graph with the upper right triangle of\n"
   "    the matrix (including the diagonal)\n"
   "  - C{ADJ_LOWER} - undirected graph with the lower left triangle of\n"
   "    the matrix (including the diagonal)\n"
   "\n"
   "  These values can also be given as strings without the C{ADJ} prefix.\n"
   },

  /* interface to igraph_asymmetric_preference_game */
  {"Asymmetric_Preference",
   (PyCFunction) igraphmodule_Graph_Asymmetric_Preference,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Asymmetric_Preference(n, type_dist_matrix, pref_matrix, attribute=None, loops=False)\n\n"
   "Generates a graph based on asymmetric vertex types and connection probabilities.\n\n"
   "This is the asymmetric variant of L{Graph.Preference}.\n"
   "A given number of vertices are generated. Every vertex is assigned to an\n"
   "\"incoming\" and an \"outgoing\" vertex type according to the given joint\n"
   "type probabilities. Finally, every vertex pair is evaluated and a\n"
   "directed edge is created between them with a probability depending on\n"
   "the \"outgoing\" type of the source vertex and the \"incoming\" type of\n"
   "the target vertex.\n\n"
   "@param n: the number of vertices in the graph\n"
   "@param type_dist_matrix: matrix giving the joint distribution of vertex\n"
   "  types\n"
   "@param pref_matrix: matrix giving the connection probabilities for\n"
   "  different vertex types.\n"
   "@param attribute: the vertex attribute name used to store the vertex\n"
   "  types. If C{None}, vertex types are not stored.\n"
   "@param loops: whether loop edges are allowed.\n"},

  // interface to igraph_atlas
  {"Atlas", (PyCFunction) igraphmodule_Graph_Atlas,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Atlas(idx)\n\n"
   "Generates a graph from the Graph Atlas.\n\n"
   "@param idx: The index of the graph to be generated.\n"
   "  Indices start from zero, graphs are listed:\n\n"
   "    1. in increasing order of number of vertices;\n"
   "    2. for a fixed number of vertices, in increasing order of the\n"
   "       number of edges;\n"
   "    3. for fixed numbers of vertices and edges, in increasing order\n"
   "       of the degree sequence, for example 111223 < 112222;\n"
   "    4. for fixed degree sequence, in increasing number of automorphisms.\n\n"
   "@newfield ref: Reference\n"
   "@ref: I{An Atlas of Graphs} by Ronald C. Read and Robin J. Wilson,\n"
   "  Oxford University Press, 1998."},

  // interface to igraph_barabasi_game
  {"Barabasi", (PyCFunction) igraphmodule_Graph_Barabasi,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Barabasi(n, m, outpref=False, directed=False, power=1,\n"
   "         zero_appeal=1, implementation=\"psumtree\", start_from=None)\n\n"
   "Generates a graph based on the Barabasi-Albert model.\n\n"
   "@param n: the number of vertices\n"
   "@param m: either the number of outgoing edges generated for\n"
   "  each vertex or a list containing the number of outgoing\n"
   "  edges for each vertex explicitly.\n"
   "@param outpref: C{True} if the out-degree of a given vertex\n"
   "  should also increase its citation probability (as well as\n"
   "  its in-degree), but it defaults to C{False}.\n"
   "@param directed: C{True} if the generated graph should be\n"
   "  directed (default: C{False}).\n"
   "@param power: the power constant of the nonlinear model.\n"
   "  It can be omitted, and in this case the usual linear model\n"
   "  will be used.\n"
   "@param zero_appeal: the attractivity of vertices with degree\n"
   "  zero.\n\n"
   "@param implementation: the algorithm to use to generate the\n"
   "  network. Possible values are:\n\n"
   "    - C{\"bag\"}: the algorithm that was the default in\n"
   "      igraph before 0.6. It works by putting the ids of the\n"
   "      vertices into a bag (multiset) exactly as many times\n"
   "      as their in-degree, plus once more. The required number\n"
   "      of cited vertices are then drawn from the bag with\n"
   "      replacement. It works only for I{power}=1 and\n"
   "      I{zero_appeal}=1.\n\n"
   "    - C{\"psumtree\"}: this algorithm uses a partial prefix-sum\n"
   "      tree to generate the graph. It does not generate multiple\n"
   "      edges and it works for any values of I{power} and\n"
   "      I{zero_appeal}.\n\n"
   "    - C{\"psumtree_multiple\"}: similar to C{\"psumtree\"}, but\n"
   "      it will generate multiple edges as well. igraph before\n"
   "      0.6 used this algorithm for I{power}s other than 1.\n\n"
   "@param start_from: if given and not C{None}, this must be another\n"
   "      L{Graph} object. igraph will use this graph as a starting\n"
   "      point for the preferential attachment model.\n\n"
   "@newfield ref: Reference\n"
   "@ref: Barabasi, A-L and Albert, R. 1999. Emergence of scaling\n"
   "  in random networks. Science, 286 509-512."},

  /* interface to igraph_create_bipartite */
  {"_Bipartite", (PyCFunction) igraphmodule_Graph_Bipartite,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "_Bipartite(types, edges, directed=False)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.Bipartite()\n\n"},

  /* interface to igraph_de_bruijn */
  {"De_Bruijn", (PyCFunction) igraphmodule_Graph_De_Bruijn,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "De_Bruijn(m, n)\n\n"
   "Generates a de Bruijn graph with parameters (m, n)\n\n"
   "A de Bruijn graph represents relationships between strings. An alphabet\n"
   "of M{m} letters are used and strings of length M{n} are considered.\n"
   "A vertex corresponds to every possible string and there is a directed edge\n"
   "from vertex M{v} to vertex M{w} if the string of M{v} can be transformed into\n"
   "the string of M{w} by removing its first letter and appending a letter to it.\n"
   "\n"
   "Please note that the graph will have M{m^n} vertices and even more edges,\n" 
   "so probably you don't want to supply too big numbers for M{m} and M{n}.\n\n"
   "@param m: the size of the alphabet\n"
   "@param n: the length of the strings\n"
  },

  // interface to igraph_establishment_game
  {"Establishment", (PyCFunction) igraphmodule_Graph_Establishment,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Establishment(n, k, type_dist, pref_matrix, directed=False)\n\n"
   "Generates a graph based on a simple growing model with vertex types.\n\n"
   "A single vertex is added at each time step. This new vertex tries to\n"
   "connect to k vertices in the graph. The probability that such a\n"
   "connection is realized depends on the types of the vertices involved.\n"
   "\n"
   "@param n: the number of vertices in the graph\n"
   "@param k: the number of connections tried in each step\n"
   "@param type_dist: list giving the distribution of vertex types\n"
   "@param pref_matrix: matrix (list of lists) giving the connection\n"
   "  probabilities for different vertex types\n"
   "@param directed: whether to generate a directed graph.\n"},

  // interface to igraph_erdos_renyi_game
  {"Erdos_Renyi", (PyCFunction) igraphmodule_Graph_Erdos_Renyi,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Erdos_Renyi(n, p, m, directed=False, loops=False)\n\n"
   "Generates a graph based on the Erdos-Renyi model.\n\n"
   "@param n: the number of vertices.\n"
   "@param p: the probability of edges. If given, C{m} must be missing.\n"
   "@param m: the number of edges. If given, C{p} must be missing.\n"
   "@param directed: whether to generate a directed graph.\n"
   "@param loops: whether self-loops are allowed.\n"},

  /* interface to igraph_famous */
	{"Famous", (PyCFunction) igraphmodule_Graph_Famous,
	 METH_VARARGS | METH_CLASS | METH_KEYWORDS,
	 "Famous(name)\n\n"
	 "Generates a famous graph based on its name.\n\n"
	 "Several famous graphs are known to C{igraph} including (but not limited to)\n"
	 "the Chvatal graph, the Petersen graph or the Tutte graph. This method\n"
	 "generates one of them based on its name (case insensitive). See the\n"
	 "documentation of the C interface of C{igraph} for the names available:\n"
	 "U{http://igraph.sourceforge.net/doc/html/igraph_famous.html}.\n\n"
	 "@param name: the name of the graph to be generated.\n"
	},

  /* interface to igraph_forest_fire_game */
  {"Forest_Fire", (PyCFunction) igraphmodule_Graph_Forest_Fire,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Forest_Fire(n, fw_prob, bw_factor=0.0, ambs=1, directed=False)\n\n"
   "Generates a graph based on the forest fire model\n\n"
   "The forest fire model is a growin graph model. In every time step, a new\n"
   "vertex is added to the graph. The new vertex chooses an ambassador (or\n"
   "more than one if M{ambs>1}) and starts a simulated forest fire at its\n"
   "ambassador(s). The fire spreads through the edges. The spreading probability\n"
   "along an edge is given by M{fw_prob}. The fire may also spread backwards\n"
   "on an edge by probability M{fw_prob * bw_factor}. When the fire ended, the\n"
   "newly added vertex connects to the vertices ``burned'' in the previous\n"
   "fire.\n\n"
   "@param n: the number of vertices in the graph\n"
   "@param fw_prob: forward burning probability\n"
   "@param bw_factor: ratio of backward and forward burning probability\n"
   "@param ambs: number of ambassadors chosen in each step\n"
   "@param directed: whether the graph will be directed\n"
  },

  /* interface to igraph_full_citation */
  {"Full_Citation", (PyCFunction) igraphmodule_Graph_Full_Citation,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Full_Citation(n, directed=False)\n\n"
   "Generates a full citation graph\n\n"
   "A full citation graph is a graph where the vertices are indexed from 0 to\n"
   "M{n-1} and vertex M{i} has a directed edge towards all vertices with an\n"
   "index less than M{i}.\n\n"
   "@param n: the number of vertices.\n"
   "@param directed: whether to generate a directed graph.\n"},

  /* interface to igraph_full */
  {"Full", (PyCFunction) igraphmodule_Graph_Full,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Full(n, directed=False, loops=False)\n\n"
   "Generates a full graph (directed or undirected, with or without loops).\n\n"
   "@param n: the number of vertices.\n"
   "@param directed: whether to generate a directed graph.\n"
   "@param loops: whether self-loops are allowed.\n"},

  /* interface to igraph_full_bipartite */
  {"_Full_Bipartite", (PyCFunction) igraphmodule_Graph_Full_Bipartite,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "_Full_Bipartite(n1, n2, directed=False, loops=False)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.Full_Bipartite()\n\n"},

  /* interface to igraph_grg_game */
  {"_GRG", (PyCFunction) igraphmodule_Graph_GRG,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "_GRG(n, radius, torus=False)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.GRG()\n\n"},

  /* interface to igraph_growing_random_game */
  {"Growing_Random", (PyCFunction) igraphmodule_Graph_Growing_Random,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Growing_Random(n, m, directed=False, citation=False)\n\n"
   "Generates a growing random graph.\n\n"
   "@param n: The number of vertices in the graph\n"
   "@param m: The number of edges to add in each step (after adding a new vertex)\n"
   "@param directed: whether the graph should be directed.\n"
   "@param citation: whether the new edges should originate from the most\n"
   "   recently added vertex.\n"},

  /* interface to igraph_incidence */
  {"_Incidence", (PyCFunction) igraphmodule_Graph_Incidence,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "_Incidence(matrix, directed=False, mode=ALL, multiple=False)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.Incidence()\n\n"},

  /* interface to igraph_kautz */
  {"Kautz", (PyCFunction) igraphmodule_Graph_Kautz,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Kautz(m, n)\n\n"
   "Generates a Kautz graph with parameters (m, n)\n\n"
   "A Kautz graph is a labeled graph, vertices are labeled by strings\n"
   "of length M{n+1} above an alphabet with M{m+1} letters, with\n"
   "the restriction that every two consecutive letters in the string\n"
   "must be different. There is a directed edge from a vertex M{v} to\n"
   "another vertex M{w} if it is possible to transform the string of\n"
   "M{v} into the string of M{w} by removing the first letter and\n"
   "appending a letter to it.\n\n"
   "@param m: the size of the alphabet minus one\n"
   "@param n: the length of the strings minus one\n"
   },

  /* interface to igraph_k_regular */
  {"K_Regular", (PyCFunction) igraphmodule_Graph_K_Regular,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "K_Regular(n, k, directed=False, multiple=False)\n\n"
   "Generates a k-regular random graph\n\n"
   "A k-regular random graph is a random graph where each vertex has degree k.\n"
   "If the graph is directed, both the in-degree and the out-degree of each vertex\n"
   "will be k.\n\n"
   "@param n: The number of vertices in the graph\n\n"
   "@param k: The degree of each vertex if the graph is undirected, or the in-degree\n"
   "  and out-degree of each vertex if the graph is directed\n"
   "@param directed: whether the graph should be directed.\n"
   "@param multiple: whether it is allowed to create multiple edges.\n"
  },

  /* interface to igraph_preference_game */
  {"Preference", (PyCFunction) igraphmodule_Graph_Preference,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Preference(n, type_dist, pref_matrix, attribute=None, directed=False, loops=False)\n\n"
   "Generates a graph based on vertex types and connection probabilities.\n\n"
   "This is practically the nongrowing variant of L{Graph.Establishment}.\n"
   "A given number of vertices are generated. Every vertex is assigned to a\n"
   "vertex type according to the given type probabilities. Finally, every\n"
   "vertex pair is evaluated and an edge is created between them with a\n"
   "probability depending on the types of the vertices involved.\n\n"
   "@param n: the number of vertices in the graph\n"
   "@param type_dist: list giving the distribution of vertex types\n"
   "@param pref_matrix: matrix giving the connection probabilities for\n"
   "  different vertex types.\n"
   "@param attribute: the vertex attribute name used to store the vertex\n"
   "  types. If C{None}, vertex types are not stored.\n"
   "@param directed: whether to generate a directed graph.\n"
   "@param loops: whether loop edges are allowed.\n"},

  /* interface to igraph_recent_degree_game */
  {"Recent_Degree", (PyCFunction) igraphmodule_Graph_Recent_Degree,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Recent_Degree(n, m, window, outpref=False, directed=False, power=1)\n\n"
   "Generates a graph based on a stochastic model where the probability\n"
   "of an edge gaining a new node is proportional to the edges gained in\n"
   "a given time window.\n\n"
   "@param n: the number of vertices\n"
   "@param m: either the number of outgoing edges generated for\n"
   "  each vertex or a list containing the number of outgoing\n"
   "  edges for each vertex explicitly.\n"
   "@param window: size of the window in time steps\n"
   "@param outpref: C{True} if the out-degree of a given vertex\n"
   "  should also increase its citation probability (as well as\n"
   "  its in-degree), but it defaults to C{False}.\n"
   "@param directed: C{True} if the generated graph should be\n"
   "  directed (default: C{False}).\n"
   "@param power: the power constant of the nonlinear model.\n"
   "  It can be omitted, and in this case the usual linear model\n"
   "  will be used.\n"},

  // interface to igraph_star
  {"Star", (PyCFunction) igraphmodule_Graph_Star,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Star(n, mode=\"undirected\", center=0)\n\n"
   "Generates a star graph.\n\n"
   "@param n: the number of vertices in the graph\n"
   "@param mode: Gives the type of the star graph to create. Should be\n"
   "  either \"in\", \"out\", \"mutual\" or \"undirected\"\n"
   "@param center: Vertex ID for the central vertex in the star.\n"},

  // interface to igraph_lattice
  {"Lattice", (PyCFunction) igraphmodule_Graph_Lattice,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Lattice(dim, nei=1, directed=False, mutual=True, circular=True)\n\n"
   "Generates a regular lattice.\n\n"
   "@param dim: list with the dimensions of the lattice\n"
   "@param nei: value giving the distance (number of steps) within which\n"
   "   two vertices will be connected.\n"
   "@param directed: whether to create a directed graph.\n"
   "@param mutual: whether to create all connections as mutual\n"
   "    in case of a directed graph.\n"
   "@param circular: whether the generated lattice is periodic.\n"},

  /* interface to igraph_lcf */
	{"LCF", (PyCFunction) igraphmodule_Graph_LCF,
	 METH_VARARGS | METH_CLASS | METH_KEYWORDS,
	 "LCF(n, shifts, repeats)\n\n"
	 "Generates a graph from LCF notation.\n\n"
	 "LCF is short for Lederberg-Coxeter-Frucht, it is a concise notation\n"
	 "for 3-regular Hamiltonian graphs. It consists of three parameters,\n"
	 "the number of vertices in the graph, a list of shifts giving\n"
	 "additional edges to a cycle backbone and another integer giving how\n"
	 "many times the shifts should be performed. See\n"
	 "U{http://mathworld.wolfram.com/LCFNotation.html} for details.\n\n" 
	 "@param n: the number of vertices\n"
	 "@param shifts: the shifts in a list or tuple\n"
	 "@param repeats: the number of repeats\n"
	},
	 
  // interface to igraph_ring
  {"Ring", (PyCFunction) igraphmodule_Graph_Ring,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Ring(n, directed=False, mutual=False, circular=True)\n\n"
   "Generates a ring graph.\n\n"
   "@param n: the number of vertices in the ring\n"
   "@param directed: whether to create a directed ring.\n"
   "@param mutual: whether to create mutual edges in a directed ring.\n"
   "@param circular: whether to create a closed ring.\n"},

  /* interface to igraph_static_fitness_game */
  {"Static_Fitness", (PyCFunction) igraphmodule_Graph_Static_Fitness,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Static_Fitness(m, fitness_out, fitness_in=None, loops=False, multiple=False)\n\n"
   "Generates a non-growing graph with edge probabilities proportional to node\n"
   "fitnesses.\n\n"
   "The algorithm randomly selects vertex pairs and connects them until the given\n"
   "number of edges are created. Each vertex is selected with a probability\n"
   "proportional to its fitness; for directed graphs, a vertex is selected as a\n"
   "source proportional to its out-fitness and as a target proportional to its\n"
   "in-fitness.\n\n"
   "@param m: the number of edges in the graph\n"
   "@param fitness_out: a numeric vector with non-negative entries, one for each\n"
   "  vertex. These values represent the fitness scores (out-fitness scores for\n"
   "  directed graphs). I{fitness} is an alias of this keyword argument.\n"
   "@param fitness_in: a numeric vector with non-negative entries, one for each\n"
   "  vertex. These values represent the in-fitness scores for directed graphs.\n"
   "  For undirected graphs, this argument must be C{None}.\n"
   "@param loops: whether loop edges are allowed.\n"
   "@param multiple: whether multiple edges are allowed.\n"
   "@return: a directed or undirected graph with the prescribed power-law\n"
   "  degree distributions.\n"
  },

  /* interface to igraph_static_power_law_game */
  {"Static_Power_Law", (PyCFunction) igraphmodule_Graph_Static_Power_Law,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Static_Power_Law(n, m, exponent_out, exponent_in=-1, loops=False,\n"
   "    multiple=False, finite_size_correction=True)\n\n"
   "Generates a non-growing graph with prescribed power-law degree distributions.\n\n"
   "@param n: the number of vertices in the graph\n"
   "@param m: the number of edges in the graph\n"
   "@param exponent_out: the exponent of the out-degree distribution, which\n"
   "  must be between 2 and infinity (inclusive). When I{exponent_in} is\n"
   "  not given or negative, the graph will be undirected and this parameter\n"
   "  specifies the degree distribution. I{exponent} is an alias to this\n"
   "  keyword argument.\n"
   "@param exponent_in: the exponent of the in-degree distribution, which\n"
   "  must be between 2 and infinity (inclusive) It can also be negative, in\n"
   "  which case an undirected graph will be generated.\n"
   "@param loops: whether loop edges are allowed.\n"
   "@param multiple: whether multiple edges are allowed.\n"
   "@param finite_size_correction: whether to apply a finite-size correction\n"
   "  to the generated fitness values for exponents less than 3. See the\n"
   "  paper of Cho et al for more details.\n"
   "@return: a directed or undirected graph with the prescribed power-law\n"
   "  degree distributions.\n"
   "\n"
   "@newfield ref: Reference\n"
   "@ref: Goh K-I, Kahng B, Kim D: Universal behaviour of load distribution\n"
   "  in scale-free networks. Phys Rev Lett 87(27):278701, 2001.\n"
   "@ref: Cho YS, Kim JS, Park J, Kahng B, Kim D: Percolation transitions in\n"
   "  scale-free networks under the Achlioptas process. Phys Rev Lett\n"
   "  103:135702, 2009.\n"
  },

  // interface to igraph_tree
  {"Tree", (PyCFunction) igraphmodule_Graph_Tree,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Tree(n, children, type=TREE_UNDIRECTED)\n\n"
   "Generates a tree in which almost all vertices have the same number of children.\n\n"
   "@param n: the number of vertices in the graph\n"
   "@param children: the number of children of a vertex in the graph\n"
   "@param type: determines whether the tree should be directed, and if\n"
   "  this is the case, also its orientation. Must be one of\n"
   "  C{TREE_IN}, C{TREE_OUT} and C{TREE_UNDIRECTED}.\n"},

  /* interface to igraph_degree_sequence_game */
  {"Degree_Sequence", (PyCFunction) igraphmodule_Graph_Degree_Sequence,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Degree_Sequence(out, in=None, method=\"simple\")\n\n"
   "Generates a graph with a given degree sequence.\n\n"
   "@param out: the out-degree sequence for a directed graph. If the\n"
   "  in-degree sequence is omitted, the generated graph\n"
   "  will be undirected, so this will be the in-degree\n"
   "  sequence as well\n"
   "@param in: the in-degree sequence for a directed graph.\n"
   "  If omitted, the generated graph will be undirected.\n"
   "@param method: the generation method to be used. One of the following:\n"
   "  \n"
   "    - C{\"simple\"} -- simple generator that sometimes generates\n"
   "      loop edges and multiple edges. The generated graph is not\n"
   "      guaranteed to be connected.\n"
   "    - C{\"no_multiple\"} -- similar to C{\"simple\"} but avoids the\n"
   "      generation of multiple and loop edges at the expense of increased\n"
   "      time complexity. The method will re-start the generation every time\n"
   "      it gets stuck in a configuration where it is not possible to insert\n"
   "      any more edges without creating loops or multiple edges, and there\n"
   "      is no upper bound on the number of iterations, but it will succeed\n"
   "      eventually if the input degree sequence is graphical and throw an\n"
   "      exception if the input degree sequence is not graphical.\n"
   "    - C{\"vl\"} -- a more sophisticated generator that can sample\n"
   "      undirected, connected simple graphs uniformly. It uses\n"
   "      Monte-Carlo methods to randomize the graphs.\n"
   "      This generator should be favoured if undirected and connected\n"
   "      graphs are to be generated and execution time is not a concern.\n"
   "      igraph uses the original implementation of Fabien Viger; see the\n"
   "      following URL and the paper cited on it for the details of the\n"
   "      algorithm: U{http://www-rp.lip6.fr/~latapy/FV/generation.html}.\n"
  },

  /* interface to igraph_isoclass_create */
  {"Isoclass", (PyCFunction) igraphmodule_Graph_Isoclass,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Isoclass(n, class, directed=False)\n\n"
   "Generates a graph with a given isomorphy class.\n\n"
   "@param n: the number of vertices in the graph (3 or 4)\n"
   "@param class: the isomorphy class\n"
   "@param directed: whether the graph should be directed.\n"},

  /* interface to igraph_watts_strogatz_game */
  {"Watts_Strogatz", (PyCFunction) igraphmodule_Graph_Watts_Strogatz,
   METH_VARARGS | METH_CLASS | METH_KEYWORDS,
   "Watts_Strogatz(dim, size, nei, p, loops=False, multiple=False)\n\n"
   "@param dim: the dimension of the lattice\n"
   "@param size: the size of the lattice along all dimensions\n"
   "@param nei: value giving the distance (number of steps) within which\n"
   "   two vertices will be connected.\n"
   "@param p: rewiring probability\n\n"
   "@param loops: specifies whether loop edges are allowed\n"
   "@param multiple: specifies whether multiple edges are allowed\n"
   "@see: L{Lattice()}, L{rewire()}, L{rewire_edges()} if more flexibility is\n"
   "  needed\n"
   "@newfield ref: Reference\n"
   "@ref: Duncan J Watts and Steven H Strogatz: I{Collective dynamics of\n"
   "  small world networks}, Nature 393, 440-442, 1998\n"},

  /* interface to igraph_weighted_adjacency */
  {"Weighted_Adjacency", (PyCFunction) igraphmodule_Graph_Weighted_Adjacency,
   METH_CLASS | METH_VARARGS | METH_KEYWORDS,
   "Weighted_Adjacency(matrix, mode=ADJ_DIRECTED, attr=\"weight\", loops=True)\n\n"
   "Generates a graph from its adjacency matrix.\n\n"
   "@param matrix: the adjacency matrix\n"
   "@param mode: the mode to be used. Possible values are:\n"
   "\n"
   "  - C{ADJ_DIRECTED} - the graph will be directed and a matrix\n"
   "    element gives the number of edges between two vertex.\n"
   "  - C{ADJ_UNDIRECTED} - alias to C{ADJ_MAX} for convenience.\n"
   "  - C{ADJ_MAX}   - undirected graph will be created and the number of\n"
   "    edges between vertex M{i} and M{j} is M{max(A(i,j), A(j,i))}\n"
   "  - C{ADJ_MIN}   - like C{ADJ_MAX}, but with M{min(A(i,j), A(j,i))}\n"
   "  - C{ADJ_PLUS}  - like C{ADJ_MAX}, but with M{A(i,j) + A(j,i)}\n"
   "  - C{ADJ_UPPER} - undirected graph with the upper right triangle of\n"
   "    the matrix (including the diagonal)\n"
   "  - C{ADJ_LOWER} - undirected graph with the lower left triangle of\n"
   "    the matrix (including the diagonal)\n"
   "\n"
   "  These values can also be given as strings without the C{ADJ} prefix.\n"
   "@param attr: the name of the edge attribute that stores the edge\n"
   "  weights.\n"
   "@param loops: whether to include loop edges. When C{False}, the diagonal\n"
   "  of the adjacency matrix will be ignored.\n"
  },

  /////////////////////////////////////
  // STRUCTURAL PROPERTIES OF GRAPHS //
  /////////////////////////////////////

  // interface to igraph_are_connected
  {"are_connected", (PyCFunction) igraphmodule_Graph_are_connected,
   METH_VARARGS | METH_KEYWORDS,
   "are_connected(v1, v2)\n\n"
   "Decides whether two given vertices are directly connected.\n\n"
   "@param v1: the first vertex\n"
   "@param v2: the second vertex\n"
   "@return: C{True} if there exists an edge from v1 to v2, C{False}\n"
   "  otherwise.\n"},

  /* interface to igraph_articulation_points */
  {"articulation_points", (PyCFunction)igraphmodule_Graph_articulation_points,
   METH_NOARGS,
   "articulation_points()\n\n"
   "Returns the list of articulation points in the graph.\n\n"
   "A vertex is an articulation point if its removal increases the number of\n"
   "connected components in the graph.\n"
  },

  /* interface to igraph_assortativity */
  {"assortativity", (PyCFunction)igraphmodule_Graph_assortativity,
   METH_VARARGS | METH_KEYWORDS,
   "assortativity(types1, types2=None, directed=True)\n\n"
   "Returns the assortativity of the graph based on numeric properties\n"
   "of the vertices.\n\n"
   "This coefficient is basically the correlation between the actual\n"
   "connectivity patterns of the vertices and the pattern expected from the\n"
   "disribution of the vertex types.\n\n"
   "See equation (21) in Newman MEJ: Mixing patterns in networks, Phys Rev E\n"
   "67:026126 (2003) for the proper definition. The actual calculation is\n"
   "performed using equation (26) in the same paper for directed graphs, and\n"
   "equation (4) in Newman MEJ: Assortative mixing in networks, Phys Rev Lett\n"
   "89:208701 (2002) for undirected graphs.\n\n"
   "@param types1: vertex types in a list or the name of a vertex attribute\n"
   "  holding vertex types. Types are ideally denoted by numeric values.\n"
   "@param types2: in directed assortativity calculations, each vertex can\n"
   "  have an out-type and an in-type. In this case, I{types1} contains the\n"
   "  out-types and this parameter contains the in-types in a list or the\n"
   "  name of a vertex attribute. If C{None}, it is assumed to be equal\n"
   "  to I{types1}.\n\n"
   "@param directed: whether to consider edge directions or not.\n"
   "@return: the assortativity coefficient\n\n"
   "@newfield ref: Reference\n"
   "@ref: Newman MEJ: Mixing patterns in networks, Phys Rev E 67:026126, 2003.\n"
   "@ref: Newman MEJ: Assortative mixing in networks, Phys Rev Lett 89:208701,\n"
   "  2002.\n"
   "@see: L{assortativity_degree()} when the types are the vertex degrees\n"
  },

  /* interface to igraph_assortativity_degree */
  {"assortativity_degree", (PyCFunction)igraphmodule_Graph_assortativity_degree,
   METH_VARARGS | METH_KEYWORDS,
   "assortativity_degree(directed=True)\n\n"
   "Returns the assortativity of a graph based on vertex degrees.\n\n"
   "See L{assortativity()} for the details. L{assortativity_degree()} simply\n"
   "calls L{assortativity()} with the vertex degrees as types.\n\n"
   "@param directed: whether to consider edge directions for directed graphs\n"
   "  or not. This argument is ignored for undirected graphs.\n"
   "@return: the assortativity coefficient\n\n"
   "@see: L{assortativity()}\n"
  },

  /* interface to igraph_assortativity_nominal */
  {"assortativity_nominal", (PyCFunction)igraphmodule_Graph_assortativity_nominal,
   METH_VARARGS | METH_KEYWORDS,
   "assortativity_nominal(types, directed=True)\n\n"
   "Returns the assortativity of the graph based on vertex categories.\n\n"
   "Assuming that the vertices belong to different categories, this\n"
   "function calculates the assortativity coefficient, which specifies\n"
   "the extent to which the connections stay within categories. The\n"
   "assortativity coefficient is one if all the connections stay within\n"
   "categories and minus one if all the connections join vertices of\n"
   "different categories. For a randomly connected network, it is\n"
   "asymptotically zero.\n\n"
   "See equation (2) in Newman MEJ: Mixing patterns in networks, Phys Rev E\n"
   "67:026126 (2003) for the proper definition.\n\n"
   "@param types: vertex types in a list or the name of a vertex attribute\n"
   "  holding vertex types. Types should be denoted by numeric values.\n"
   "@param directed: whether to consider edge directions or not.\n"
   "@return: the assortativity coefficient\n\n"
   "@newfield ref: Reference\n"
   "@ref: Newman MEJ: Mixing patterns in networks, Phys Rev E 67:026126, 2003.\n"
  },

  /* interface to igraph_average_path_length */
  {"average_path_length",
   (PyCFunction) igraphmodule_Graph_average_path_length,
   METH_VARARGS | METH_KEYWORDS,
   "average_path_length(directed=True, unconn=True)\n\n"
   "Calculates the average path length in a graph.\n\n"
   "@param directed: whether to consider directed paths in case of a\n"
   "  directed graph. Ignored for undirected graphs.\n"
   "@param unconn: what to do when the graph is unconnected. If C{True},\n"
   "  the average of the geodesic lengths in the components is\n"
   "  calculated. Otherwise for all unconnected vertex pairs,\n"
   "  a path length equal to the number of vertices is used.\n"
   "@return: the average path length in the graph\n"},

  /* interface to igraph_authority_score */
  {"authority_score", (PyCFunction)igraphmodule_Graph_authority_score,
   METH_VARARGS | METH_KEYWORDS,
   "authority_score(weights=None, scale=True, arpack_options=None, return_eigenvalue=False)\n\n"
   "Calculates Kleinberg's authority score for the vertices of the graph\n\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@param scale: whether to normalize the scores so that the largest one\n"
   "  is 1.\n"
   "@param arpack_options: an L{ARPACKOptions} object used to fine-tune\n"
   "  the ARPACK eigenvector calculation. If omitted, the module-level\n"
   "  variable called C{arpack_options} is used.\n"
   "@param return_eigenvalue: whether to return the largest eigenvalue\n"
   "@return: the authority scores in a list and optionally the largest eigenvalue\n"
   "  as a second member of a tuple\n\n"
   "@see: hub_score()\n"
  },

  /* interface to igraph_betweenness[_estimate] */
  {"betweenness", (PyCFunction) igraphmodule_Graph_betweenness,
   METH_VARARGS | METH_KEYWORDS,
   "betweenness(vertices=None, directed=True, cutoff=None, weights=None, nobigint=True)\n\n"
   "Calculates or estimates the betweenness of vertices in a graph.\n\n"
   "Keyword arguments:\n"
   "@param vertices: the vertices for which the betweennesses must be returned.\n"
   "  If C{None}, assumes all of the vertices in the graph.\n"
   "@param directed: whether to consider directed paths.\n"
   "@param cutoff: if it is an integer, only paths less than or equal to this\n"
   "  length are considered, effectively resulting in an estimation of the\n"
   "  betweenness for the given vertices. If C{None}, the exact betweenness is\n"
   "  returned.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@param nobigint: if C{True}, igraph uses the longest available integer\n"
   "  type on the current platform to count shortest paths. For some large\n"
   "  networks that have a specific structure, the counters may overflow.\n"
   "  To prevent this, use C{nobigint=False}, which forces igraph to use\n"
   "  arbitrary precision integers at the expense of increased computation\n"
   "  time.\n"
   "@return: the (possibly estimated) betweenness of the given vertices in a list\n"},

  /* interface to biconnected_components */
  {"biconnected_components", (PyCFunction) igraphmodule_Graph_biconnected_components,
   METH_VARARGS | METH_KEYWORDS,
   "biconnected_components(return_articulation_points=True)\n\n"
   "Calculates the biconnected components of the graph.\n\n"
   "@param return_articulation_points: whether to return the articulation\n"
   "  points as well\n"
   "@return: a list of lists containing edge indices making up spanning trees\n"
   "  of the biconnected components (one spanning tree for each component)\n"
   "  and optionally the list of articulation points"
  },

  /* interface to igraph_bipartite_projection */
  {"bipartite_projection", (PyCFunction) igraphmodule_Graph_bipartite_projection,
   METH_VARARGS | METH_KEYWORDS,
   "bipartite_projection(types, multiplicity=True, probe1=-1)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.bipartite_projection()\n"},

  /* interface to igraph_bipartite_projection_size */
  {"bipartite_projection_size", (PyCFunction) igraphmodule_Graph_bipartite_projection_size,
   METH_VARARGS | METH_KEYWORDS,
   "bipartite_projection_size(types)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.bipartite_projection_size()\n"},

  /* interface to igraph_closeness */
  {"closeness", (PyCFunction) igraphmodule_Graph_closeness,
   METH_VARARGS | METH_KEYWORDS,
   "closeness(vertices=None, mode=ALL, cutoff=None, weights=None)\n\n"
   "Calculates the closeness centralities of given vertices in a graph.\n\n"
   "The closeness centerality of a vertex measures how easily other\n"
   "vertices can be reached from it (or the other way: how easily it\n"
   "can be reached from the other vertices). It is defined as the\n"
   "number of the number of vertices minus one divided by the sum of\n"
   "the lengths of all geodesics from/to the given vertex.\n\n"
   "If the graph is not connected, and there is no path between two\n"
   "vertices, the number of vertices is used instead the length of\n"
   "the geodesic. This is always longer than the longest possible\n"
   "geodesic.\n\n"
   "@param vertices: the vertices for which the closenesses must\n"
   "  be returned. If C{None}, uses all of the vertices in the graph.\n"
   "@param mode: must be one of L{IN}, L{OUT} and L{ALL}. L{IN} means\n"
   "  that the length of the incoming paths, L{OUT} means that the\n"
   "  length of the outgoing paths must be calculated. L{ALL} means\n"
   "  that both of them must be calculated.\n"
   "@param cutoff: if it is an integer, only paths less than or equal to this\n"
   "  length are considered, effectively resulting in an estimation of the\n"
   "  closeness for the given vertices (which is always an underestimation of the\n"
   "  real closeness, since some vertex pairs will appear as disconnected even\n"
   "  though they are connected).. If C{None}, the exact closeness is\n"
   "  returned.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@return: the calculated closenesses in a list\n"},

  /* interface to igraph_clusters */
  {"clusters", (PyCFunction) igraphmodule_Graph_clusters,
   METH_VARARGS | METH_KEYWORDS,
   "clusters(mode=STRONG)\n\n"
   "Calculates the (strong or weak) clusters for a given graph.\n\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a L{VertexClustering} object.\n"
   "  It is advised to use that.\n"
   "@param mode: must be either C{STRONG} or C{WEAK}, depending on\n"
   "  the clusters being sought. Optional, defaults to C{STRONG}.\n"
   "@return: the component index for every node in the graph.\n"},
  {"copy", (PyCFunction) igraphmodule_Graph_copy,
   METH_NOARGS,
   "copy()\n\n" "Creates an exact deep copy of the graph."},
  {"decompose", (PyCFunction) igraphmodule_Graph_decompose,
   METH_VARARGS | METH_KEYWORDS,
   "decompose(mode=STRONG, maxcompno=None, minelements=1)\n\n"
   "Decomposes the graph into subgraphs.\n\n"
   "@param mode: must be either STRONG or WEAK, depending on the\n"
   "  clusters being sought.\n"
   "@param maxcompno: maximum number of components to return.\n"
   "  C{None} means all possible components.\n"
   "@param minelements: minimum number of vertices in a component.\n"
   "  By setting this to 2, isolated vertices are not returned\n"
   "  as separate components.\n"
   "@return: a list of the subgraphs. Every returned subgraph is a\n"
   "  copy of the original.\n"},
  /* interface to igraph_contract_vertices */
  {"contract_vertices", (PyCFunction) igraphmodule_Graph_contract_vertices,
   METH_VARARGS | METH_KEYWORDS,
   "contract_vertices(mapping, combine_attrs=None)\n\n"
   "Contracts some vertices in the graph, i.e. replaces groups of vertices\n"
   "with single vertices. Edges are not affected.\n\n"
   "@param mapping: numeric vector which gives the mapping between old and\n"
   "  new vertex IDs. Vertices having the same new vertex ID in this vector\n"
   "  will be remapped into a single new vertex. It is safe to pass the\n"
   "  membership vector of a L{VertexClustering} object here.\n"
   "@param combine_attrs: specifies how to combine the attributes of\n"
   "  the vertices being collapsed into a single one. If it is C{None},\n"
   "  all the attributes will be lost. If it is a function, the\n"
   "  attributes of the vertices will be collected and passed on to\n"
   "  that function which will return the new attribute value that has to\n"
   "  be assigned to the single collapsed vertex. It can also be one of\n"
   "  the following string constants which define built-in collapsing\n"
   "  functions: C{sum}, C{prod}, C{mean}, C{median}, C{max}, C{min},\n"
   "  C{first}, C{last}, C{random}. You can also specify different\n"
   "  combination functions for different attributes by passing a dict\n"
   "  here which maps attribute names to functions. See\n"
   "  L{Graph.simplify()} for more details.\n"
   "@return: C{None}.\n"
   "@see: L{Graph.simplify()}\n"
  },
  /* interface to igraph_constraint */
  {"constraint", (PyCFunction) igraphmodule_Graph_constraint,
   METH_VARARGS | METH_KEYWORDS,
   "constraint(vertices=None, weights=None)\n\n"
   "Calculates Burt's constraint scores for given vertices in a graph.\n\n"
   "Burt's constraint is higher if ego has less, or mutually stronger\n"
   "related (i.e. more redundant) contacts. Burt's measure of\n"
   "constraint, C[i], of vertex i's ego network V[i], is defined for\n"
   "directed and valued graphs as follows:\n\n"
   "C[i] = sum( sum( (p[i,q] p[q,j])^2, q in V[i], q != i,j ), j in V[], j != i)\n\n"
   "for a graph of order (ie. number od vertices) N, where proportional\n"
   "tie strengths are defined as follows:\n\n"
   "p[i,j]=(a[i,j]+a[j,i]) / sum(a[i,k]+a[k,i], k in V[i], k != i),\n"
   "a[i,j] are elements of A and the latter being the graph adjacency matrix.\n\n"
   "For isolated vertices, constraint is undefined.\n\n"
   "@param vertices: the vertices to be analysed or C{None} for all vertices.\n"
   "@param weights: weights associated to the edges. Can be an attribute name\n"
   "  as well. If C{None}, every edge will have the same weight.\n"
   "@return: constraint scores for all given vertices in a matrix."},

  /* interface to igraph_density */
  {"density", (PyCFunction) igraphmodule_Graph_density,
   METH_VARARGS | METH_KEYWORDS,
   "density(loops=False)\n\n"
   "Calculates the density of the graph.\n\n"
   "@param loops: whether to take loops into consideration. If C{True},\n"
   "  the algorithm assumes that there might be some loops in the graph\n"
   "  and calculates the density accordingly. If C{False}, the algorithm\n"
   "  assumes that there can't be any loops.\n"
   "@return: the reciprocity of the graph."},

  /* interfaces to igraph_diameter */
  {"diameter", (PyCFunction) igraphmodule_Graph_diameter,
   METH_VARARGS | METH_KEYWORDS,
   "diameter(directed=True, unconn=True, weights=None)\n\n"
   "Calculates the diameter of the graph.\n\n"
   "@param directed: whether to consider directed paths.\n"
   "@param unconn: if C{True} and the graph is unconnected, the\n"
   "  longest geodesic within a component will be returned. If\n"
   "  C{False} and the graph is unconnected, the result is the\n"
   "  number of vertices if there are no weights or infinity\n"
   "  if there are weights.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@return: the diameter"},
  {"get_diameter", (PyCFunction) igraphmodule_Graph_get_diameter,
   METH_VARARGS | METH_KEYWORDS,
   "get_diameter(directed=True, unconn=True, weights=None)\n\n"
   "Returns a path with the actual diameter of the graph.\n\n"
   "If there are many shortest paths with the length of the diameter,\n"
   "it returns the first one it founds.\n\n"
   "@param directed: whether to consider directed paths.\n"
   "@param unconn: if C{True} and the graph is unconnected, the\n"
   "  longest geodesic within a component will be returned. If\n"
   "  C{False} and the graph is unconnected, the result is the\n"
   "  number of vertices if there are no weights or infinity\n"
   "  if there are weights.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@return: the vertices in the path in order."},
  {"farthest_points", (PyCFunction) igraphmodule_Graph_farthest_points,
   METH_VARARGS | METH_KEYWORDS,
   "farthest_points(directed=True, unconn=True, weights=None)\n\n"
   "Returns two vertex IDs whose distance equals the actual diameter\n"
   "of the graph.\n\n"
   "If there are many shortest paths with the length of the diameter,\n"
   "it returns the first one it found.\n\n"
   "@param directed: whether to consider directed paths.\n"
   "@param unconn: if C{True} and the graph is unconnected, the\n"
   "  longest geodesic within a component will be returned. If\n"
   "  C{False} and the graph is unconnected, the result contains the\n"
   "  number of vertices if there are no weights or infinity\n"
   "  if there are weights.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@return: a triplet containing the two vertex IDs and their distance.\n"
   "  The IDs are C{None} if the graph is unconnected and C{unconn}\n"
   "  is C{False}."},

  /* interface to igraph_diversity */
  {"diversity", (PyCFunction) igraphmodule_Graph_diversity,
   METH_VARARGS | METH_KEYWORDS,
   "diversity(vertices=None, weights=None)\n\n"
   "Calculates the structural diversity index of the vertices.\n\n"
   "The structural diversity index of a vertex is simply the (normalized)\n"
   "Shannon entropy of the weights of the edges incident on the vertex.\n\n"
   "The measure is defined for undirected graphs only; edge directions are\n"
   "ignored.\n\n"
   "@param vertices: the vertices for which the diversity indices must\n"
   "  be returned. If C{None}, uses all of the vertices in the graph.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@return: the calculated diversity indices in a list, or a single number if\n"
   "  a single vertex was supplied.\n"
   "@newfield ref: Reference\n"
   "@ref: Eagle N, Macy M and Claxton R: Network diversity and economic\n"
   "  development, Science 328, 1029--1031, 2010."
  },

  /* interface to igraph_eccentricity */
  {"eccentricity", (PyCFunction) igraphmodule_Graph_eccentricity,
   METH_VARARGS | METH_KEYWORDS,
   "eccentricity(vertices=None, mode=ALL)\n\n"
   "Calculates the eccentricities of given vertices in a graph.\n\n"
   "The eccentricity of a vertex is calculated by measuring the\n"
   "shortest distance from (or to) the vertex, to (or from) all other\n"
   "vertices in the graph, and taking the maximum.\n\n"
   "@param vertices: the vertices for which the eccentricity scores must\n"
   "  be returned. If C{None}, uses all of the vertices in the graph.\n"
   "@param mode: must be one of L{IN}, L{OUT} and L{ALL}. L{IN} means\n"
   "  that edge directions are followed; C{OUT} means that edge directions\n"
   "  are followed the opposite direction; C{ALL} means that directions are\n"
   "  ignored. The argument has no effect for undirected graphs.\n"
   "@return: the calculated eccentricities in a list, or a single number if\n"
   "  a single vertex was supplied.\n"},

  /* interface to igraph_edge_betweenness[_estimate] */
  {"edge_betweenness", (PyCFunction) igraphmodule_Graph_edge_betweenness,
   METH_VARARGS | METH_KEYWORDS,
   "edge_betweenness(directed=True, cutoff=None, weights=None)\n\n"
   "Calculates or estimates the edge betweennesses in a graph.\n\n"
   "@param directed: whether to consider directed paths.\n"
   "@param cutoff: if it is an integer, only paths less than or equal to this\n"
   "  length are considered, effectively resulting in an estimation of the\n"
   "  betweenness values. If C{None}, the exact betweennesses are\n"
   "  returned.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@return: a list with the (exact or estimated) edge betweennesses of all\n"
   "  edges.\n"},

  /* interface to igraph_[st_]edge_connectivity */
  {"edge_connectivity", (PyCFunction) igraphmodule_Graph_edge_connectivity,
   METH_VARARGS | METH_KEYWORDS,
   "edge_connectivity(source=-1, target=-1, checks=True)\n\n"
   "Calculates the edge connectivity of the graph or between some vertices.\n\n"
   "The edge connectivity between two given vertices is the number of edges\n"
   "that have to be removed in order to disconnect the two vertices into two\n"
   "separate components. This is also the number of edge disjoint directed\n"
   "paths between the vertices. The edge connectivity of the graph is the minimal\n"
   "edge connectivity over all vertex pairs.\n\n"
   "This method calculates the edge connectivity of a given vertex pair if both\n"
   "the source and target vertices are given. If none of them is given (or they\n"
   "are both negative), the overall edge connectivity is returned.\n\n"
   "@param source: the source vertex involved in the calculation.\n"
   "@param target: the target vertex involved in the calculation.\n"
   "@param checks: if the whole graph connectivity is calculated and this is\n"
   "  C{True}, igraph performs some basic checks before calculation. If the\n"
   "  graph is not strongly connected, then the connectivity is obviously\n"
   "  zero. If the minimum degree is one, then the connectivity is\n"
   "  also one. These simple checks are much faster than checking the entire\n"
   "  graph, therefore it is advised to set this to C{True}. The parameter\n"
   "  is ignored if the connectivity between two given vertices is computed.\n"
   "@return: the edge connectivity\n"
  },

  /* interface to igraph_eigenvector_centrality */
  {"eigenvector_centrality",
   (PyCFunction) igraphmodule_Graph_eigenvector_centrality,
   METH_VARARGS | METH_KEYWORDS,
   "eigenvector_centrality(directed=True, scale=True, weights=None, return_eigenvalue=False, arpack_options=None)\n\n"
   "Calculates the eigenvector centralities of the vertices in a graph.\n\n"
   "@param directed: whether to consider edge directions in a directed\n"
   "  graph. Ignored for undirected graphs.\n"
   "@param scale: whether to normalize the centralities so the largest\n"
   "  one will always be 1.\n"
   "@param weights: edge weights given as a list or an edge attribute. If\n"
   "  C{None}, all edges have equal weight.\n"
   "@param return_eigenvalue: whether to return the actual largest\n"
   "  eigenvalue along with the centralities\n"
   "@param arpack_options: an L{ARPACKOptions} object that can be used\n"
   "  to fine-tune the calculation. If it is omitted, the module-level\n"
   "  variable called C{arpack_options} is used.\n"
   "@return: the eigenvector centralities in a list and optionally the\n"
   "  largest eigenvalue (as a second member of a tuple)"
  },

  /* interface to igraph_feedback_arc_set */
  {"feedback_arc_set", (PyCFunction) igraphmodule_Graph_feedback_arc_set,
   METH_VARARGS | METH_KEYWORDS,
   "feedback_arc_set(weights=None, method=\"eades\")\n\n"
   "Calculates an approximately or exactly minimal feedback arc set.\n\n"
   "A feedback arc set is a set of edges whose removal makes the graph acyclic.\n"
   "Since this is always possible by removing all the edges, we are in general\n"
   "interested in removing the smallest possible number of edges, or an edge set\n"
   "with as small total weight as possible. This method calculates one such edge\n"
   "set. Note that the task is trivial for an undirected graph as it is enough\n"
   "to find a spanning tree and then remove all the edges not in the spanning\n"
   "tree. Of course it is more complicated for directed graphs.\n\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name. When given, the algorithm will strive to\n"
   "  remove lightweight edges in order to minimize the total weight of the\n"
   "  feedback arc set.\n"
   "@param method: the algorithm to use. C{\"eades\"} uses the greedy cycle\n"
   "  breaking heuristic of Eades, Lin and Smyth, which is linear in the number\n"
   "  of edges but not necessarily optimal; however, it guarantees that the\n"
   "  number of edges to be removed is smaller than |E|/2 - |V|/6. C{\"ip\"} uses\n"
   "  an integer programming formulation which is guaranteed to yield an optimal\n"
   "  result, but is too slow for large graphs.\n"
   "@return: the IDs of the edges to be removed, in a list.\n\n"
   "@newfield ref: Reference\n"
   "@ref: Eades P, Lin X and Smyth WF: A fast and effective heuristic for the\n"
   "  feedback arc set problem. In: Proc Inf Process Lett 319-323, 1993.\n"
  },

  // interface to igraph_get_shortest_paths
  {"get_shortest_paths", (PyCFunction) igraphmodule_Graph_get_shortest_paths,
   METH_VARARGS | METH_KEYWORDS,
   "get_shortest_paths(v, to=None, weights=None, mode=OUT, output=\"vpath\")\n\n"
   "Calculates the shortest paths from/to a given node in a graph.\n\n"
   "@param v: the source/destination for the calculated paths\n"
   "@param to: a vertex selector describing the destination/source for\n"
   "  the calculated paths. This can be a single vertex ID, a list of\n"
   "  vertex IDs, a single vertex name, a list of vertex names or a\n"
   "  L{VertexSeq} object. C{None} means all the vertices.\n"
   "@param weights: edge weights in a list or the name of an edge attribute\n"
   "  holding edge weights. If C{None}, all edges are assumed to have\n"
   "  equal weight.\n"
   "@param mode: the directionality of the paths. L{IN} means to\n"
   "  calculate incoming paths, L{OUT} means to calculate outgoing\n"
   "  paths, L{ALL} means to calculate both ones.\n"
   "@param output: determines what should be returned. If this is\n"
   "  C{\"vpath\"}, a list of vertex IDs will be returned, one path\n"
   "  for each target vertex. For unconnected graphs, some of the list\n"
   "  elements may be empty. Note that in case of mode=L{IN}, the vertices\n"
   "  in a path are returned in reversed order. If C{output=\"epath\"},\n"
   "  edge IDs are returned instead of vertex IDs.\n"
   "@return: see the documentation of the C{output} parameter.\n"},

  // interface to igraph_get_all_shortest_paths
  {"get_all_shortest_paths",
   (PyCFunction) igraphmodule_Graph_get_all_shortest_paths,
   METH_VARARGS | METH_KEYWORDS,
   "get_all_shortest_paths(v, to=None, weights=None, mode=OUT)\n\n"
   "Calculates all of the shortest paths from/to a given node in a graph.\n\n"
   "@param v: the source/destination for the calculated paths\n"
   "@param to: a vertex selector describing the destination/source for\n"
   "  the calculated paths. This can be a single vertex ID, a list of\n"
   "  vertex IDs, a single vertex name, a list of vertex names or a\n"
   "  L{VertexSeq} object. C{None} means all the vertices.\n"
   "@param weights: edge weights in a list or the name of an edge attribute\n"
   "  holding edge weights. If C{None}, all edges are assumed to have\n"
   "  equal weight.\n"
   "@param mode: the directionality of the paths. L{IN} means to calculate\n"
   "  incoming paths, L{OUT} means to calculate outgoing paths,\n"
   "  L{ALL} means to calculate both ones.\n"
   "@return: all of the shortest path from the given node to every other\n"
   "reachable node in the graph in a list. Note that in case of mode=L{IN},\n"
   "the vertices in a path are returned in reversed order!"},

  /* interface to igraph_girth */
  {"girth", (PyCFunction)igraphmodule_Graph_girth,
   METH_VARARGS | METH_KEYWORDS,
   "girth(return_shortest_circle=False)\n\n"
   "Returns the girth of the graph.\n\n"
   "The girth of a graph is the length of the shortest circle in it.\n\n"
   "@param return_shortest_circle: whether to return one of the shortest\n"
   "  circles found in the graph.\n"
   "@return: the length of the shortest circle or (if C{return_shortest_circle})\n"
   "  is true, the shortest circle itself as a list\n"
  },

  /* interface to igraph_convergence_degree */
  {"convergence_degree", (PyCFunction)igraphmodule_Graph_convergence_degree,
    METH_NOARGS,
    "convergence_degree()\n\n"
	"Undocumented (yet)."
  },

  /* interface to igraph_convergence_field_size */
  {"convergence_field_size", (PyCFunction)igraphmodule_Graph_convergence_field_size,
    METH_NOARGS,
    "convergence_field_size()\n\n"
	"Undocumented (yet)."
  },

  /* interface to igraph_hub_score */
  {"hub_score", (PyCFunction)igraphmodule_Graph_hub_score,
   METH_VARARGS | METH_KEYWORDS,
   "hub_score(weights=None, scale=True, arpack_options=None, return_eigenvalue=False)\n\n"
   "Calculates Kleinberg's hub score for the vertices of the graph\n\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@param scale: whether to normalize the scores so that the largest one\n"
   "  is 1.\n"
   "@param arpack_options: an L{ARPACKOptions} object used to fine-tune\n"
   "  the ARPACK eigenvector calculation. If omitted, the module-level\n"
   "  variable called C{arpack_options} is used.\n"
   "@param return_eigenvalue: whether to return the largest eigenvalue\n"
   "@return: the hub scores in a list and optionally the largest eigenvalue\n"
   "  as a second member of a tuple\n\n"
   "@see: authority_score()\n"
  },

  /* interface to igraph_induced_subgraph */
  {"induced_subgraph", (PyCFunction) igraphmodule_Graph_induced_subgraph,
   METH_VARARGS | METH_KEYWORDS,
   "induced_subgraph(vertices, implementation=\"auto\")\n\n"
   "Returns a subgraph spanned by the given vertices.\n\n"
   "@param vertices: a list containing the vertex IDs which\n"
   "  should be included in the result.\n"
   "@param implementation: the implementation to use when constructing\n"
   "  the new subgraph. igraph includes two implementations at the\n"
   "  moment. C{\"copy_and_delete\"} copies the original graph and\n"
   "  removes those vertices that are not in the given set. This is more\n"
   "  efficient if the size of the subgraph is comparable to the original\n"
   "  graph. The other implementation (C{\"create_from_scratch\"})\n"
   "  constructs the result graph from scratch and then copies the\n"
   "  attributes accordingly. This is a better solution if the subgraph\n"
   "  is relatively small, compared to the original graph. C{\"auto\"}\n"
   "  selects between the two implementations automatically, based on\n"
   "  the ratio of the size of the subgraph and the size of the original\n"
   "  graph.\n"
   "@return: the subgraph\n"},

  /* interface to igraph_is_bipartite */
  {"is_bipartite", (PyCFunction) igraphmodule_Graph_is_bipartite,
   METH_VARARGS | METH_KEYWORDS,
   "is_bipartite(return_types=False)\n\n"
   "Decides whether the graph is bipartite or not.\n\n"
   "Vertices of a bipartite graph can be partitioned into two groups A\n"
   "and B in a way that all edges go between the two groups.\n\n"
   "@param return_types: if C{False}, the method will simply\n"
   "  return C{True} or C{False} depending on whether the graph is\n"
   "  bipartite or not. If C{True}, the actual group assignments\n"
   "  are also returned as a list of boolean values. (Note that\n"
   "  the group assignment is not unique, especially if the graph\n"
   "  consists of multiple components, since the assignments of\n"
   "  components are independent from each other).\n"
   "@return: C{True} if the graph is bipartite, C{False} if not.\n"
   "  If C{return_types} is C{True}, the group assignment is also\n"
   "  returned.\n"
  },

  /* interface to igraph_avg_nearest_neighbor_degree */
  {"knn", (PyCFunction) igraphmodule_Graph_knn,
   METH_VARARGS | METH_KEYWORDS,
   "knn(vids=None, weights=None)\n\n"
   "Calculates the average degree of the neighbors for each vertex, and\n"
   "the same quantity as the function of vertex degree.\n\n"
   "@param vids: the vertices for which the calculation is performed.\n"
   "  C{None} means all vertices.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name. If this is given, the vertex strength\n"
   "  will be used instead of the vertex degree in the calculations, but\n"
   "  the \"ordinary\" vertex degree will be used for the second (degree-\n"
   "  dependent) list in the result.\n"
   "@return: two lists in a tuple. The first list contains the average\n"
   "  degree of neighbors for each vertex, the second contains the average\n"
   "  degree of neighbors as a function of vertex degree. The zeroth element\n"
   "  of this list corresponds to vertices of degree 1.\n"
  },

  /* interface to igraph_is_connected */
  {"is_connected", (PyCFunction) igraphmodule_Graph_is_connected,
   METH_VARARGS | METH_KEYWORDS,
   "is_connected(mode=STRONG)\n\n"
   "Decides whether the graph is connected.\n\n"
   "@param mode: whether we should calculate strong or weak connectivity.\n"
   "@return: C{True} if the graph is connected, C{False} otherwise.\n"},

  /* interface to igraph_linegraph */
  {"linegraph", (PyCFunction) igraphmodule_Graph_linegraph,
   METH_VARARGS | METH_KEYWORDS,
   "linegraph()\n\n"
   "Returns the line graph of the graph.\n\n"
   "The line graph M{L(G)} of an undirected graph is defined as follows:\n"
   "M{L(G)} has one vertex for each edge in G and two vertices in M{L(G)}\n"
   "are connected iff their corresponding edges in the original graph\n"
   "share an end point.\n\n"
   "The line graph of a directed graph is slightly different: two vertices\n"
   "are connected by a directed edge iff the target of the first vertex's\n"
   "corresponding edge is the same as the source of the second vertex's\n"
   "corresponding edge.\n"
  },

  /* interface to igraph_maxdegree */
  {"maxdegree", (PyCFunction) igraphmodule_Graph_maxdegree,
   METH_VARARGS | METH_KEYWORDS,
   "maxdegree(vertices=None, mode=ALL, loops=False)\n\n"
   "Returns the maximum degree of a vertex set in the graph.\n\n"
   "This method accepts a single vertex ID or a list of vertex IDs as a\n"
   "parameter, and returns the degree of the given vertices (in the\n"
   "form of a single integer or a list, depending on the input\n"
   "parameter).\n"
   "\n"
   "@param vertices: a single vertex ID or a list of vertex IDs, or\n"
   "  C{None} meaning all the vertices in the graph.\n"
   "@param mode: the type of degree to be returned (L{OUT} for\n"
   "  out-degrees, L{IN} IN for in-degrees or L{ALL} for the sum of\n"
   "  them).\n"
   "@param loops: whether self-loops should be counted.\n"},

  /* interface to igraph_neighborhood */
  {"neighborhood", (PyCFunction) igraphmodule_Graph_neighborhood,
   METH_VARARGS | METH_KEYWORDS,
   "neighborhood(vertices=None, order=1, mode=ALL)\n\n"
   "For each vertex specified by I{vertices}, returns the\n"
   "vertices reachable from that vertex in at most I{order} steps.\n\n"
   "@param vertices: a single vertex ID or a list of vertex IDs, or\n"
   "  C{None} meaning all the vertices in the graph.\n"
   "@param mode: specifies how to take into account the direction of\n"
   "  the edges if a directed graph is analyzed. C{\"out\"} means that\n"
   "  only the outgoing edges are followed, so all vertices reachable\n"
   "  from the source vertex in at most I{order} steps are counted.\n"
   "  C{\"in\"} means that only the incoming edges are followed (in\n"
   "  reverse direction of course), so all vertices from which the source\n"
   "  vertex is reachable in at most I{order} steps are counted. C{\"all\"}\n"
   "  treats directed edges as undirected.\n"
   "@return: a single list specifying the neighborhood if I{vertices}\n"
   "  was an integer specifying a single vertex index, or a list of lists\n"
   "  if I{vertices} was a list or C{None}.\n"
  },

  /* interface to igraph_neighborhood_size */
  {"neighborhood_size", (PyCFunction) igraphmodule_Graph_neighborhood_size,
   METH_VARARGS | METH_KEYWORDS,
   "neighborhood_size(vertices=None, order=1, mode=ALL)\n\n"
   "For each vertex specified by I{vertices}, returns the number of\n"
   "vertices reachable from that vertex in at most I{order} steps.\n\n"
   "@param vertices: a single vertex ID or a list of vertex IDs, or\n"
   "  C{None} meaning all the vertices in the graph.\n"
   "@param mode: specifies how to take into account the direction of\n"
   "  the edges if a directed graph is analyzed. C{\"out\"} means that\n"
   "  only the outgoing edges are followed, so all vertices reachable\n"
   "  from the source vertex in at most I{order} steps are counted.\n"
   "  C{\"in\"} means that only the incoming edges are followed (in\n"
   "  reverse direction of course), so all vertices from which the source\n"
   "  vertex is reachable in at most I{order} steps are counted. C{\"all\"}\n"
   "  treats directed edges as undirected.\n"
   "@return: a single number specifying the neighborhood size if I{vertices}\n"
   "  was an integer specifying a single vertex index, or a list of sizes\n"
   "  if I{vertices} was a list or C{None}.\n"
  },

  /* interface to igraph_personalized_pagerank */
  {"personalized_pagerank", (PyCFunction) igraphmodule_Graph_personalized_pagerank,
   METH_VARARGS | METH_KEYWORDS,
   "personalized_pagerank(vertices=None, directed=True, damping=0.85, "
   "reset=None, reset_vertices=None, weights=None, arpack_options=None)\n\n"
   "Calculates the personalized PageRank values of a graph.\n\n"
   "The personalized PageRank calculation is similar to the PageRank\n"
   "calculation, but the random walk is reset to a non-uniform distribution\n"
   "over the vertices in every step with probability M{1-damping} instead of a\n"
   "uniform distribution.\n\n"
   "@param vertices: the indices of the vertices being queried.\n"
   "  C{None} means all of the vertices.\n"
   "@param directed: whether to consider directed paths.\n"
   "@param damping: the damping factor.\n"
   "  M{1-damping} is the PageRank value for vertices with no\n"
   "  incoming links.\n"
   "@param reset: the distribution over the vertices to be used when resetting\n"
   "  the random walk. Can be a sequence, an iterable or a vertex attribute\n"
   "  name as long as they return a list of floats whose length is equal to\n"
   "  the number of vertices. If C{None}, a uniform distribution is assumed,\n"
   "  which makes the method equivalent to the original PageRank algorithm.\n"
   "@param reset_vertices: an alternative way to specify the distribution\n"
   "  over the vertices to be used when resetting the random walk. Simply\n"
   "  supply a list of vertex IDs here, or a L{VertexSeq} or a L{Vertex}.\n"
   "  Resetting will take place using a uniform distribution over the specified\n"
   "  vertices.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@param arpack_options: an L{ARPACKOptions} object used to fine-tune\n"
   "  the ARPACK eigenvector calculation. If omitted, the module-level\n"
   "  variable called C{arpack_options} is used.\n"
   "@return: a list with the personalized PageRank values of the specified\n"
   "  vertices.\n"},

  /* interface to igraph_path_length_hist */
  {"path_length_hist", (PyCFunction) igraphmodule_Graph_path_length_hist,
   METH_VARARGS | METH_KEYWORDS,
   "path_length_hist(directed=True)\n\n"
   "Calculates the path length histogram of the graph\n"
   "@attention: this function is wrapped in a more convenient syntax in the\n"
   "  derived class L{Graph}. It is advised to use that instead of this version.\n\n"
   "@param directed: whether to consider directed paths\n"
   "@return: a tuple. The first item of the tuple is a list of path lengths,\n"
   "  the M{i}th element of the list contains the number of paths with length\n"
   "  M{i+1}. The second item contains the number of unconnected vertex pairs\n"
   "  as a float (since it might not fit into an integer)\n"
  },

  /* interface to igraph_permute_vertices */
  {"permute_vertices", (PyCFunction) igraphmodule_Graph_permute_vertices,
   METH_VARARGS | METH_KEYWORDS,
   "permute_vertices(permutation)\n\n"
   "Permutes the vertices of the graph according to the given permutation\n"
   "and returns the new graph.\n\n"
   "Vertex M{k} of the original graph will become vertex M{permutation[k]}\n"
   "in the new graph. No validity checks are performed on the permutation\n"
   "vector.\n\n"
   "@return: the new graph\n"
  },

  /* interfaces to igraph_radius */
  {"radius", (PyCFunction) igraphmodule_Graph_radius,
   METH_VARARGS | METH_KEYWORDS,
   "radius(mode=OUT)\n\n"
   "Calculates the radius of the graph.\n\n"
   "The radius of a graph is defined as the minimum eccentricity of\n"
   "its vertices (see L{eccentricity()}).\n"
   "@param mode: what kind of paths to consider for the calculation\n"
   "  in case of directed graphs. C{OUT} considers paths that follow\n"
   "  edge directions, C{IN} considers paths that follow the opposite\n"
   "  edge directions, C{ALL} ignores edge directions. The argument is\n"
   "  ignored for undirected graphs.\n"
   "@return: the radius\n"
   "@see: L{Graph.eccentricity()}"
  },

  /* interface to igraph_reciprocity */
  {"reciprocity", (PyCFunction) igraphmodule_Graph_reciprocity,
   METH_VARARGS | METH_KEYWORDS,
   "reciprocity(ignore_loops=True, mode=\"default\")\n\n"
   "Reciprocity defines the proportion of mutual connections in a\n"
   "directed graph. It is most commonly defined as the probability\n"
   "that the opposite counterpart of a directed edge is also included\n"
   "in the graph. This measure is calculated if C{mode} is C{\"default\"}.\n"
   "\n"
   "Prior to igraph 0.6, another measure was implemented, defined as\n"
   "the probability of mutual connection between a vertex pair if we\n"
   "know that there is a (possibly non-mutual) connection between them.\n"
   "In other words, (unordered) vertex pairs are classified into three\n"
   "groups: (1) disconnected, (2) non-reciprocally connected and (3)\n"
   "reciprocally connected. The result is the size of group (3), divided\n"
   "by the sum of sizes of groups (2) and (3). This measure is calculated\n"
   "if C{mode} is C{\"ratio\"}.\n"
   "\n"
   "@param ignore_loops: whether loop edges should be ignored.\n"
   "@param mode: the algorithm to use to calculate the reciprocity; see\n"
   "  above for more details.\n"
   "@return: the reciprocity of the graph\n"
  },

  /* interface to igraph_rewire */
  {"rewire", (PyCFunction) igraphmodule_Graph_rewire,
   METH_VARARGS | METH_KEYWORDS,
   "rewire(n=1000, mode=\"simple\")\n\n"
   "Randomly rewires the graph while preserving the degree distribution.\n\n"
   "Please note that the rewiring is done \"in-place\", so the original\n"
   "graph will be modified. If you want to preserve the original graph,\n"
   "use the L{copy} method before.\n\n"
   "@param n: the number of rewiring trials.\n"
   "@param mode: the rewiring algorithm to use. It can either be C{\"simple\"} or\n"
   "  C{\"loops\"}; the former does not create or destroy loop edges while the\n"
   "  latter does.\n"},

  /* interface to igraph_rewire_edges */
  {"rewire_edges", (PyCFunction) igraphmodule_Graph_rewire_edges,
   METH_VARARGS | METH_KEYWORDS,
   "rewire_edges(prob, loops=False, multiple=False)\n\n"
   "Rewires the edges of a graph with constant probability.\n\n"
   "Each endpoint of each edge of the graph will be rewired with a constant\n"
   "probability, given in the first argument.\n\n"
   "Please note that the rewiring is done \"in-place\", so the original\n"
   "graph will be modified. If you want to preserve the original graph,\n"
   "use the L{copy} method before.\n\n"
   "@param prob: rewiring probability\n"
   "@param loops: whether the algorithm is allowed to create loop edges\n"
   "@param multiple: whether the algorithm is allowed to create multiple\n"
   "  edges.\n"},

  /* interface to igraph_shortest_paths */
  {"shortest_paths", (PyCFunction) igraphmodule_Graph_shortest_paths,
   METH_VARARGS | METH_KEYWORDS,
   "shortest_paths(source=None, target=None, weights=None, mode=OUT)\n\n"
   "Calculates shortest path lengths for given vertices in a graph.\n\n"
   "The algorithm used for the calculations is selected automatically:\n"
   "a simple BFS is used for unweighted graphs, Dijkstra's algorithm is\n"
   "used when all the weights are positive. Otherwise, the Bellman-Ford\n"
   "algorithm is used if the number of requested source vertices is larger\n"
   "than 100 and Johnson's algorithm is used otherwise.\n\n"
   "@param source: a list containing the source vertex IDs which should be\n"
   "  included in the result. If C{None}, all vertices will be considered.\n"
   "@param target: a list containing the target vertex IDs which should be\n"
   "  included in the result. If C{None}, all vertices will be considered.\n"
   "@param weights: a list containing the edge weights. It can also be\n"
   "  an attribute name (edge weights are retrieved from the given\n"
   "  attribute) or C{None} (all edges have equal weight).\n"
   "@param mode: the type of shortest paths to be used for the\n"
   "  calculation in directed graphs. L{OUT} means only outgoing,\n"
   "  L{IN} means only incoming paths. L{ALL} means to consider\n"
   "  the directed graph as an undirected one.\n"
   "@return: the shortest path lengths for given vertices in a matrix\n"},

  /* interface to igraph_simplify */
  {"simplify", (PyCFunction) igraphmodule_Graph_simplify,
   METH_VARARGS | METH_KEYWORDS,
   "simplify(multiple=True, loops=True, combine_edges=None)\n\n"
   "Simplifies a graph by removing self-loops and/or multiple edges.\n\n"
   "\n"
   "For example, suppose you have a graph with an edge attribute named\n"
   "C{weight}. C{graph.simplify(combine_edges=max)} will take the\n"
   "maximum of the weights of multiple edges and assign that weight to\n"
   "the collapsed edge. C{graph.simplify(combine_edges=sum)} will\n"
   "take the sum of the weights. You can also write\n"
   "C{graph.simplify(combine_edges=dict(weight=\"sum\"))} or\n"
   "C{graph.simplify(combine_edges=dict(weight=sum))}, since\n"
   "C{sum} is recognised both as a Python built-in function and as\n"
   "a string constant.\n\n"
   "@param multiple: whether to remove multiple edges.\n"
   "@param loops: whether to remove loops.\n"
   "@param combine_edges: specifies how to combine the attributes of\n"
   "  multiple edges between the same pair of vertices into a single\n"
   "  attribute. If it is C{None}, only one of the edges will be kept\n"
   "  and all the attributes will be lost. If it is a function, the\n"
   "  attributes of multiple edges will be collected and passed on to\n"
   "  that function which will return the new attribute value that has to\n"
   "  be assigned to the single collapsed edge. It can also be one of\n"
   "  the following string constants:\n\n"
   "    - C{\"ignore\"}: all the edge attributes will be ignored.\n\n"
   "    - C{\"sum\"}: the sum of the edge attribute values will be used for\n"
   "      the new edge.\n\n"
   "    - C{\"product\"}: the product of the edge attribute values will be used for\n"
   "      the new edge.\n"
   "    - C{\"mean\"}: the mean of the edge attribute values will be used for\n"
   "      the new edge.\n\n"
   "    - C{\"median\"}: the median of the edge attribute values will be used for\n"
   "      the new edge.\n\n"
   "    - C{\"min\"}: the minimum of the edge attribute values will be used for\n"
   "      the new edge.\n\n"
   "    - C{\"max\"}: the maximum of the edge attribute values will be used for\n"
   "      the new edge.\n\n"
   "    - C{\"first\"}: the attribute value of the first edge in the collapsed set\n"
   "      will be used for the new edge.\n\n"
   "    - C{\"last\"}: the attribute value of the last edge in the collapsed set\n"
   "      will be used for the new edge.\n\n"
   "    - C{\"random\"}: a randomly selected value will be used for the new edge\n\n"
   "    - C{\"concat\"}: the attribute values will be concatenated for the new\n"
   "      edge.\n\n"
   "  You can also use a dict mapping edge attribute names to functions or\n"
   "  the above string constants if you want to make the behaviour of the\n"
   "  simplification process depend on the name of the attribute.\n"
   "  C{None} is a special key in this dict, its value will be used for all\n"
   "  the attributes not specified explicitly in the dictionary.\n"
  },

  /* interface to igraph_minimum_spanning_tree */
  {"_spanning_tree", (PyCFunction) igraphmodule_Graph_spanning_tree,
   METH_VARARGS | METH_KEYWORDS,
   "_spanning_tree(weights=None)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.spanning_tree()"},

  // interface to igraph_subcomponent
  {"subcomponent", (PyCFunction) igraphmodule_Graph_subcomponent,
   METH_VARARGS | METH_KEYWORDS,
   "subcomponent(v, mode=ALL)\n\n"
   "Determines the indices of vertices which are in the same component as a given vertex.\n\n"
   "@param v: the index of the vertex used as the source/destination\n"
   "@param mode: if equals to L{IN}, returns the vertex IDs from\n"
   "  where the given vertex can be reached. If equals to L{OUT},\n"
   "  returns the vertex IDs which are reachable from the given\n"
   "  vertex. If equals to L{ALL}, returns all vertices within the\n"
   "  same component as the given vertex, ignoring edge directions.\n"
   "  Note that this is not equal to calculating the union of the \n"
   "  results of L{IN} and L{OUT}.\n"
   "@return: the indices of vertices which are in the same component as a given vertex.\n"},

  /* interface to igraph_subgraph_edges */
  {"subgraph_edges", (PyCFunction) igraphmodule_Graph_subgraph_edges,
   METH_VARARGS | METH_KEYWORDS,
   "subgraph_edges(edges, delete_vertices=True)\n\n"
   "Returns a subgraph spanned by the given edges.\n\n"
   "@param edges: a list containing the edge IDs which should\n"
   "  be included in the result.\n"
   "@param delete_vertices: if C{True}, vertices not incident on\n"
   "  any of the specified edges will be deleted from the result.\n"
   "  If C{False}, all vertices will be kept.\n"
   "@return: the subgraph\n"},

  /* interface to igraph_topological_sorting */
  {"topological_sorting",
   (PyCFunction) igraphmodule_Graph_topological_sorting,
   METH_VARARGS | METH_KEYWORDS,
   "topological_sorting(mode=OUT)\n\n"
   "Calculates a possible topological sorting of the graph.\n\n"
   "Returns a partial sorting and issues a warning if the graph is not\n"
   "a directed acyclic graph.\n\n"
   "@param mode: if L{OUT}, vertices are returned according to the\n"
   "  forward topological order -- all vertices come before their\n"
   "  successors. If L{IN}, all vertices come before their ancestors.\n"
   "@return: a possible topological ordering as a list"},

  // interface to igraph_transitivity_undirected
  {"transitivity_undirected",
   (PyCFunction) igraphmodule_Graph_transitivity_undirected,
   METH_VARARGS | METH_KEYWORDS,
   "transitivity_undirected(mode=\"nan\")\n\n"
   "Calculates the global transitivity (clustering coefficient) of the\n"
   "graph.\n\n"
   "The transitivity measures the probability that two neighbors of a\n"
   "vertex are connected. More precisely, this is the ratio of the\n"
   "triangles and connected triplets in the graph. The result is a\n"
   "single real number. Directed graphs are considered as undirected\n"
   "ones.\n\n"
   "Note that this measure is different from the local transitivity\n"
   "measure (see L{transitivity_local_undirected()}) as it calculates\n"
   "a single value for the whole graph.\n\n"
   "@param mode: if C{TRANSITIVITY_ZERO} or C{\"zero\"}, the result will\n"
   "  be zero if the graph does not have any triplets. If C{\"nan\"} or\n"
   "  C{TRANSITIVITY_NAN}, the result will be C{NaN} (not a number).\n"
   "@return: the transitivity\n"
   "@see: L{transitivity_local_undirected()}, L{transitivity_avglocal_undirected()}\n"
   "@newfield ref: Reference\n"
   "@ref: S. Wasserman and K. Faust: I{Social Network Analysis: Methods and\n"
   "  Applications}. Cambridge: Cambridge University Press, 1994."
  },

  /* interface to igraph_transitivity_local_undirected and
   * igraph_transitivity_barrat */
  {"transitivity_local_undirected",
   (PyCFunction) igraphmodule_Graph_transitivity_local_undirected,
   METH_VARARGS | METH_KEYWORDS,
   "transitivity_local_undirected(vertices=None, mode=\"nan\", weights=None)\n\n"
   "Calculates the local transitivity (clustering coefficient) of the\n"
   "given vertices in the graph.\n\n"
   "The transitivity measures the probability that two neighbors of a\n"
   "vertex are connected. In case of the local transitivity, this\n"
   "probability is calculated separately for each vertex.\n\n"
   "Note that this measure is different from the global transitivity\n"
   "measure (see L{transitivity_undirected()}) as it calculates\n"
   "a transitivity value for each vertex individually.\n\n"
   "The traditional local transitivity measure applies for unweighted graphs\n"
   "only. When the C{weights} argument is given, this function calculates\n"
   "the weighted local transitivity proposed by Barrat et al (see references).\n\n"
   "@param vertices: a list containing the vertex IDs which should be\n"
   "  included in the result. C{None} means all of the vertices.\n"
   "@param mode: defines how to treat vertices with degree less than two.\n"
   "  If C{TRANSITIVITT_ZERO} or C{\"zero\"}, these vertices will have\n"
   "  zero transitivity. If C{TRANSITIVITY_NAN} or C{\"nan\"}, these\n"
   "  vertices will have C{NaN} (not a number) as their transitivity.\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@return: the transitivities for the given vertices in a list\n"
   "@see: L{transitivity_undirected()}, L{transitivity_avglocal_undirected()}\n"
   "@newfield ref: Reference\n"
   "@ref: Watts DJ and Strogatz S: I{Collective dynamics of small-world\n"
   "  networks}. Nature 393(6884):440-442, 1998.\n"
   "@ref: Barrat A, Barthelemy M, Pastor-Satorras R and Vespignani A:\n"
   "  I{The architecture of complex weighted networks}. PNAS 101, 3747 (2004).\n"
   "  U{http://arxiv.org/abs/cond-mat/0311416}."
  },

  /* interface to igraph_transitivity_avglocal_undirected */
  {"transitivity_avglocal_undirected",
   (PyCFunction) igraphmodule_Graph_transitivity_avglocal_undirected,
   METH_VARARGS | METH_KEYWORDS,
   "transitivity_avglocal_undirected(mode=\"nan\")\n\n"
   "Calculates the average of the vertex transitivities of the graph.\n\n"
   "The transitivity measures the probability that two neighbors of a\n"
   "vertex are connected. In case of the average local transitivity,\n"
   "this probability is calculated for each vertex and then the average\n"
   "is taken. Vertices with less than two neighbors require special\n"
   "treatment, they will either be left out from the calculation or\n"
   "they will be considered as having zero transitivity, depending on\n"
   "the I{mode} parameter.\n\n"
   "Note that this measure is different from the global transitivity measure\n"
   "(see L{transitivity_undirected()}) as it simply takes the average local\n"
   "transitivity across the whole network.\n\n"
   "@param mode: defines how to treat vertices with degree less than two.\n"
   "  If C{TRANSITIVITT_ZERO} or C{\"zero\"}, these vertices will have\n"
   "  zero transitivity. If C{TRANSITIVITY_NAN} or C{\"nan\"}, these\n"
   "  vertices will be excluded from the average.\n"
   "@see: L{transitivity_undirected()}, L{transitivity_local_undirected()}\n"
   "@newfield ref: Reference\n"
   "@ref: D. J. Watts and S. Strogatz: I{Collective dynamics of small-world\n"
   "  networks}. Nature 393(6884):440-442, 1998."
  },

  /* interface to igraph_unfold_tree */
  {"unfold_tree", (PyCFunction) igraphmodule_Graph_unfold_tree,
   METH_VARARGS | METH_KEYWORDS,
   "unfold_tree(sources=None, mode=OUT)\n\n"
   "Unfolds the graph using a BFS to a tree by duplicating vertices as necessary.\n\n"
   "@param sources: the source vertices to start the unfolding from. It should be a\n"
   "  list of vertex indices, preferably one vertex from each connected component.\n"
   "  You can use L{Graph.topological_sorting()} to determine a suitable set. A single\n"
   "  vertex index is also accepted.\n"
   "@param mode: which edges to follow during the BFS. C{OUT} follows outgoing edges,\n"
   "  C{IN} follows incoming edges, C{ALL} follows both. Ignored for undirected\n"
   "  graphs.\n"
   "@return: the unfolded tree graph and a mapping from the new vertex indices to the\n"
   "  old ones.\n"
  },

  /* interface to igraph_[st_]vertex_connectivity */
  {"vertex_connectivity", (PyCFunction) igraphmodule_Graph_vertex_connectivity,
   METH_VARARGS | METH_KEYWORDS,
   "vertex_connectivity(source=-1, target=-1, checks=True, neighbors=\"error\")\n\n"
   "Calculates the vertex connectivity of the graph or between some vertices.\n\n"
   "The vertex connectivity between two given vertices is the number of vertices\n"
   "that have to be removed in order to disconnect the two vertices into two\n"
   "separate components. This is also the number of vertex disjoint directed\n"
   "paths between the vertices (apart from the source and target vertices of\n"
   "course). The vertex connectivity of the graph is the minimal vertex\n"
   "connectivity over all vertex pairs.\n\n"
   "This method calculates the vertex connectivity of a given vertex pair if both\n"
   "the source and target vertices are given. If none of them is given (or they\n"
   "are both negative), the overall vertex connectivity is returned.\n\n"
   "@param source: the source vertex involved in the calculation.\n"
   "@param target: the target vertex involved in the calculation.\n"
   "@param checks: if the whole graph connectivity is calculated and this is\n"
   "  C{True}, igraph performs some basic checks before calculation. If the\n"
   "  graph is not strongly connected, then the connectivity is obviously\n"
   "  zero. If the minimum degree is one, then the connectivity is\n"
   "  also one. These simple checks are much faster than checking the entire\n"
   "  graph, therefore it is advised to set this to C{True}. The parameter\n"
   "  is ignored if the connectivity between two given vertices is computed.\n"
   "@param neighbors: tells igraph what to do when the two vertices are\n"
   "  connected. C{\"error\"} raises an exception, C{\"infinity\"} returns\n"
   "  infinity, C{\"ignore\"} ignores the edge.\n"
   "@return: the vertex connectivity\n"
  },

  /***********************/
  /* SIMILARITY MEASURES */
  /***********************/

  /* interface to igraph_bibcoupling */
  {"bibcoupling", (PyCFunction) igraphmodule_Graph_bibcoupling,
   METH_VARARGS | METH_KEYWORDS,
   "bibcoupling(vertices=None)\n\n"
   "Calculates bibliographic coupling scores for given vertices in a graph.\n\n"
   "@param vertices: the vertices to be analysed. If C{None}, all vertices\n"
   "  will be considered.\n"
   "@return: bibliographic coupling scores for all given vertices in a matrix."},
  /* interface to igraph_cocitation */
  {"cocitation", (PyCFunction) igraphmodule_Graph_cocitation,
   METH_VARARGS | METH_KEYWORDS,
   "cocitation(vertices=None)\n\n"
   "Calculates cocitation scores for given vertices in a graph.\n\n"
   "@param vertices: the vertices to be analysed. If C{None}, all vertices\n"
   "  will be considered.\n"
   "@return: cocitation scores for all given vertices in a matrix."},
  /* interface to igraph_similarity_dice */
  {"similarity_dice", (PyCFunction) igraphmodule_Graph_similarity_dice,
   METH_VARARGS | METH_KEYWORDS,
   "similarity_dice(vertices=None, pairs=None, mode=IGRAPH_ALL, loops=True)\n\n"
   "Dice similarity coefficient of vertices.\n\n"
   "The Dice similarity coefficient of two vertices is twice the number of\n"
   "their common neighbors divided by the sum of their degrees. This\n"
   "coefficient is very similar to the Jaccard coefficient, but usually\n"
   "gives higher similarities than its counterpart.\n\n"
   "@param vertices: the vertices to be analysed. If C{None} and I{pairs} is also\n"
   "  C{None}, all vertices will be considered.\n"
   "@param pairs: the vertex pairs to be analysed. If this is given, I{vertices}\n"
   "  must be C{None}, and the similarity values will be calculated only for the\n"
   "  given pairs. Vertex pairs must be specified as tuples of vertex IDs.\n"
   "@param mode: which neighbors should be considered for directed graphs.\n"
   "  Can be L{ALL}, L{IN} or L{OUT}, ignored for undirected graphs.\n"
   "@param loops: whether vertices should be considered adjacent to\n"
   "  themselves. Setting this to C{True} assumes a loop edge for all vertices\n"
   "  even if none is present in the graph. Setting this to C{False} may\n"
   "  result in strange results: nonadjacent vertices may have larger\n"
   "  similarities compared to the case when an edge is added between them --\n"
   "  however, this might be exactly the result you want to get.\n"
   "@return: the pairwise similarity coefficients for the vertices specified,\n"
   "  in the form of a matrix if C{pairs} is C{None} or in the form of a list\n"
   "  if C{pairs} is not C{None}.\n"
  },
  /* interface to igraph_similarity_inverse_log_weighted */
  {"similarity_inverse_log_weighted",
    (PyCFunction) igraphmodule_Graph_similarity_inverse_log_weighted,
   METH_VARARGS | METH_KEYWORDS,
   "similarity_inverse_log_weighted(vertices=None, mode=IGRAPH_ALL)\n\n"
   "Inverse log-weighted similarity coefficient of vertices.\n\n"
   "Each vertex is assigned a weight which is 1 / log(degree). The\n"
   "log-weighted similarity of two vertices is the sum of the weights\n"
   "of their common neighbors.\n\n"
   "@param vertices: the vertices to be analysed. If C{None}, all vertices\n"
   "  will be considered.\n"
   "@param mode: which neighbors should be considered for directed graphs.\n"
   "  Can be L{ALL}, L{IN} or L{OUT}, ignored for undirected graphs.\n"
   "  L{IN} means that the weights are determined by the out-degrees, L{OUT}\n"
   "  means that the weights are determined by the in-degrees.\n"
   "@return: the pairwise similarity coefficients for the vertices specified,\n"
   "  in the form of a matrix (list of lists).\n"
  },
  /* interface to igraph_similarity_jaccard */
  {"similarity_jaccard", (PyCFunction) igraphmodule_Graph_similarity_jaccard,
   METH_VARARGS | METH_KEYWORDS,
   "similarity_jaccard(vertices=None, pairs=None, mode=IGRAPH_ALL, loops=True)\n\n"
   "Jaccard similarity coefficient of vertices.\n\n"
   "The Jaccard similarity coefficient of two vertices is the number of their\n"
   "common neighbors divided by the number of vertices that are adjacent to\n"
   "at least one of them.\n\n"
   "@param vertices: the vertices to be analysed. If C{None} and I{pairs} is also\n"
   "  C{None}, all vertices will be considered.\n"
   "@param pairs: the vertex pairs to be analysed. If this is given, I{vertices}\n"
   "  must be C{None}, and the similarity values will be calculated only for the\n"
   "  given pairs. Vertex pairs must be specified as tuples of vertex IDs.\n"
   "@param mode: which neighbors should be considered for directed graphs.\n"
   "  Can be L{ALL}, L{IN} or L{OUT}, ignored for undirected graphs.\n"
   "@param loops: whether vertices should be considered adjacent to\n"
   "  themselves. Setting this to C{True} assumes a loop edge for all vertices\n"
   "  even if none is present in the graph. Setting this to C{False} may\n"
   "  result in strange results: nonadjacent vertices may have larger\n"
   "  similarities compared to the case when an edge is added between them --\n"
   "  however, this might be exactly the result you want to get.\n"
   "@return: the pairwise similarity coefficients for the vertices specified,\n"
   "  in the form of a matrix if C{pairs} is C{None} or in the form of a list\n"
   "  if C{pairs} is not C{None}.\n"
  },

  /******************/
  /* MOTIF COUNTING */
  /******************/
  {"motifs_randesu", (PyCFunction) igraphmodule_Graph_motifs_randesu,
   METH_VARARGS | METH_KEYWORDS,
   "motifs_randesu(size=3, cut_prob=None, callback=None)\n\n"
   "Counts the number of motifs in the graph\n\n"
   "Motifs are small subgraphs of a given structure in a graph. It is\n"
   "argued that the motif profile (ie. the number of different motifs in\n"
   "the graph) is characteristic for different types of networks and\n"
   "network function is related to the motifs in the graph.\n\n"
   "This function is able to find the different motifs of size three\n"
   "and four (ie. the number of different subgraphs with three and four\n"
   "vertices) in the network.\n\n"
   "In a big network the total number of motifs can be very large, so\n"
   "it takes a lot of time to find all of them. In such cases, a sampling\n"
   "method can be used. This function is capable of doing sampling via\n"
   "the I{cut_prob} argument. This argument gives the probability that\n"
   "a branch of the motif search tree will not be explored.\n\n"
   "@newfield ref: Reference\n"
   "@ref: S. Wernicke and F. Rasche: FANMOD: a tool for fast network\n"
   "  motif detection, Bioinformatics 22(9), 1152--1153, 2006.\n\n"
   "@param size: the size of the motifs (3 or 4).\n"
   "@param cut_prob: the cut probabilities for different levels of the search\n"
   "  tree. This must be a list of length I{size} or C{None} to find all\n"
   "  motifs.\n"
   "@param callback: C{None} or a callable that will be called for every motif\n"
   "  found in the graph. The callable must accept three parameters: the graph\n"
   "  itself, the list of vertices in the motif and the isomorphy class of the\n"
   "  motif (see L{Graph.isoclass()}). The search will stop when the callback\n"
   "  returns an object with a non-zero truth value or raises an exception.\n"
   "@return: the list of motifs if I{callback} is C{None}, or C{None} otherwise\n"
   "@see: Graph.motifs_randesu_no()\n"
  },
  {"motifs_randesu_no", (PyCFunction) igraphmodule_Graph_motifs_randesu_no,
   METH_VARARGS | METH_KEYWORDS,
   "motifs_randesu_no(size=3, cut_prob=None)\n\n"
   "Counts the total number of motifs in the graph\n\n"
   "Motifs are small subgraphs of a given structure in a graph.\n"
   "This function counts the total number of motifs in a graph without\n"
   "assigning isomorphism classes to them.\n\n"
   "@newfield ref: Reference\n"
   "@ref: S. Wernicke and F. Rasche: FANMOD: a tool for fast network\n"
   "  motif detection, Bioinformatics 22(9), 1152--1153, 2006.\n\n"
   "@param size: the size of the motifs (3 or 4).\n"
   "@param cut_prob: the cut probabilities for different levels of the search\n"
   "  tree. This must be a list of length I{size} or C{None} to find all\n"
   "  motifs.\n"
   "@see: Graph.motifs_randesu()\n"
  },
  {"motifs_randesu_estimate",
   (PyCFunction) igraphmodule_Graph_motifs_randesu_estimate,
   METH_VARARGS | METH_KEYWORDS,
   "motifs_randesu_estimate(size=3, cut_prob=None, sample)\n\n"
   "Counts the total number of motifs in the graph\n\n"
   "Motifs are small subgraphs of a given structure in a graph.\n"
   "This function estimates the total number of motifs in a graph without\n"
   "assigning isomorphism classes to them by extrapolating from a random\n"
   "sample of vertices.\n\n"
   "@newfield ref: Reference\n"
   "@ref: S. Wernicke and F. Rasche: FANMOD: a tool for fast network\n"
   "  motif detection, Bioinformatics 22(9), 1152--1153, 2006.\n\n"
   "@param size: the size of the motifs (3 or 4).\n"
   "@param cut_prob: the cut probabilities for different levels of the search\n"
   "  tree. This must be a list of length I{size} or C{None} to find all\n"
   "  motifs.\n"
   "@param sample: the size of the sample or the vertex IDs of the vertices\n"
   "  to be used for sampling.\n"
   "@see: Graph.motifs_randesu()\n"
  },
  {"dyad_census", (PyCFunction) igraphmodule_Graph_dyad_census,
   METH_NOARGS,
   "dyad_census()\n\n"
   "Dyad census, as defined by Holland and Leinhardt\n\n"
   "Dyad census means classifying each pair of vertices of a directed\n"
   "graph into three categories: mutual, there is an edge from I{a} to\n"
   "I{b} and also from I{b} to I{a}; asymmetric, there is an edge\n"
   "either from I{a} to I{b} or from I{b} to I{a} but not the other way\n"
   "and null, no edges between I{a} and I{b}.\n\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a L{DyadCensus} object.\n"
   "  It is advised to use that.\n\n"
   "@return: the number of mutual, asymmetric and null connections in a\n"
   "  3-tuple."
  },
  {"triad_census", (PyCFunction) igraphmodule_Graph_triad_census,
   METH_NOARGS,
   "triad_census()\n\n"
   "Triad census, as defined by Davis and Leinhardt\n\n"
   "Calculating the triad census means classifying every triplets of\n"
   "vertices in a directed graph. A triplet can be in one of 16 states,\n"
   "these are listed in the documentation of the C interface of igraph.\n"
   "\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a L{TriadCensus} object.\n"
   "  It is advised to use that. The name of the triplet classes are\n"
   "  also documented there.\n\n"
  },

  /********************/
  /* LAYOUT FUNCTIONS */
  /********************/

  /* interface to igraph_layout_circle */
  {"layout_circle", (PyCFunction) igraphmodule_Graph_layout_circle,
   METH_VARARGS | METH_KEYWORDS,
   "layout_circle(dim=2)\n\n"
   "Places the vertices of the graph uniformly on a circle or a sphere.\n\n"
   "@param dim: the desired number of dimensions for the layout. dim=2\n"
   "  means a 2D layout, dim=3 means a 3D layout.\n"
   "@return: the calculated layout."},

  /* interface to igraph_layout_grid */
  {"layout_grid", (PyCFunction) igraphmodule_Graph_layout_grid,
   METH_VARARGS | METH_KEYWORDS,
   "layout_grid(width=0, height=0, dim=2)\n\n"
   "Places the vertices of a graph in a 2D or 3D grid.\n\n"
   "@param width: the number of vertices in a single row of the layout.\n"
   "  Zero or negative numbers mean that the width should be determined\n"
   "  automatically.\n"
   "@param height: the number of vertices in a single column of the layout.\n"
   "  Zero or negative numbers mean that the height should be determined\n"
   "  automatically. It must not be given if the number of dimensions is 2.\n"
   "@param dim: the desired number of dimensions for the layout. dim=2\n"
   "  means a 2D layout, dim=3 means a 3D layout.\n"
   "@return: the calculated layout."},

  /* interface to igraph_layout_star */
  {"layout_star", (PyCFunction) igraphmodule_Graph_layout_star,
   METH_VARARGS | METH_KEYWORDS,
   "layout_star(center=0, order=None)\n\n"
   "Calculates a star-like layout for the graph.\n\n"
   "@param center: the ID of the vertex to put in the center\n"
   "@param order: a numeric vector giving the order of the vertices\n"
   "  (including the center vertex!). If it is C{None}, the vertices\n"
   "  will be placed in increasing vertex ID order.\n"
   "@return: the calculated layout."
  },

  /* interface to igraph_layout_kamada_kawai */
  {"layout_kamada_kawai",
   (PyCFunction) igraphmodule_Graph_layout_kamada_kawai,
   METH_VARARGS | METH_KEYWORDS,
   "layout_kamada_kawai(maxiter=1000, sigma=None, initemp=10, coolexp=0.99,\n"
   "  kkconst=None, seed=None, minx=None, maxx=None, miny=None, maxy=None, \n"
   "  minz=None, maxz=None, dim=2)\n\n"
   "Places the vertices on a plane according to the Kamada-Kawai algorithm.\n\n"
   "This is a force directed layout, see Kamada, T. and Kawai, S.:\n"
   "An Algorithm for Drawing General Undirected Graphs.\n"
   "Information Processing Letters, 31/1, 7--15, 1989.\n\n"
   "@param maxiter: the number of iterations to perform.\n"
   "@param sigma: the standard base deviation of the position\n"
   "  change proposals. C{None} means the number of vertices / 4\n"
   "@param initemp: initial temperature of the simulated annealing.\n"
   "@param coolexp: cooling exponent of the simulated annealing.\n"
   "@param kkconst: the Kamada-Kawai vertex attraction constant.\n"
   "  C{None} means the square of the number of vertices.\n"
   "@param minx: if not C{None}, it must be a vector with exactly as many\n"
   "  elements as there are vertices in the graph. Each element is a\n"
   "  minimum constraint on the X value of the vertex in the layout.\n"
   "@param maxx: similar to I{minx}, but with maximum constraints\n"
   "@param miny: similar to I{minx}, but with the Y coordinates\n"
   "@param maxy: similar to I{maxx}, but with the Y coordinates\n"
   "@param minz: similar to I{minx}, but with the Z coordinates. Use only\n"
   "  for 3D layouts (C{dim}=3).\n"
   "@param maxz: similar to I{maxx}, but with the Z coordinates. Use only\n"
   "  for 3D layouts (C{dim}=3).\n"
   "@param seed: if C{None}, uses a random starting layout for the\n"
   "  algorithm. If a matrix (list of lists), uses the given matrix\n"
   "  as the starting position.\n"
   "@param dim: the desired number of dimensions for the layout. dim=2\n"
   "  means a 2D layout, dim=3 means a 3D layout.\n"
   "@return: the calculated layout."
  },

  /* interface to igraph_layout_drl */
  {"layout_drl",
   (PyCFunction) igraphmodule_Graph_layout_drl,
   METH_VARARGS | METH_KEYWORDS,
   "layout_drl(weights=None, fixed=None, seed=None, options=None, dim=2)\n\n"
   "Places the vertices on a 2D plane or in the 3D space ccording to the DrL\n"
   "layout algorithm.\n\n"
   "This is an algorithm suitable for quite large graphs, but it can be\n"
   "surprisingly slow for small ones (where the simpler force-based layouts\n"
   "like C{layout_kamada_kawai()} or C{layout_fruchterman_reingold()} are\n"
   "more useful.\n\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@param seed: if C{None}, uses a random starting layout for the\n"
   "  algorithm. If a matrix (list of lists), uses the given matrix\n"
   "  as the starting position.\n"
   "@param fixed: if a seed is given, you can specify some vertices to be\n"
   "  kept fixed at their original position in the seed by passing an\n"
   "  appropriate list here. The list must have exactly as many items as\n"
   "  the number of vertices in the graph. Items of the list that evaluate\n"
   "  to C{True} denote vertices that will not be moved.\n"
   "@param options: if you give a string argument here, you can select from\n"
   "  five default preset parameterisations: C{default}, C{coarsen} for a\n"
   "  coarser layout, C{coarsest} for an even coarser layout, C{refine} for\n"
   "  refining an existing layout and C{final} for finalizing a layout. If\n"
   "  you supply an object that is not a string, the DrL layout parameters\n"
   "  are retrieved from the respective keys of the object (so it should\n"
   "  be a dict or something else that supports the mapping protocol).\n"
   "  The following keys can be used:\n"
   "  \n"
   "    - C{edge_cut}: edge cutting is done in the late stages of the\n"
   "      algorithm in order to achieve less dense layouts. Edges are\n"
   "      cut if there is a lot of stress on them (a large value in the\n"
   "      objective function sum). The edge cutting parameter is a value\n"
   "      between 0 and 1 with 0 representing no edge cutting and 1\n"
   "      representing maximal edge cutting.\n\n"
   "    - C{init_iterations}: number of iterations in the initialization\n"
   "      phase\n\n"
   "    - C{init_temperature}: start temperature during initialization\n\n"
   "    - C{init_attraction}: attraction during initialization\n\n"
   "    - C{init_damping_mult}: damping multiplier during initialization\n\n"
   "    - C{liquid_iterations}, C{liquid_temperature}, C{liquid_attraction},\n"
   "      C{liquid_damping_mult}: same parameters for the liquid phase\n\n"
   "    - C{expansion_iterations}, C{expansion_temperature},\n"
   "      C{expansion_attraction}, C{expansion_damping_mult}:\n"
   "      parameters for the expansion phase\n\n"
   "    - C{cooldown_...}: parameters for the cooldown phase\n\n"
   "    - C{crunch_...}: parameters for the crunch phase\n\n"
   "    - C{simmer_...}: parameters for the simmer phase\n\n"
   "  \n"
   "  Instead of a mapping, you can also use an arbitrary Python object\n"
   "  here: if the object does not support the mapping protocol, an\n"
   "  attribute of the object with the same name is looked up instead. If\n"
   "  a parameter cannot be found either as a key or an attribute, the\n"
   "  default from the C{default} preset will be used.\n\n"
   "@param dim: the desired number of dimensions for the layout. dim=2\n"
   "  means a 2D layout, dim=3 means a 3D layout.\n"
   "@return: the calculated layout."
  },

  /* interface to igraph_layout_fruchterman_reingold */
  {"layout_fruchterman_reingold",
   (PyCFunction) igraphmodule_Graph_layout_fruchterman_reingold,
   METH_VARARGS | METH_KEYWORDS,
   "layout_fruchterman_reingold(weights=None, maxiter=500, maxdelta=None, area=None,\n"
   "  coolexp=1.5, repulserad=None, seed=None, minx=None, maxx=None, miny=None, \n"
   "  maxy=None, minz=None, maxz=None, dim=2)\n\n"
   "Places the vertices on a 2D plane or in the 3D space according to the\n"
   "Fruchterman-Reingold algorithm.\n\n"
   "This is a force directed layout, see Fruchterman, T. M. J. and Reingold, E. M.:\n"
   "Graph Drawing by Force-directed Placement.\n"
   "Software -- Practice and Experience, 21/11, 1129--1164, 1991\n\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@param maxiter: the number of iterations to perform. The default\n"
   "  is 500.\n"
   "@param maxdelta: the maximum distance to move a vertex in\n"
   "  an iteration. The default is the number of vertices.\n"
   "@param area: the area of the square on which the vertices\n"
   "  will be placed. The default is the square of the number of\n"
   "  vertices.\n"
   "@param coolexp: the cooling exponent of the simulated annealing.\n"
   "  The default is 1.5.\n"
   "@param repulserad: determines the radius at which vertex-vertex\n"
   "  repulsion cancels out attraction of adjacent vertices.\n"
   "  The default is the number of vertices^3.\n"
   "@param minx: if not C{None}, it must be a vector with exactly as many\n"
   "  elements as there are vertices in the graph. Each element is a\n"
   "  minimum constraint on the X value of the vertex in the layout.\n"
   "@param maxx: similar to I{minx}, but with maximum constraints\n"
   "@param miny: similar to I{minx}, but with the Y coordinates\n"
   "@param maxy: similar to I{maxx}, but with the Y coordinates\n"
   "@param minz: similar to I{minx}, but with the Z coordinates. Use only\n"
   "  for 3D layouts (C{dim}=3).\n"
   "@param maxz: similar to I{maxx}, but with the Z coordinates. Use only\n"
   "  for 3D layouts (C{dim}=3).\n"
   "@param seed: if C{None}, uses a random starting layout for the\n"
   "  algorithm. If a matrix (list of lists), uses the given matrix\n"
   "  as the starting position.\n"
   "@param dim: the desired number of dimensions for the layout. dim=2\n"
   "  means a 2D layout, dim=3 means a 3D layout.\n"
   "@return: the calculated layout."
  },

  /* interface to igraph_layout_graphopt */
  {"layout_graphopt",
   (PyCFunction) igraphmodule_Graph_layout_graphopt,
   METH_VARARGS | METH_KEYWORDS,
   "layout_graphopt(niter=500, node_charge=0.001, node_mass=30, spring_length=0, spring_constant=1, max_sa_movement=5, seed=None)\n\n"
   "This is a port of the graphopt layout algorithm by Michael Schmuhl.\n"
   "graphopt version 0.4.1 was rewritten in C and the support for layers\n"
   "was removed.\n\n"
   "graphopt uses physical analogies for defining attracting and repelling\n"
   "forces among the vertices and then the physical system is simulated\n"
   "until it reaches an equilibrium or the maximal number of iterations is\n"
   "reached.\n\n"
   "See U{http://www.schmuhl.org/graphopt/} for the original graphopt.\n\n"
   "@param niter: the number of iterations to perform. Should be a couple\n"
   "  of hundred in general.\n\n"
   "@param node_charge: the charge of the vertices, used to calculate electric\n"
   "  repulsion.\n"
   "@param node_mass: the mass of the vertices, used for the spring forces\n"
   "@param spring_length: the length of the springs\n"
   "@param spring_constant: the spring constant\n"
   "@param max_sa_movement: the maximum amount of movement allowed in a single\n"
   "  step along a single axis.\n"
   "@param seed: a matrix containing a seed layout from which the algorithm\n"
   "  will be started. If C{None}, a random layout will be used.\n"
   "@return: the calculated layout."
  },

  /* interface to igraph_layout_grid_fruchterman_reingold */
  {"layout_grid_fruchterman_reingold",
   (PyCFunction) igraphmodule_Graph_layout_grid_fruchterman_reingold,
   METH_VARARGS | METH_KEYWORDS,
   "layout_grid_fruchterman_reingold(maxiter=500, maxdelta=None, area=None, coolexp=0.99, repulserad=maxiter*maxdelta, cellsize=None, seed=None)\n\n"
   "Places the vertices on a 2D plane according to the Fruchterman-Reingold\n"
   "grid algorithm.\n\n"
   "This is a modified version of a force directed layout, see\n"
   "Fruchterman, T. M. J. and Reingold, E. M.:\n"
   "Graph Drawing by Force-directed Placement.\n"
   "Software -- Practice and Experience, 21/11, 1129--1164, 1991.\n"
   "The algorithm partitions the 2D space to a grid and vertex\n"
   "repulsion is then calculated only for vertices nearby.\n\n"
   "@param maxiter: the number of iterations to perform.\n"
   "@param maxdelta: the maximum distance to move a vertex in\n"
   "  an iteration. C{None} means the number of vertices.\n"
   "@param area: the area of the square on which the vertices\n"
   "  will be placed. C{None} means the square of M{maxdelta}.\n"
   "@param coolexp: the cooling exponent of the simulated annealing.\n"
   "@param repulserad: determines the radius at which vertex-vertex\n"
   "  repulsion cancels out attraction of adjacent vertices.\n"
   "  C{None} means M{maxiter*maxdelta}.\n"
   "@param cellsize: the size of the grid cells. When calculating the\n"
   "  repulsion forces, only vertices in the same or neighboring\n"
   "  grid cells are taken into account. Defaults to the fourth\n"
   "  root of M{area}.\n"
   "@param seed: if C{None}, uses a random starting layout for the\n"
   "  algorithm. If a matrix (list of lists), uses the given matrix\n"
   "  as the starting position.\n"
   "@return: the calculated layout."
  },

  /* interface to igraph_layout_lgl */
  {"layout_lgl", (PyCFunction) igraphmodule_Graph_layout_lgl,
   METH_VARARGS | METH_KEYWORDS,
   "layout_lgl(maxiter=150, maxdelta=-1, area=-1, coolexp=1.5, repulserad=-1, cellsize=-1, root=None)\n\n"
   "Places the vertices on a 2D plane according to the Large Graph Layout.\n\n"
   "@param maxiter: the number of iterations to perform.\n"
   "@param maxdelta: the maximum distance to move a vertex in\n"
   "  an iteration. If negative, defaults to the number of vertices.\n"
   "@param area: the area of the square on which the vertices\n"
   "  will be placed. If negative, defaults to the number of vertices\n"
   "  squared.\n"
   "@param coolexp: the cooling exponent of the simulated annealing.\n"
   "@param repulserad: determines the radius at which vertex-vertex\n"
   "  repulsion cancels out attraction of adjacent vertices.\n"
   "  If negative, defaults to M{area} times the number of vertices.\n"
   "@param cellsize: the size of the grid cells. When calculating the\n"
   "  repulsion forces, only vertices in the same or neighboring\n"
   "  grid cells are taken into account. Defaults to the fourth\n"
   "  root of M{area}.\n"
   "@param root: the root vertex, this is placed first, its neighbors\n"
   "  in the first iteration, second neighbors in the second,\n"
   "  etc. C{None} means that a random vertex will be chosen.\n"
   "@return: the calculated layout."
  },

  /* interface to igraph_layout_mds */
  {"layout_mds",
   (PyCFunction) igraphmodule_Graph_layout_mds,
   METH_VARARGS | METH_KEYWORDS,
   "layout_mds(dist=None, dim=2, arpack_options=None)\n"
   "Places the vertices in an Euclidean space with the given number of\n"
   "dimensions using multidimensional scaling.\n\n"
   "This layout requires a distance matrix, where the intersection of\n"
   "row M{i} and column M{j} specifies the desired distance between\n"
   "vertex M{i} and vertex M{j}. The algorithm will try to place the\n"
   "vertices in a way that approximates the distance relations\n"
   "prescribed in the distance matrix. igraph uses the classical\n"
   "multidimensional scaling by Torgerson (see reference below).\n\n"
   "For unconnected graphs, the method will decompose the graph into\n"
   "weakly connected components and then lay out the components\n"
   "individually using the appropriate parts of the distance matrix.\n\n"
   "@param dist: the distance matrix. It must be symmetric and the\n"
   "  symmetry is not checked -- results are unspecified when a\n"
   "  non-symmetric distance matrix is used. If this parameter is\n"
   "  C{None}, the shortest path lengths will be used as distances.\n"
   "  Directed graphs are treated as undirected when calculating\n"
   "  the shortest path lengths to ensure symmetry.\n"
   "@param dim: the number of dimensions. For 2D layouts, supply\n"
   "  2 here; for 3D layouts, supply 3.\n"
   "@param arpack_options: an L{ARPACKOptions} object used to fine-tune\n"
   "  the ARPACK eigenvector calculation. If omitted, the module-level\n"
   "  variable called C{arpack_options} is used.\n"
   "@return: the calculated layout.\n\n"
   "@newfield ref: Reference\n"
   "@ref: Cox & Cox: Multidimensional Scaling (1994), Chapman and\n"
   "  Hall, London.\n"
  },

  /* interface to igraph_layout_reingold_tilford */
  {"layout_reingold_tilford",
   (PyCFunction) igraphmodule_Graph_layout_reingold_tilford,
   METH_VARARGS | METH_KEYWORDS,
   "layout_reingold_tilford(mode=\"out\", root=None, rootlevel=None)\n"
   "Places the vertices on a 2D plane according to the Reingold-Tilford\n"
   "layout algorithm.\n\n"
   "This is a tree layout. If the given graph is not a tree, a breadth-first\n"
   "search is executed first to obtain a possible spanning tree.\n\n"
   "@param mode: specifies which edges to consider when builing the tree.\n"
   "  If it is C{OUT} then only the outgoing, if it is C{IN} then only the\n"
   "  incoming edges of a parent are considered. If it is C{ALL} then all\n"
   "  edges are used (this was the behaviour in igraph 0.5 and before).\n"
   "  This parameter also influences how the root vertices are calculated\n"
   "  if they are not given. See the I{root} parameter.\n"
   "@param root: the index of the root vertex or root vertices.\n"
   "  if this is a non-empty vector then the supplied vertex IDs are\n"
   "  used as the roots of the trees (or a single tree if the graph is\n"
   "  connected. If this is C{None} or an empty list, the root vertices\n"
   "  are automatically calculated based on topological sorting,\n"
   "  performed with the opposite of the I{mode} argument.\n"
   "@param rootlevel: this argument is useful when drawing forests which are\n"
   "  not trees. It specifies the level of the root vertices for every tree\n"
   "  in the forest.\n"
   "@return: the calculated layout.\n\n"
   "@see: layout_reingold_tilford_circular\n"
   "@newfield ref: Reference\n"
   "@ref: EM Reingold, JS Tilford: I{Tidier Drawings of Trees.}\n"
   "IEEE Transactions on Software Engineering 7:22, 223-228, 1981."},

  /* interface to igraph_layout_reingold_tilford_circular */
  {"layout_reingold_tilford_circular",
   (PyCFunction) igraphmodule_Graph_layout_reingold_tilford_circular,
   METH_VARARGS | METH_KEYWORDS,
   "layout_reingold_tilford_circular(mode=\"out\", root=None, rootlevel=None)\n"
   "Circular Reingold-Tilford layout for trees.\n\n"
   "This layout is similar to the Reingold-Tilford layout, but the vertices\n"
   "are placed in a circular way, with the root vertex in the center.\n\n"
   "See L{layout_reingold_tilford} for the explanation of the parameters.\n\n"
   "@return: the calculated layout.\n\n"
   "@see: layout_reingold_tilford\n"
   "@newfield ref: Reference\n"
   "@ref: EM Reingold, JS Tilford: I{Tidier Drawings of Trees.}\n"
   "IEEE Transactions on Software Engineering 7:22, 223-228, 1981."},

  /* interface to igraph_layout_random */
  {"layout_random", (PyCFunction) igraphmodule_Graph_layout_random,
   METH_VARARGS | METH_KEYWORDS,
   "layout_random(dim=2)\n"
   "Places the vertices of the graph randomly.\n\n"
   "@param dim: the desired number of dimensions for the layout. dim=2\n"
   "  means a 2D layout, dim=3 means a 3D layout.\n"
   "@return: the coordinate pairs in a list."},

  /* interface to igraph_layout_sugiyama */
  {"_layout_sugiyama",
   (PyCFunction) igraphmodule_Graph_layout_sugiyama,
   METH_VARARGS | METH_KEYWORDS,
   "Internal function, undocumented.\n\n"
   "@see: Graph.layout_sugiyama()\n\n"},

  ////////////////////////////
  // VISITOR-LIKE FUNCTIONS //
  ////////////////////////////
  {"bfs", (PyCFunction) igraphmodule_Graph_bfs,
   METH_VARARGS | METH_KEYWORDS,
   "bfs(vid, mode=OUT)\n\n"
   "Conducts a breadth first search (BFS) on the graph.\n\n"
   "@param vid: the root vertex ID\n"
   "@param mode: either L{IN} or L{OUT} or L{ALL}, ignored\n"
   "  for undirected graphs.\n"
   "@return: a tuple with the following items:\n"
   "   - The vertex IDs visited (in order)\n"
   "   - The start indices of the layers in the vertex list\n"
   "   - The parent of every vertex in the BFS\n"},
  {"bfsiter", (PyCFunction) igraphmodule_Graph_bfsiter,
   METH_VARARGS | METH_KEYWORDS,
   "bfsiter(vid, mode=OUT, advanced=False)\n\n"
   "Constructs a breadth first search (BFS) iterator of the graph.\n\n"
   "@param vid: the root vertex ID\n"
   "@param mode: either L{IN} or L{OUT} or L{ALL}.\n"
   "@param advanced: if C{False}, the iterator returns the next\n"
   "  vertex in BFS order in every step. If C{True}, the iterator\n"
   "  returns the distance of the vertex from the root and the\n"
   "  parent of the vertex in the BFS tree as well.\n"
   "@return: the BFS iterator as an L{igraph.BFSIter} object.\n"},

  /////////////////
  // CONVERSIONS //
  /////////////////

  // interface to igraph_get_adjacency
  {"get_adjacency", (PyCFunction) igraphmodule_Graph_get_adjacency,
   METH_VARARGS | METH_KEYWORDS,
   "get_adjacency(type=GET_ADJACENCY_BOTH, eids=False)\n\n"
   "Returns the adjacency matrix of a graph.\n\n"
   "@param type: either C{GET_ADJACENCY_LOWER} (uses the\n"
   "  lower triangle of the matrix) or C{GET_ADJACENCY_UPPER}\n"
   "  (uses the upper triangle) or C{GET_ADJACENCY_BOTH}\n"
   "  (uses both parts). Ignored for directed graphs.\n"
   "@param eids: if C{True}, the result matrix will contain\n"
   "  zeros for non-edges and the ID of the edge plus one\n"
   "  for edges in the appropriate cell. If C{False}, the\n"
   "  result matrix will contain the number of edges for\n"
   "  each vertex pair.\n"
   "@return: the adjacency matrix.\n"},

  // interface to igraph_get_edgelist
  {"get_edgelist", (PyCFunction) igraphmodule_Graph_get_edgelist,
   METH_NOARGS,
   "get_edgelist()\n\n" "Returns the edge list of a graph."},

  /* interface to igraph_get_incidence */
  {"get_incidence", (PyCFunction) igraphmodule_Graph_get_incidence,
   METH_VARARGS | METH_KEYWORDS,
   "get_incidence(types)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: Graph.get_incidence()\n\n"},

  // interface to igraph_to_directed
  {"to_directed", (PyCFunction) igraphmodule_Graph_to_directed,
   METH_VARARGS | METH_KEYWORDS,
   "to_directed(mutual=True)\n\n"
   "Converts an undirected graph to directed.\n\n"
   "@param mutual: C{True} if mutual directed edges should be\n"
   "  created for every undirected edge. If C{False}, a directed\n"
   "  edge with arbitrary direction is created.\n"},

  // interface to igraph_to_undirected
  {"to_undirected", (PyCFunction) igraphmodule_Graph_to_undirected,
   METH_VARARGS | METH_KEYWORDS,
   "to_undirected(mode=\"collapse\", combine_edges=None)\n\n"
   "Converts a directed graph to undirected.\n\n"
   "@param mode: specifies what to do with multiple directed edges\n"
   "  going between the same vertex pair. C{True} or C{\"collapse\"}\n"
   "  means that only a single edge should be created from multiple\n"
   "  directed edges. C{False} or C{\"each\"} means that every edge\n"
   "  will be kept (with the arrowheads removed). C{\"mutual\"}\n"
   "  creates one undirected edge for each mutual directed edge pair.\n"
   "@param combine_edges: specifies how to combine the attributes of\n"
   "  multiple edges between the same pair of vertices into a single\n"
   "  attribute. See L{Graph.simplify()} for more details.\n"
  },

  /* interface to igraph_laplacian */
  {"laplacian", (PyCFunction) igraphmodule_Graph_laplacian,
   METH_VARARGS | METH_KEYWORDS,
   "laplacian(weights=None, normalized=False)\n\n"
   "Returns the Laplacian matrix of a graph.\n\n"
   "The Laplacian matrix is similar to the adjacency matrix, but the edges\n"
   "are denoted with -1 and the diagonal contains the node degrees.\n\n"
   "Normalized Laplacian matrices have 1 or 0 in their diagonals (0 for vertices\n"
   "with no edges), edges are denoted by 1 / sqrt(d_i * d_j) where d_i is the\n"
   "degree of node i.\n\n"
   "Multiple edges and self-loops are silently ignored. Although it is\n"
   "possible to calculate the Laplacian matrix of a directed graph, it does\n"
   "not make much sense.\n\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name. When edge weights are used, the degree\n"
   "  of a node is considered to be the weight of its incident edges.\n"
   "@param normalized: whether to return the normalized Laplacian matrix.\n"
   "@return: the Laplacian matrix.\n"},

  ///////////////////////////////
  // LOADING AND SAVING GRAPHS //
  ///////////////////////////////

  // interface to igraph_read_graph_dimacs
  {"Read_DIMACS", (PyCFunction) igraphmodule_Graph_Read_DIMACS,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_DIMACS(f, directed=False)\n\n"
   "Reads a graph from a file conforming to the DIMACS minimum-cost flow file format.\n\n"
   "For the exact description of the format, see\n"
   "U{http://lpsolve.sourceforge.net/5.5/DIMACS.htm}\n\n"
   "Restrictions compared to the official description of the format:\n\n"
   "  - igraph's DIMACS reader requires only three fields in an arc definition,\n"
   "    describing the edge's source and target node and its capacity.\n"
   "  - Source vertices are identified by 's' in the FLOW field, target vertices are\n"
   "    identified by 't'.\n"
   "  - Node indices start from 1. Only a single source and target node is allowed.\n\n"
   "@param f: the name of the file or a Python file handle\n"
   "@param directed: whether the generated graph should be directed.\n"
   "@return: the generated graph, the source and the target of the flow and the edge\n"
   "  capacities in a tuple\n"},

  /* interface to igraph_read_graph_dl */
  {"Read_DL", (PyCFunction) igraphmodule_Graph_Read_DL,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_DL(f, directed=True)\n\n"
   "Reads an UCINET DL file and creates a graph based on it.\n\n"
   "@param f: the name of the file or a Python file handle\n"
   "@param directed: whether the generated graph should be directed.\n"},

  /* interface to igraph_read_graph_edgelist */
  {"Read_Edgelist", (PyCFunction) igraphmodule_Graph_Read_Edgelist,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_Edgelist(f, directed=True)\n\n"
   "Reads an edge list from a file and creates a graph based on it.\n\n"
   "Please note that the vertex indices are zero-based.\n\n"
   "@param f: the name of the file or a Python file handle\n"
   "@param directed: whether the generated graph should be directed.\n"},
  /* interface to igraph_read_graph_graphdb */
  {"Read_GraphDB", (PyCFunction) igraphmodule_Graph_Read_GraphDB,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_GraphDB(f, directed=False)\n\n"
   "Reads a GraphDB format file and creates a graph based on it.\n\n"
   "GraphDB is a binary format, used in the graph database for\n"
   "isomorphism testing (see U{http://amalfi.dis.unina.it/graph/}).\n\n"
   "@param f: the name of the file or a Python file handle\n"
   "@param directed: whether the generated graph should be directed.\n"},
  /* interface to igraph_read_graph_graphml */
  {"Read_GraphML", (PyCFunction) igraphmodule_Graph_Read_GraphML,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_GraphML(f, directed=True, index=0)\n\n"
   "Reads a GraphML format file and creates a graph based on it.\n\n"
   "@param f: the name of the file or a Python file handle\n"
   "@param index: if the GraphML file contains multiple graphs,\n"
   "  specifies the one that should be loaded. Graph indices\n"
   "  start from zero, so if you want to load the first graph,\n"
   "  specify 0 here.\n"},
  /* interface to igraph_read_graph_gml */
  {"Read_GML", (PyCFunction) igraphmodule_Graph_Read_GML,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_GML(f)\n\n"
   "Reads a GML file and creates a graph based on it.\n\n"
   "@param f: the name of the file or a Python file handle\n"
  },
  /* interface to igraph_read_graph_ncol */
  {"Read_Ncol", (PyCFunction) igraphmodule_Graph_Read_Ncol,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_Ncol(f, names=True, weights=\"if_present\", directed=True)\n\n"
   "Reads an .ncol file used by LGL.\n\n"
   "It is also useful for creating graphs from \"named\" (and\n"
   "optionally weighted) edge lists.\n\n"
   "This format is used by the Large Graph Layout program. See the\n"
   "U{documentation of LGL <http://bioinformatics.icmb.utexas.edu/lgl/>}\n"
   "regarding the exact format description.\n\n"
   "LGL originally cannot deal with graphs containing multiple or loop\n"
   "edges, but this condition is not checked here, as igraph is happy\n"
   "with these.\n\n"
   "@param f: the name of the file or a Python file handle\n"
   "@param names: If C{True}, the vertex names are added as a\n"
   "  vertex attribute called 'name'.\n"
   "@param weights: If True, the edge weights are added as an\n"
   "  edge attribute called 'weight', even if there are no\n"
   "  weights in the file. If False, the edge weights are never\n"
   "  added, even if they are present. C{\"auto\"} or C{\"if_present\"}\n"
   "  means that weights are added if there is at least one weighted\n"
   "  edge in the input file, but they are not added otherwise.\n"
   "@param directed: whether the graph being created should be\n"
   "  directed\n"
  },
  /* interface to igraph_read_graph_lgl */
  {"Read_Lgl", (PyCFunction) igraphmodule_Graph_Read_Lgl,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_Lgl(f, names=True, weights=\"if_present\", directed=True)\n\n"
   "Reads an .lgl file used by LGL.\n\n"
   "It is also useful for creating graphs from \"named\" (and\n"
   "optionally weighted) edge lists.\n\n"
   "This format is used by the Large Graph Layout program. See the\n"
   "U{documentation of LGL <http://bioinformatics.icmb.utexas.edu/lgl/>}\n"
   "regarding the exact format description.\n\n"
   "LGL originally cannot deal with graphs containing multiple or loop\n"
   "edges, but this condition is not checked here, as igraph is happy\n"
   "with these.\n\n"
   "@param f: the name of the file or a Python file handle\n"
   "@param names: If C{True}, the vertex names are added as a\n"
   "  vertex attribute called 'name'.\n"
   "@param weights: If True, the edge weights are added as an\n"
   "  edge attribute called 'weight', even if there are no\n"
   "  weights in the file. If False, the edge weights are never\n"
   "  added, even if they are present. C{\"auto\"} or C{\"if_present\"}\n"
   "  means that weights are added if there is at least one weighted\n"
   "  edge in the input file, but they are not added otherwise.\n"
   "@param directed: whether the graph being created should be\n"
   "  directed\n"
  },
  /* interface to igraph_read_graph_pajek */
  {"Read_Pajek", (PyCFunction) igraphmodule_Graph_Read_Pajek,
   METH_VARARGS | METH_KEYWORDS | METH_CLASS,
   "Read_Pajek(f)\n\n"
   "Reads a Pajek format file and creates a graph based on it.\n\n"
   "@param f: the name of the file or a Python file handle\n"},
  /* interface to igraph_write_graph_dimacs */
  {"write_dimacs", (PyCFunction) igraphmodule_Graph_write_dimacs,
   METH_VARARGS | METH_KEYWORDS,
   "write_dimacs(f, source, target, capacity=None)\n\n"
   "Writes the graph in DIMACS format to the given file.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
   "@param source: the source vertex ID\n"
   "@param target: the target vertex ID\n"
   "@param capacity: the capacities of the edges in a list. If it is not a\n"
   "  list, the corresponding edge attribute will be used to retrieve\n"
   "  capacities."},
  /* interface to igraph_write_graph_dot */
  {"write_dot", (PyCFunction) igraphmodule_Graph_write_dot,
   METH_VARARGS | METH_KEYWORDS,
   "write_dot(f)\n\n"
   "Writes the graph in DOT format to the given file.\n\n"
   "DOT is the format used by the U{GraphViz <http://www.graphviz.org>}\n"
   "software package.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
   },
  /* interface to igraph_write_graph_edgelist */
  {"write_edgelist", (PyCFunction) igraphmodule_Graph_write_edgelist,
   METH_VARARGS | METH_KEYWORDS,
   "write_edgelist(f)\n\n"
   "Writes the edge list of a graph to a file.\n\n"
   "Directed edges are written in (from, to) order.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"},
  /* interface to igraph_write_graph_gml */
  {"write_gml", (PyCFunction) igraphmodule_Graph_write_gml,
   METH_VARARGS | METH_KEYWORDS,
   "write_gml(f, creator=None, ids=None)\n\n"
   "Writes the graph in GML format to the given file.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
   "@param creator: optional creator information to be written to the file.\n"
   "  If C{None}, the current date and time is added.\n"
   "@param ids: optional numeric vertex IDs to use in the file. This must\n"
   "  be a list of integers or C{None}. If C{None}, the C{id} attribute of\n"
   "  the vertices are used, or if they don't exist, numeric vertex IDs\n"
   "  will be generated automatically."},
  /* interface to igraph_write_graph_ncol */
  {"write_ncol", (PyCFunction) igraphmodule_Graph_write_ncol,
   METH_VARARGS | METH_KEYWORDS,
   "write_ncol(f, names=\"name\", weights=\"weights\")\n\n"
   "Writes the edge list of a graph to a file in .ncol format.\n\n"
   "Note that multiple edges and/or loops break the LGL software,\n"
   "but igraph does not check for this condition. Unless you know\n"
   "that the graph does not have multiple edges and/or loops, it\n"
   "is wise to call L{simplify()} before saving.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
   "@param names: the name of the vertex attribute containing the name\n"
   "  of the vertices. If you don't want to store vertex names,\n"
   "  supply C{None} here.\n"
   "@param weights: the name of the edge attribute containing the weight\n"
   "  of the vertices. If you don't want to store weights,\n"
   "  supply C{None} here.\n"},
  /* interface to igraph_write_graph_lgl */
  {"write_lgl", (PyCFunction) igraphmodule_Graph_write_lgl,
   METH_VARARGS | METH_KEYWORDS,
   "write_lgl(f, names=\"name\", weights=\"weights\", isolates=True)\n\n"
   "Writes the edge list of a graph to a file in .lgl format.\n\n"
   "Note that multiple edges and/or loops break the LGL software,\n"
   "but igraph does not check for this condition. Unless you know\n"
   "that the graph does not have multiple edges and/or loops, it\n"
   "is wise to call L{simplify()} before saving.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
   "@param names: the name of the vertex attribute containing the name\n"
   "  of the vertices. If you don't want to store vertex names,\n"
   "  supply C{None} here.\n"
   "@param weights: the name of the edge attribute containing the weight\n"
   "  of the vertices. If you don't want to store weights,\n"
   "  supply C{None} here.\n"
   "@param isolates: whether to include isolated vertices in the output.\n"},
  /* interface to igraph_write_graph_pajek */
  {"write_pajek", (PyCFunction) igraphmodule_Graph_write_pajek,
   METH_VARARGS | METH_KEYWORDS,
   "write_pajek(f)\n\n"
   "Writes the graph in Pajek format to the given file.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
   },
  /* interface to igraph_write_graph_edgelist */
  {"write_graphml", (PyCFunction) igraphmodule_Graph_write_graphml,
   METH_VARARGS | METH_KEYWORDS,
   "write_graphml(f)\n\n"
   "Writes the graph to a GraphML file.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
  },
  /* interface to igraph_write_graph_leda */
  {"write_leda", (PyCFunction) igraphmodule_Graph_write_leda,
   METH_VARARGS | METH_KEYWORDS,
   "write_leda(f, names=\"name\", weights=\"weights\")\n\n"
   "Writes the graph to a file in LEDA native format.\n\n"
   "The LEDA format supports at most one attribute per vertex and edge. You can\n"
   "specify which vertex and edge attribute you want to use. Note that the\n"
   "name of the attribute is not saved in the LEDA file.\n\n"
   "@param f: the name of the file to be written or a Python file handle\n"
   "@param names: the name of the vertex attribute to be stored along with\n"
   "  the vertices. It is usually used to store the vertex names (hence the\n"
   "  name of the keyword argument), but you may also use a numeric attribute.\n"
   "  If you don't want to store any vertex attributes, supply C{None} here.\n"
   "@param weights: the name of the edge attribute to be stored along with\n"
   "  the edges. It is usually used to store the edge weights (hence the\n"
   "  name of the keyword argument), but you may also use a string attribute.\n"
   "  If you don't want to store any edge attributes, supply C{None} here.\n"},

  /////////////////
  // ISOMORPHISM //
  /////////////////
  {"isoclass", (PyCFunction) igraphmodule_Graph_isoclass,
   METH_VARARGS | METH_KEYWORDS,
   "isoclass(vertices)\n\n"
   "Returns the isomorphy class of the graph or its subgraph.\n\n"
   "Isomorphy class calculations are implemented only for graphs with\n"
   "3 or 4 vertices.\n\n"
   "@param vertices: a list of vertices if we want to calculate the\n"
   "  isomorphy class for only a subset of vertices. C{None} means to\n"
   "  use the full graph.\n"
   "@return: the isomorphy class of the (sub)graph\n\n"},
  {"isomorphic", (PyCFunction) igraphmodule_Graph_isomorphic,
   METH_VARARGS | METH_KEYWORDS,
   "isomorphic(other)\n\n"
   "Checks whether the graph is isomorphic to another graph.\n\n"
   "The algorithm being used is selected using a simple heuristic:\n\n"
   "  - If one graph is directed and the other undirected, an exception\n"
   "    is thrown.\n\n"
   "  - If the two graphs does not have the same number of vertices and\n"
   "    edges, it returns with C{False}\n\n"
   "  - If the graphs have three or four vertices, then an O(1) algorithm\n"
   "    is used with precomputed data.\n\n"
   "  - Otherwise if the graphs are directed, then the VF2 isomorphism\n"
   "    algorithm is used (see L{Graph.isomorphic_vf2}).\n\n"
   "  - Otherwise the BLISS isomorphism algorithm is used, see\n"
   "    L{Graph.isomorphic_bliss}.\n\n"
   "@return: C{True} if the graphs are isomorphic, C{False} otherwise.\n"
  },
  {"isomorphic_bliss", (PyCFunction) igraphmodule_Graph_isomorphic_bliss,
   METH_VARARGS | METH_KEYWORDS,
   "isomorphic_bliss(other, return_mapping_12=False, return_mapping_21=False,\n"
   "  sh1=\"fm\", sh2=\"fm\")\n\n"
   "Checks whether the graph is isomorphic to another graph, using the\n"
   "BLISS isomorphism algorithm.\n\n"
   "@param other: the other graph with which we want to compare the graph.\n"
   "@param return_mapping_12: if C{True}, calculates the mapping which maps\n"
   "  the vertices of the first graph to the second.\n"
   "@param return_mapping_21: if C{True}, calculates the mapping which maps\n"
   "  the vertices of the second graph to the first.\n"
   "@param sh1: splitting heuristics for the first graph as a\n"
   "  case-insensitive string, with the following possible values:\n\n"
   "    - C{\"f\"}: first non-singleton cell\n\n"
   "    - C{\"fl\"}: first largest non-singleton cell\n\n"
   "    - C{\"fs\"}: first smallest non-singleton cell\n\n"
   "    - C{\"fm\"}: first maximally non-trivially connected non-singleton\n"
   "      cell\n\n"
   "    - C{\"flm\"}: largest maximally non-trivially connected\n"
   "      non-singleton cell\n\n"
   "    - C{\"fsm\"}: smallest maximally non-trivially connected\n"
   "      non-singleton cell\n\n"
   "@param sh2: splitting heuristics to be used for the second graph.\n"
   "  Accepted values are as above.\n"
   "@return: if no mapping is calculated, the result is C{True} if the graphs\n"
   "  are isomorphic, C{False} otherwise. If any or both mappings are\n"
   "  calculated, the result is a 3-tuple, the first element being the\n"
   "  above mentioned boolean, the second element being the 1 -> 2 mapping\n"
   "  and the third element being the 2 -> 1 mapping. If the corresponding\n"
   "  mapping was not calculated, C{None} is returned in the appropriate\n"
   "  element of the 3-tuple.\n"},
  {"isomorphic_vf2", (PyCFunction) igraphmodule_Graph_isomorphic_vf2,
   METH_VARARGS | METH_KEYWORDS,
   "isomorphic_vf2(other=None, color1=None, color2=None, edge_color1=None,\n"
   "  edge_color2=None, return_mapping_12=False, return_mapping_21=False,\n"
   "  node_compat_fn=None, edge_compat_fn=None, callback=None)\n\n"
   "Checks whether the graph is isomorphic to another graph, using the\n"
   "VF2 isomorphism algorithm.\n\n"
   "Vertex and edge colors may be used to restrict the isomorphisms, as only\n"
   "vertices and edges with the same color will be allowed to match each other.\n\n"
   "@param other: the other graph with which we want to compare the graph.\n"
   "  If C{None}, the automorphjisms of the graph will be tested.\n"
   "@param color1: optional vector storing the coloring of the vertices of\n"
   "  the first graph. If C{None}, all vertices have the same color.\n"
   "@param color2: optional vector storing the coloring of the vertices of\n"
   "  the second graph. If C{None}, all vertices have the same color.\n"
   "@param edge_color1: optional vector storing the coloring of the edges of\n"
   "  the first graph. If C{None}, all edges have the same color.\n"
   "@param edge_color2: optional vector storing the coloring of the edges of\n"
   "  the second graph. If C{None}, all edges have the same color.\n"
   "@param return_mapping_12: if C{True}, calculates the mapping which maps\n"
   "  the vertices of the first graph to the second.\n"
   "@param return_mapping_21: if C{True}, calculates the mapping which maps\n"
   "  the vertices of the second graph to the first.\n"
   "@param callback: if not C{None}, the isomorphism search will not stop at\n"
   "  the first match; it will call this callback function instead for every\n"
   "  isomorphism found. The callback function must accept four arguments:\n"
   "  the first graph, the second graph, a mapping from the nodes of the\n"
   "  first graph to the second, and a mapping from the nodes of the second\n"
   "  graph to the first. The function must return C{True} if the search\n"
   "  should continue or C{False} otherwise.\n"
   "@param node_compat_fn: a function that receives the two graphs and two\n"
   "  node indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the nodes given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on node-specific\n"
   "  criteria that are too complicated to be represented by node color\n"
   "  vectors (i.e. the C{color1} and C{color2} parameters). C{None} means\n"
   "  that every node is compatible with every other node.\n"
   "@param edge_compat_fn: a function that receives the two graphs and two\n"
   "  edge indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the edges given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on edge-specific\n"
   "  criteria that are too complicated to be represented by edge color\n"
   "  vectors (i.e. the C{edge_color1} and C{edge_color2} parameters). C{None}\n"
   "  means that every edge is compatible with every other node.\n"
   "@return: if no mapping is calculated, the result is C{True} if the graphs\n"
   "  are isomorphic, C{False} otherwise. If any or both mappings are\n"
   "  calculated, the result is a 3-tuple, the first element being the\n"
   "  above mentioned boolean, the second element being the 1 -> 2 mapping\n"
   "  and the third element being the 2 -> 1 mapping. If the corresponding\n"
   "  mapping was not calculated, C{None} is returned in the appropriate\n"
   "  element of the 3-tuple.\n"},
  {"count_isomorphisms_vf2",
   (PyCFunction) igraphmodule_Graph_count_isomorphisms_vf2,
   METH_VARARGS | METH_KEYWORDS,
   "count_isomorphisms_vf2(other=None, color1=None, color2=None, edge_color1=None,\n"
   "  edge_color2=None, node_compat_fn=None, edge_compat_fn=None)\n\n"
   "Determines the number of isomorphisms between the graph and another one\n\n"
   "Vertex and edge colors may be used to restrict the isomorphisms, as only\n"
   "vertices and edges with the same color will be allowed to match each other.\n\n"
   "@param other: the other graph. If C{None}, the number of automorphisms\n"
   "  will be returned.\n"
   "@param color1: optional vector storing the coloring of the vertices of\n"
   "  the first graph. If C{None}, all vertices have the same color.\n"
   "@param color2: optional vector storing the coloring of the vertices of\n"
   "  the second graph. If C{None}, all vertices have the same color.\n"
   "@param edge_color1: optional vector storing the coloring of the edges of\n"
   "  the first graph. If C{None}, all edges have the same color.\n"
   "@param edge_color2: optional vector storing the coloring of the edges of\n"
   "  the second graph. If C{None}, all edges have the same color.\n"
   "@param node_compat_fn: a function that receives the two graphs and two\n"
   "  node indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the nodes given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on node-specific\n"
   "  criteria that are too complicated to be represented by node color\n"
   "  vectors (i.e. the C{color1} and C{color2} parameters). C{None} means\n"
   "  that every node is compatible with every other node.\n"
   "@param edge_compat_fn: a function that receives the two graphs and two\n"
   "  edge indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the edges given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on edge-specific\n"
   "  criteria that are too complicated to be represented by edge color\n"
   "  vectors (i.e. the C{edge_color1} and C{edge_color2} parameters). C{None}\n"
   "  means that every edge is compatible with every other node.\n"
   "@return: the number of isomorphisms between the two given graphs (or the\n"
   "  number of automorphisms if C{other} is C{None}.\n"},
  {"get_isomorphisms_vf2", (PyCFunction) igraphmodule_Graph_get_isomorphisms_vf2,
   METH_VARARGS | METH_KEYWORDS,
   "get_isomorphisms_vf2(other=None, color1=None, color2=None, edge_color1=None,\n"
   "  edge_color2=None, node_compat_fn=None, edge_compat_fn=None)\n\n"
   "Returns all isomorphisms between the graph and another one\n\n"
   "Vertex and edge colors may be used to restrict the isomorphisms, as only\n"
   "vertices and edges with the same color will be allowed to match each other.\n\n"
   "@param other: the other graph. If C{None}, the automorphisms\n"
   "  will be returned.\n"
   "@param color1: optional vector storing the coloring of the vertices of\n"
   "  the first graph. If C{None}, all vertices have the same color.\n"
   "@param color2: optional vector storing the coloring of the vertices of\n"
   "  the second graph. If C{None}, all vertices have the same color.\n"
   "@param edge_color1: optional vector storing the coloring of the edges of\n"
   "  the first graph. If C{None}, all edges have the same color.\n"
   "@param edge_color2: optional vector storing the coloring of the edges of\n"
   "  the second graph. If C{None}, all edges have the same color.\n"
   "@param node_compat_fn: a function that receives the two graphs and two\n"
   "  node indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the nodes given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on node-specific\n"
   "  criteria that are too complicated to be represented by node color\n"
   "  vectors (i.e. the C{color1} and C{color2} parameters). C{None} means\n"
   "  that every node is compatible with every other node.\n"
   "@param edge_compat_fn: a function that receives the two graphs and two\n"
   "  edge indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the edges given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on edge-specific\n"
   "  criteria that are too complicated to be represented by edge color\n"
   "  vectors (i.e. the C{edge_color1} and C{edge_color2} parameters). C{None}\n"
   "  means that every edge is compatible with every other node.\n"
   "@return: a list of lists, each item of the list containing the mapping\n"
   "  from vertices of the second graph to the vertices of the first one\n"},

  {"subisomorphic_vf2", (PyCFunction) igraphmodule_Graph_subisomorphic_vf2,
   METH_VARARGS | METH_KEYWORDS,
   "subisomorphic_vf2(other, color1=None, color2=None, edge_color1=None,\n"
   "  edge_color2=None, return_mapping_12=False, return_mapping_21=False,\n"
   "  callback=None, node_compat_fn=None, edge_compat_fn=None)\n\n"
   "Checks whether a subgraph of the graph is isomorphic to another graph.\n\n"
   "Vertex and edge colors may be used to restrict the isomorphisms, as only\n"
   "vertices and edges with the same color will be allowed to match each other.\n\n"
   "@param other: the other graph with which we want to compare the graph.\n"
   "@param color1: optional vector storing the coloring of the vertices of\n"
   "  the first graph. If C{None}, all vertices have the same color.\n"
   "@param color2: optional vector storing the coloring of the vertices of\n"
   "  the second graph. If C{None}, all vertices have the same color.\n"
   "@param edge_color1: optional vector storing the coloring of the edges of\n"
   "  the first graph. If C{None}, all edges have the same color.\n"
   "@param edge_color2: optional vector storing the coloring of the edges of\n"
   "  the second graph. If C{None}, all edges have the same color.\n"
   "@param return_mapping_12: if C{True}, calculates the mapping which maps\n"
   "  the vertices of the first graph to the second. The mapping can contain\n"
   "  -1 if a given node is not mapped.\n"
   "@param return_mapping_21: if C{True}, calculates the mapping which maps\n"
   "  the vertices of the second graph to the first. The mapping can contain\n"
   "  -1 if a given node is not mapped.\n"
   "@param callback: if not C{None}, the subisomorphism search will not stop at\n"
   "  the first match; it will call this callback function instead for every\n"
   "  subisomorphism found. The callback function must accept four arguments:\n"
   "  the first graph, the second graph, a mapping from the nodes of the\n"
   "  first graph to the second, and a mapping from the nodes of the second\n"
   "  graph to the first. The function must return C{True} if the search\n"
   "  should continue or C{False} otherwise.\n"
   "@param node_compat_fn: a function that receives the two graphs and two\n"
   "  node indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the nodes given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on node-specific\n"
   "  criteria that are too complicated to be represented by node color\n"
   "  vectors (i.e. the C{color1} and C{color2} parameters). C{None} means\n"
   "  that every node is compatible with every other node.\n"
   "@param edge_compat_fn: a function that receives the two graphs and two\n"
   "  edge indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the edges given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on edge-specific\n"
   "  criteria that are too complicated to be represented by edge color\n"
   "  vectors (i.e. the C{edge_color1} and C{edge_color2} parameters). C{None}\n"
   "  means that every edge is compatible with every other node.\n"
   "@return: if no mapping is calculated, the result is C{True} if the graph\n"
   "  contains a subgraph that's isomorphic to the given one, C{False}\n"
   "  otherwise. If any or both mappings are calculated, the result is a\n"
   "  3-tuple, the first element being the above mentioned boolean, the\n"
   "  second element being the 1 -> 2 mapping and the third element being\n"
   "  the 2 -> 1 mapping. If the corresponding mapping was not calculated,\n"
   "  C{None} is returned in the appropriate element of the 3-tuple.\n"},
  {"count_subisomorphisms_vf2",
   (PyCFunction) igraphmodule_Graph_count_subisomorphisms_vf2,
   METH_VARARGS | METH_KEYWORDS,
   "count_subisomorphisms_vf2(other, color1=None, color2=None,\n"
   "  edge_color1=None, edge_color2=None, node_compat_fn=None,\n"
   "  edge_compat_fn=None)\n\n"
   "Determines the number of subisomorphisms between the graph and another one\n\n"
   "Vertex and edge colors may be used to restrict the isomorphisms, as only\n"
   "vertices and edges with the same color will be allowed to match each other.\n\n"
   "@param other: the other graph.\n"
   "@param color1: optional vector storing the coloring of the vertices of\n"
   "  the first graph. If C{None}, all vertices have the same color.\n"
   "@param color2: optional vector storing the coloring of the vertices of\n"
   "  the second graph. If C{None}, all vertices have the same color.\n"
   "@param edge_color1: optional vector storing the coloring of the edges of\n"
   "  the first graph. If C{None}, all edges have the same color.\n"
   "@param edge_color2: optional vector storing the coloring of the edges of\n"
   "  the second graph. If C{None}, all edges have the same color.\n"
   "@param node_compat_fn: a function that receives the two graphs and two\n"
   "  node indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the nodes given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on node-specific\n"
   "  criteria that are too complicated to be represented by node color\n"
   "  vectors (i.e. the C{color1} and C{color2} parameters). C{None} means\n"
   "  that every node is compatible with every other node.\n"
   "@param edge_compat_fn: a function that receives the two graphs and two\n"
   "  edge indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the edges given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on edge-specific\n"
   "  criteria that are too complicated to be represented by edge color\n"
   "  vectors (i.e. the C{edge_color1} and C{edge_color2} parameters). C{None}\n"
   "  means that every edge is compatible with every other node.\n"
   "@return: the number of subisomorphisms between the two given graphs\n"},
  {"get_subisomorphisms_vf2",
   (PyCFunction) igraphmodule_Graph_get_subisomorphisms_vf2,
   METH_VARARGS | METH_KEYWORDS,
   "get_subisomorphisms_vf2(other, color1=None, color2=None,\n"
   "  edge_color1=None, edge_color2=None, node_compat_fn=None,\n"
   "  edge_compat_fn=None)\n\n"
   "Returns all subisomorphisms between the graph and another one\n\n"
   "Vertex and edge colors may be used to restrict the isomorphisms, as only\n"
   "vertices and edges with the same color will be allowed to match each other.\n\n"
   "@param other: the other graph.\n"
   "@param color1: optional vector storing the coloring of the vertices of\n"
   "  the first graph. If C{None}, all vertices have the same color.\n"
   "@param color2: optional vector storing the coloring of the vertices of\n"
   "  the second graph. If C{None}, all vertices have the same color.\n"
   "@param edge_color1: optional vector storing the coloring of the edges of\n"
   "  the first graph. If C{None}, all edges have the same color.\n"
   "@param edge_color2: optional vector storing the coloring of the edges of\n"
   "  the second graph. If C{None}, all edges have the same color.\n"
   "@param node_compat_fn: a function that receives the two graphs and two\n"
   "  node indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the nodes given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on node-specific\n"
   "  criteria that are too complicated to be represented by node color\n"
   "  vectors (i.e. the C{color1} and C{color2} parameters). C{None} means\n"
   "  that every node is compatible with every other node.\n"
   "@param edge_compat_fn: a function that receives the two graphs and two\n"
   "  edge indices (one from the first graph, one from the second graph) and\n"
   "  returns C{True} if the edges given by the two indices are compatible\n"
   "  (i.e. they could be matched to each other) or C{False} otherwise. This\n"
   "  can be used to restrict the set of isomorphisms based on edge-specific\n"
   "  criteria that are too complicated to be represented by edge color\n"
   "  vectors (i.e. the C{edge_color1} and C{edge_color2} parameters). C{None}\n"
   "  means that every edge is compatible with every other node.\n"
   "@return: a list of lists, each item of the list containing the mapping\n"
   "  from vertices of the second graph to the vertices of the first one\n"},

  ////////////////////////
  // ATTRIBUTE HANDLING //
  ////////////////////////
  {"attributes", (PyCFunction) igraphmodule_Graph_attributes,
   METH_NOARGS,
   "attributes()\n\n" "@return: the attribute name list of the graph\n"},
  {"vertex_attributes", (PyCFunction) igraphmodule_Graph_vertex_attributes,
   METH_NOARGS,
   "vertex_attributes()\n\n"
   "@return: the attribute name list of the graph's vertices\n"},
  {"edge_attributes", (PyCFunction) igraphmodule_Graph_edge_attributes,
   METH_NOARGS,
   "edge_attributes()\n\n"
   "@return: the attribute name list of the graph's edges\n"},

  ///////////////
  // OPERATORS //
  ///////////////
  {"complementer", (PyCFunction) igraphmodule_Graph_complementer,
   METH_VARARGS,
   "complementer(loops=False)\n\n"
   "Returns the complementer of the graph\n\n"
   "@param loops: whether to include loop edges in the complementer.\n"
   "@return: the complementer of the graph\n"},
  {"compose", (PyCFunction) igraphmodule_Graph_compose,
   METH_O, "compose(other)\n\nReturns the composition of two graphs."},
  {"difference", (PyCFunction) igraphmodule_Graph_difference,
   METH_O,
   "difference(other)\n\nSubtracts the given graph from the original"},
  {"disjoint_union", (PyCFunction) igraphmodule_Graph_disjoint_union,
   METH_O,
   "disjoint_union(graphs)\n\n"
   "Creates the disjoint union of two (or more) graphs.\n\n"
   "@param graphs: the list of graphs to be united with the current one.\n"},
  {"intersection", (PyCFunction) igraphmodule_Graph_intersection,
   METH_O,
   "intersection(graphs)\n\n"
   "Creates the intersection of two (or more) graphs.\n\n"
   "@param graphs: the list of graphs to be intersected with\n"
   "  the current one.\n"},
  {"union", (PyCFunction) igraphmodule_Graph_union,
   METH_O,
   "union(graphs)\n\n"
   "Creates the union of two (or more) graphs.\n\n"
   "@param graphs: the list of graphs to be united with\n"
   "  the current one.\n"},

  /*****************/
  /* MAXIMUM FLOWS */
  /*****************/
  {"maxflow_value", (PyCFunction) igraphmodule_Graph_maxflow_value,
   METH_VARARGS | METH_KEYWORDS,
   "maxflow_value(source, target, capacity=None)\n\n"
   "Returns the value of the maximum flow between the source and target vertices.\n\n"
   "@param source: the source vertex ID\n"
   "@param target: the target vertex ID\n"
   "@param capacity: the capacity of the edges. It must be a list or a valid\n"
   "  attribute name or C{None}. In the latter case, every edge will have the\n"
   "  same capacity.\n"
   "@return: the value of the maximum flow between the given vertices\n"},

  {"maxflow", (PyCFunction) igraphmodule_Graph_maxflow,
   METH_VARARGS | METH_KEYWORDS,
   "maxflow(source, target, capacity=None)\n\n"
   "Returns the maximum flow between the source and target vertices.\n\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a L{Flow} object. It is advised\n"
   "  to use that.\n"
   "@param source: the source vertex ID\n"
   "@param target: the target vertex ID\n"
   "@param capacity: the capacity of the edges. It must be a list or a valid\n"
   "  attribute name or C{None}. In the latter case, every edge will have the\n"
   "  same capacity.\n"
   "@return: a tuple containing the following: the value of the maximum flow\n"
   "  between the given vertices, the flow value on all the edges, the edge\n"
   "  IDs that are part of the corresponding minimum cut, and the vertex IDs\n"
   "  on one side of the cut. For directed graphs, the flow value vector gives\n"
   "  the flow value on each edge. For undirected graphs, the flow value is\n"
   "  positive if the flow goes from the smaller vertex ID to the bigger one\n"
   "  and negative if the flow goes from the bigger vertex ID to the smaller."
  },

  /**********************/
  /* CUTS, MINIMUM CUTS */
  /**********************/
  {"all_st_cuts", (PyCFunction) igraphmodule_Graph_all_st_cuts,
   METH_VARARGS | METH_KEYWORDS,
   "all_st_cuts(source, target)\n\n"
   "Returns all the cuts between the source and target vertices in a\n"
   "directed graph.\n\n"
   "This function lists all edge-cuts between a source and a target vertex.\n"
   "Every cut is listed exactly once.\n\n"
   "@param source: the source vertex ID\n"
   "@param target: the target vertex ID\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a list of L{Cut} objects. It is\n"
   "  advised to use that.\n"
   "@return: a tuple where the first element is a list of lists of edge IDs\n"
   "  representing a cut and the second element is a list of lists of vertex\n"
   "  IDs representing the sets of vertices that were separated by the cuts.\n"
  },
  {"all_st_mincuts", (PyCFunction) igraphmodule_Graph_all_st_mincuts,
   METH_VARARGS | METH_KEYWORDS,
   "all_st_mincuts(source, target)\n\n"
   "Returns all minimum cuts between the source and target vertices in a\n"
   "directed graph.\n\n"
   "@param source: the source vertex ID\n"
   "@param target: the target vertex ID\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a list of L{Cut} objects. It is\n"
   "  advised to use that.\n"
  },
  {"mincut_value", (PyCFunction) igraphmodule_Graph_mincut_value,
   METH_VARARGS | METH_KEYWORDS,
   "mincut_value(source=-1, target=-1, capacity=None)\n\n"
   "Returns the minimum cut between the source and target vertices.\n\n"
   "@param source: the source vertex ID. If negative, the calculation is\n"
   "  done for every vertex except the target and the minimum is returned.\n"
   "@param target: the target vertex ID. If negative, the calculation is\n"
   "  done for every vertex except the source and the minimum is returned.\n"
   "@param capacity: the capacity of the edges. It must be a list or a valid\n"
   "  attribute name or C{None}. In the latter case, every edge will have the\n"
   "  same capacity.\n"
   "@return: the value of the minimum cut between the given vertices\n"},

  {"mincut", (PyCFunction) igraphmodule_Graph_mincut,
   METH_VARARGS | METH_KEYWORDS,
   "mincut(capacity=None)\n\n"
   "Calculates the minimum cut in a graph.\n\n"
   "The minimum cut is the minimum set of edges which needs to be removed\n"
   "to disconnect the graph. The minimum is calculated using the weights\n"
   "(capacities) of the edges, so the cut with the minimum total capacity\n"
   "is calculated.\n"
   "For undirected graphs, it uses the Stoer-Wagner algorithm as described\n"
   "in the reference given below.\n\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a L{Cut} object. It is advised\n"
   "  to use that.\n"
   "@return: the value of the minimum cut, the IDs of vertices in the\n"
   "  first and second partition, and the IDs of edges in the cut,\n"
   "  packed in a 4-tuple\n\n"
   "@newfield ref: Reference\n"
   "@ref: M. Stoer, F. Wagner: A simple min-cut algorithm. Journal of\n"
   "  the ACM 44(4):585-591, 1997.\n"
   },

  /*********************/
  /* VERTEX SEPARATORS */
  /*********************/
  {"all_minimal_st_separators",
   (PyCFunction) igraphmodule_Graph_all_minimal_st_separators, METH_NOARGS,
   "all_minimal_st_separators()\n\n"
   "Returns a list containing all the minimal s-t separators of a graph.\n\n"
   "A minimal separator is a set of vertices whose removal disconnects the graph,\n"
   "while the removal of any subset of the set keeps the graph connected.\n\n"
   "@return: a list where each item lists the vertex indices of a given\n"
   "  minimal s-t separator.\n"
   "@newfield ref: Reference\n"
   "@ref: Anne Berry, Jean-Paul Bordat and Olivier Cogis: Generating all the\n"
   "  minimal separators of a graph. In: Peter Widmayer, Gabriele Neyer and\n"
   "  Stephan Eidenbenz (eds.): Graph-theoretic concepts in computer science,\n"
   "  1665, 167--172, 1999. Springer.\n"},

  {"is_minimal_separator", (PyCFunction) igraphmodule_Graph_is_minimal_separator,
   METH_VARARGS | METH_KEYWORDS,
   "is_minimal_separator(vertices)\n\n"
   "Decides whether the given vertex set is a minimal separator.\n\n"
   "A minimal separator is a set of vertices whose removal disconnects the graph,\n"
   "while the removal of any subset of the set keeps the graph connected.\n\n"
   "@param vertices: a single vertex ID or a list of vertex IDs\n"
   "@return: C{True} is the given vertex set is a minimal separator, C{False}\n"
   "  otherwise.\n"},

  {"is_separator", (PyCFunction) igraphmodule_Graph_is_separator,
   METH_VARARGS | METH_KEYWORDS,
   "is_separator(vertices)\n\n"
   "Decides whether the removal of the given vertices disconnects the graph.\n\n"
   "@param vertices: a single vertex ID or a list of vertex IDs\n"
   "@return: C{True} is the given vertex set is a separator, C{False} if not.\n"},

  {"minimum_size_separators",
   (PyCFunction) igraphmodule_Graph_minimum_size_separators, METH_NOARGS,
   "minimum_size_separators()\n\n"
   "Returns a list containing all separator vertex sets of minimum size.\n\n"
   "A vertex set is a separator if its removal disconnects the graph. This method\n"
   "lists all the separators for which no smaller separator set exists in the\n"
   "given graph.\n\n"
   "@return: a list where each item lists the vertex indices of a given\n"
   "  separator of minimum size.\n"
   "@newfield ref: Reference\n"
   "@ref: Arkady Kanevsky: Finding all minimum-size separating vertex sets\n"
   "  in a graph. Networks 23:533--541, 1993.\n"},

  /*******************/
  /* COHESIVE BLOCKS */
  /*******************/
  {"cohesive_blocks", (PyCFunction) igraphmodule_Graph_cohesive_blocks,
   METH_NOARGS,
   "cohesive_blocks()\n\n"
   "Calculates the cohesive block structure of the graph.\n\n"
   "@attention: this function has a more convenient interface in class\n"
   "  L{Graph} which wraps the result in a L{CohesiveBlocks} object.\n"
   "  It is advised to use that.\n"
  },

  /********************************/
  /* CLIQUES AND INDEPENDENT SETS */
  /********************************/
  {"cliques", (PyCFunction) igraphmodule_Graph_cliques,
   METH_VARARGS | METH_KEYWORDS,
   "cliques(min=0, max=0)\n\n"
   "Returns some or all cliques of the graph as a list of tuples.\n\n"
   "A clique is a complete subgraph -- a set of vertices where an edge\n"
   "is present between any two of them (excluding loops)\n\n"
   "@param min: the minimum size of cliques to be returned. If zero or\n"
   "  negative, no lower bound will be used.\n"
   "@param max: the maximum size of cliques to be returned. If zero or\n"
   "  negative, no upper bound will be used."},
  {"largest_cliques", (PyCFunction) igraphmodule_Graph_largest_cliques,
   METH_NOARGS,
   "largest_cliques()\n\n"
   "Returns the largest cliques of the graph as a list of tuples.\n\n"
   "Quite intuitively a clique is considered largest if there is no clique\n"
   "with more vertices in the whole graph. All largest cliques are maximal\n"
   "(i.e. nonextendable) but not all maximal cliques are largest.\n\n"
   "@see: L{clique_number()} for the size of the largest cliques or\n"
   "  L{maximal_cliques()} for the maximal cliques"},
  {"maximal_cliques", (PyCFunction) igraphmodule_Graph_maximal_cliques,
   METH_VARARGS | METH_KEYWORDS,
   "maximal_cliques(min=0, max=0)\n\n"
   "Returns the maximal cliques of the graph as a list of tuples.\n\n"
   "A maximal clique is a clique which can't be extended by adding any other\n"
   "vertex to it. A maximal clique is not necessarily one of the largest\n"
   "cliques in the graph.\n\n"
   "@param min: the minimum size of maximal cliques to be returned. If zero\n"
   "  or negative, no lower bound will be used.\n\n"
   "@param max: the maximum size of maximal cliques to be returned. If zero\n"
   "  or negative, no upper bound will be used. If nonzero, the size of every\n"
   "  maximal clique found will be compared to this value and a clique will\n"
   "  be returned only if its size is smaller than this limit.\n\n"
   "@see: L{largest_cliques()} for the largest cliques."},
  {"clique_number", (PyCFunction) igraphmodule_Graph_clique_number,
   METH_NOARGS,
   "clique_number()\n\n"
   "Returns the clique number of the graph.\n\n"
   "The clique number of the graph is the size of the largest clique.\n\n"
   "@see: L{largest_cliques()} for the largest cliques."},
  {"independent_vertex_sets",
   (PyCFunction) igraphmodule_Graph_independent_vertex_sets,
   METH_VARARGS | METH_KEYWORDS,
   "independent_vertex_sets(min=0, max=0)\n\n"
   "Returns some or all independent vertex sets of the graph as a list of tuples.\n\n"
   "Two vertices are independent if there is no edge between them. Members\n"
   "of an independent vertex set are mutually independent.\n\n"
   "@param min: the minimum size of sets to be returned. If zero or\n"
   "  negative, no lower bound will be used.\n"
   "@param max: the maximum size of sets to be returned. If zero or\n"
   "  negative, no upper bound will be used."},
  {"largest_independent_vertex_sets",
   (PyCFunction) igraphmodule_Graph_largest_independent_vertex_sets,
   METH_NOARGS,
   "largest_independent_vertex_sets()\n\n"
   "Returns the largest independent vertex sets of the graph as a list of tuples.\n\n"
   "Quite intuitively an independent vertex set is considered largest if\n"
   "there is no other set with more vertices in the whole graph. All largest\n"
   "sets are maximal (i.e. nonextendable) but not all maximal sets\n"
   "are largest.\n\n"
   "@see: L{independence_number()} for the size of the largest independent\n"
   "  vertex sets or L{maximal_independent_vertex_sets()} for the maximal\n"
   "  (nonextendable) independent vertex sets"},
  {"maximal_independent_vertex_sets",
   (PyCFunction) igraphmodule_Graph_maximal_independent_vertex_sets,
   METH_NOARGS,
   "maximal_independent_vertex_sets()\n\n"
   "Returns the maximal independent vertex sets of the graph as a list of tuples.\n\n"
   "A maximal independent vertex set is an independent vertex set\n"
   "which can't be extended by adding any other vertex to it. A maximal\n"
   "independent vertex set is not necessarily one of the largest\n"
   "independent vertex sets in the graph.\n\n"
   "@see: L{largest_independent_vertex_sets()} for the largest independent\n"
   "  vertex sets\n\n"
   "@newfield ref: Reference\n"
   "@ref: S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka: I{A new\n"
   "  algorithm for generating all the maximal independent sets}.\n"
   "  SIAM J Computing, 6:505--517, 1977."},
  {"independence_number",
   (PyCFunction) igraphmodule_Graph_independence_number,
   METH_NOARGS,
   "independence_number()\n\n"
   "Returns the independence number of the graph.\n\n"
   "The independence number of the graph is the size of the largest\n"
   "independent vertex set.\n\n"
   "@see: L{largest_independent_vertex_sets()} for the largest independent\n"
   "  vertex sets"},

  /*********************************/
  /* COMMUNITIES AND DECOMPOSITION */
  /*********************************/
  {"modularity", (PyCFunction) igraphmodule_Graph_modularity,
   METH_VARARGS | METH_KEYWORDS,
   "modularity(membership, weights=None)\n\n"
   "Calculates the modularity of the graph with respect to some vertex types.\n\n"
   "The modularity of a graph w.r.t. some division measures how good the\n"
   "division is, or how separated are the different vertex types from each\n"
   "other. It is defined as M{Q=1/(2m) * sum(Aij-ki*kj/(2m)delta(ci,cj),i,j)}.\n"
   "M{m} is the number of edges, M{Aij} is the element of the M{A} adjacency\n"
   "matrix in row M{i} and column M{j}, M{ki} is the degree of node M{i},\n"
   "M{kj} is the degree of node M{j}, and M{Ci} and C{cj} are the types of\n"
   "the two vertices (M{i} and M{j}). M{delta(x,y)} is one iff M{x=y}, 0\n"
   "otherwise.\n\n"
   "If edge weights are given, the definition of modularity is modified as\n"
   "follows: M{Aij} becomes the weight of the corresponding edge, M{ki}\n"
   "is the total weight of edges incident on vertex M{i}, M{kj} is the\n"
   "total weight of edges incident on vertex M{j} and M{m} is the total\n"
   "edge weight in the graph.\n\n"
   "@attention: method overridden in L{Graph} to allow L{VertexClustering}\n"
   "  objects as a parameter. This method is not strictly necessary, since\n"
   "  the L{VertexClustering} class provides a variable called C{modularity}.\n"
   "@param membership: the membership vector, e.g. the vertex type index for\n"
   "  each vertex.\n"
   "@param weights: optional edge weights or C{None} if all edges are weighed\n"
   "  equally.\n"
   "@return: the modularity score. Score larger than 0.3 usually indicates\n"
   "  strong community structure.\n"
   "@newfield ref: Reference\n"
   "@ref: MEJ Newman and M Girvan: Finding and evaluating community structure\n"
   "  in networks. Phys Rev E 69 026113, 2004.\n"
  },
  {"coreness", (PyCFunction) igraphmodule_Graph_coreness,
   METH_VARARGS | METH_KEYWORDS,
   "coreness(mode=ALL)\n\n"
   "Finds the coreness (shell index) of the vertices of the network.\n\n"
   "The M{k}-core of a graph is a maximal subgraph in which each vertex\n"
   "has at least degree k. (Degree here means the degree in the\n"
   "subgraph of course). The coreness of a vertex is M{k} if it\n"
   "is a member of the M{k}-core but not a member of the M{k+1}-core.\n\n"
   "@param mode: whether to compute the in-corenesses (L{IN}), the\n"
   "  out-corenesses (L{OUT}) or the undirected corenesses (L{ALL}).\n"
   "  Ignored and assumed to be L{ALL} for undirected graphs.\n"
   "@return: the corenesses for each vertex.\n\n"
   "@newfield ref: Reference\n"
   "@ref: Vladimir Batagelj, Matjaz Zaversnik: I{An M{O(m)} Algorithm\n"
   "  for Core Decomposition of Networks.}"},
  {"community_fastgreedy",
   (PyCFunction) igraphmodule_Graph_community_fastgreedy,
   METH_VARARGS | METH_KEYWORDS,
   "community_fastgreedy(weights=None)\n\n"
   "Finds the community structure of the graph according to the algorithm of\n"
   "Clauset et al based on the greedy optimization of modularity.\n\n"
   "This is a bottom-up algorithm: initially every vertex belongs to a separate\n"
   "community, and communities are merged one by one. In every step, the two\n"
   "communities being merged are the ones which result in the maximal increase\n"
   "in modularity.\n\n"
   "@attention: this function is wrapped in a more convenient syntax in the\n"
   "  derived class L{Graph}. It is advised to use that instead of this version.\n\n"
   "@param weights: name of an edge attribute or a list containing\n"
   "  edge weights\n"
   "@return: a tuple with the following elements:\n"
   "  1. The list of merges\n"
   "  2. The modularity scores before each merge\n"
   "\n"
   "@newfield ref: Reference\n"
   "@ref: A. Clauset, M. E. J. Newman and C. Moore: I{Finding community\n"
   "  structure in very large networks.} Phys Rev E 70, 066111 (2004).\n"
   "@see: modularity()\n"
  },
  {"community_infomap",
   (PyCFunction) igraphmodule_Graph_community_infomap,
   METH_VARARGS | METH_KEYWORDS,
   "community_infomap(edge_weights=None, vertex_weights=None, trials=10)\n\n"
   "Finds the community structure of the network according to the Infomap\n"
   "method of Martin Rosvall and Carl T. Bergstrom.\n\n"
   "See U{http://www.mapequation.org} for a visualization of the algorithm\n"
   "or one of the references provided below.\n\n"
   "@param edge_weights: name of an edge attribute or a list containing\n"
   "  edge weights.\n"
   "@param vertex_weights: name of an vertex attribute or a list containing\n"
   "  vertex weights.\n"
   "@param trials: the number of attempts to partition the network.\n"
   "@return: the calculated membership vector and the corresponding\n"
   "  codelength in a tuple.\n"
   "\n"
   "@newfield ref: Reference\n"
   "@ref: M. Rosvall and C. T. Bergstrom: I{Maps of information flow reveal\n"
   "  community structure in complex networks}. PNAS 105, 1118 (2008).\n"
   "  U{http://arxiv.org/abs/0707.0609}\n"
   "@ref: M. Rosvall, D. Axelsson and C. T. Bergstrom: I{The map equation}.\n"
   "  Eur Phys J Special Topics 178, 13 (2009). U{http://arxiv.org/abs/0906.1405}\n"
  },
  {"community_label_propagation",
   (PyCFunction) igraphmodule_Graph_community_label_propagation,
   METH_VARARGS | METH_KEYWORDS,
   "community_label_propagation(weights=None, initial=None, fixed=None)\n\n"
   "Finds the community structure of the graph according to the label\n"
   "propagation method of Raghavan et al.\n\n"
   "Initially, each vertex is assigned a different label. After that,\n"
   "each vertex chooses the dominant label in its neighbourhood in each\n"
   "iteration. Ties are broken randomly and the order in which the\n"
   "vertices are updated is randomized before every iteration. The algorithm\n"
   "ends when vertices reach a consensus.\n\n"
   "Note that since ties are broken randomly, there is no guarantee that\n"
   "the algorithm returns the same community structure after each run.\n"
   "In fact, they frequently differ. See the paper of Raghavan et al\n"
   "on how to come up with an aggregated community structure.\n\n"
   "@param weights: name of an edge attribute or a list containing\n"
   "  edge weights\n"
   "@param initial: name of a vertex attribute or a list containing\n"
   "  the initial vertex labels. Labels are identified by integers from\n"
   "  zero to M{n-1} where M{n} is the number of vertices. Negative\n"
   "  numbers may also be present in this vector, they represent unlabeled\n"
   "  vertices.\n"
   "@param fixed: a list of booleans for each vertex. C{True} corresponds\n"
   "  to vertices whose labeling should not change during the algorithm.\n"
   "  It only makes sense if initial labels are also given. Unlabeled\n"
   "  vertices cannot be fixed. Note that vertex attribute names are not\n"
   "  accepted here.\n"
   "@return: the resulting membership vector\n"
   "\n"
   "@newfield ref: Reference\n"
   "@ref: Raghavan, U.N. and Albert, R. and Kumara, S. Near linear\n"
   "  time algorithm to detect community structures in large-scale\n"
   "  networks. Phys Rev E 76:036106, 2007. U{http://arxiv.org/abs/0709.2938}.\n"
  },
  {"community_leading_eigenvector", (PyCFunction) igraphmodule_Graph_community_leading_eigenvector,
   METH_VARARGS | METH_KEYWORDS,
   "community_leading_eigenvector(n=-1)\n\n"
   "A proper implementation of Newman's eigenvector community structure\n"
   "detection. Each split is done by maximizing the modularity regarding\n"
   "the original network. See the reference for details.\n\n"
   "@attention: this function is wrapped in a more convenient syntax in the\n"
   "  derived class L{Graph}. It is advised to use that instead of this version.\n\n"
   "@param n: the desired number of communities. If negative, the algorithm\n"
   "  tries to do as many splits as possible. Note that the algorithm\n"
   "  won't split a community further if the signs of the leading eigenvector\n"
   "  are all the same.\n"
   "@return: a tuple where the first element is the membership vector of the\n"
   "  clustering and the second element is the merge matrix.\n\n"
   "@newfield ref: Reference\n"
  "@ref: MEJ Newman: Finding community structure in networks using the\n"
  "  eigenvectors of matrices, arXiv:physics/0605087\n"
  },
  {"community_multilevel",
   (PyCFunction) igraphmodule_Graph_community_multilevel,
   METH_VARARGS | METH_KEYWORDS,
   "community_multilevel(weights=None, return_levels=True)\n\n"
   "Finds the community structure of the graph according to the multilevel\n"
   "algorithm of Blondel et al. This is a bottom-up algorithm: initially\n"
   "every vertex belongs to a separate community, and vertices are moved\n"
   "between communities iteratively in a way that maximizes the vertices'\n"
   "local contribution to the overall modularity score. When a consensus is\n"
   "reached (i.e. no single move would increase the modularity score), every\n"
   "community in the original graph is shrank to a single vertex (while\n"
   "keeping the total weight of the incident edges) and the process continues\n"
   "on the next level. The algorithm stops when it is not possible to increase\n"
   "the modularity any more after shrinking the communities to vertices.\n\n"
   "@attention: this function is wrapped in a more convenient syntax in the\n"
   "  derived class L{Graph}. It is advised to use that instead of this version.\n\n"
   "@param weights: name of an edge attribute or a list containing\n"
   "  edge weights\n"
   "@param return_levels: if C{True}, returns the multilevel result. If\n"
   "  C{False}, only the best level (corresponding to the best modularity)\n"
   "  is returned.\n"
   "@return: either a single list describing the community membership of each\n"
   "  vertex (if C{return_levels} is C{False}), or a list of community membership\n"
   "  vectors, one corresponding to each level and a list of corresponding\n"
   "  modularities (if C{return_levels} is C{True}).\n"
   "\n"
   "@newfield ref: Reference\n"
   "@ref: VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast\n"
   "  unfolding of community hierarchies in large networks. J Stat Mech\n"
   "  P10008 (2008), http://arxiv.org/abs/0803.0476\n"
   "@see: modularity()\n"
  },
  {"community_edge_betweenness",
  (PyCFunction)igraphmodule_Graph_community_edge_betweenness,
  METH_VARARGS | METH_KEYWORDS,
  "community_edge_betweenness(directed=True, weights=None)\n\n"
  "Community structure detection based on the betweenness of the edges in\n"
  "the network. This algorithm was invented by M Girvan and MEJ Newman,\n"
  "see: M Girvan and MEJ Newman: Community structure in social and biological\n"
  "networks, Proc. Nat. Acad. Sci. USA 99, 7821-7826 (2002).\n\n"
  "The idea is that the betweenness of the edges connecting two communities\n"
  "is typically high. So we gradually remove the edge with the highest\n"
  "betweenness from the network and recalculate edge betweenness after every\n"
  "removal, as long as all edges are removed.\n\n"
  "@attention: this function is wrapped in a more convenient syntax in the\n"
  "  derived class L{Graph}. It is advised to use that instead of this version.\n\n"
  "@param directed: whether to take into account the directedness of the edges\n"
  "  when we calculate the betweenness values.\n"
  "@param weights: name of an edge attribute or a list containing\n"
  "  edge weights.\n\n"
  "@return: a tuple with the merge matrix that describes the dendrogram\n"
  "  and the modularity scores before each merge. The modularity scores\n"
  "  use the weights if the original graph was weighted.\n"
  },
  {"community_optimal_modularity",
   (PyCFunction) igraphmodule_Graph_community_optimal_modularity,
   METH_VARARGS | METH_KEYWORDS,
   "community_optimal_modularity(verbose=False)\n\n"
   "Calculates the optimal modularity score of the graph and the\n"
   "corresponding community structure.\n\n"
   "This function uses the GNU Linear Programming Kit to solve a large\n"
   "integer optimization problem in order to find the optimal modularity\n"
   "score and the corresponding community structure, therefore it is\n"
   "unlikely to work for graphs larger than a few (less than a hundred)\n"
   "vertices. Consider using one of the heuristic approaches instead if\n"
   "you have such a large graph.\n\n"
   "@return: the calculated membership vector and the corresponding\n"
   "  modularity in a tuple.\n"
  },
  {"community_spinglass",
   (PyCFunction) igraphmodule_Graph_community_spinglass,
   METH_VARARGS | METH_KEYWORDS,
   "community_spinglass(weights=None, spins=25, parupdate=False, "
   "start_temp=1, stop_temp=0.01, cool_fact=0.99, update_rule=\"config\", "
   "gamma=1, implementation=\"orig\", lambda=1)\n\n"
   "Finds the community structure of the graph according to the spinglass\n"
   "community detection method of Reichardt & Bornholdt.\n\n"
   "@param weights: edge weights to be used. Can be a sequence or iterable or\n"
   "  even an edge attribute name.\n"
   "@param spins: integer, the number of spins to use. This is the upper limit\n"
   "  for the number of communities. It is not a problem to supply a\n"
   "  (reasonably) big number here, in which case some spin states will be\n"
   "  unpopulated.\n"
   "@param parupdate: whether to update the spins of the vertices in parallel\n"
   "  (synchronously) or not\n"
   "@param start_temp: the starting temperature\n"
   "@param stop_temp: the stop temperature\n"
   "@param cool_fact: cooling factor for the simulated annealing\n"
   "@param update_rule: specifies the null model of the simulation. Possible\n"
   "  values are C{\"config\"} (a random graph with the same vertex degrees\n"
   "  as the input graph) or C{\"simple\"} (a random graph with the same number\n"
   "  of edges)\n"
   "@param gamma: the gamma argument of the algorithm, specifying the balance\n"
   "  between the importance of present and missing edges within a community.\n"
   "  The default value of 1.0 assigns equal importance to both of them.\n"
   "@param implementation: currently igraph contains two implementations for\n"
   "  the spinglass community detection algorithm. The faster original\n"
   "  implementation is the default. The other implementation is able to take\n"
   "  into account negative weights, this can be chosen by setting\n"
   "  C{implementation} to C{\"neg\"}.\n"
   "@param lambda: the lambda argument of the algorithm, which specifies the\n"
   "  balance between the importance of present and missing negatively\n"
   "  weighted edges within a community. Smaller values of lambda lead\n"
   "  to communities with less negative intra-connectivity. If the argument\n"
   "  is zero, the algorithm reduces to a graph coloring algorithm, using\n"
   "  the number of spins as colors. This argument is ignored if the\n"
   "  original implementation is used.\n"
   "@return: the community membership vector.\n"
  },
  {"community_walktrap",
   (PyCFunction) igraphmodule_Graph_community_walktrap,
   METH_VARARGS | METH_KEYWORDS,
   "community_walktrap(weights=None, steps=None)\n\n"
   "Finds the community structure of the graph according to the random walk\n"
   "method of Latapy & Pons.\n\n"
   "The basic idea of the algorithm is that short random walks tend to stay\n"
   "in the same community. The method provides a dendrogram.\n\n"
   "@attention: this function is wrapped in a more convenient syntax in the\n"
   "  derived class L{Graph}. It is advised to use that instead of this version.\n\n"
   "@param weights: name of an edge attribute or a list containing\n"
   "  edge weights\n"
   "@return: a tuple with the list of merges and the modularity scores corresponding\n"
   "  to each merge\n"
   "\n"
   "@newfield ref: Reference\n"
   "@ref: Pascal Pons, Matthieu Latapy: Computing communities in large networks\n"
   "  using random walks, U{http://arxiv.org/abs/physics/0512106}.\n"
   "@see: modularity()\n"
  },

  /*************/
  /* MATCHINGS */
  /*************/
  {"_is_matching", (PyCFunction)igraphmodule_Graph_is_matching,
   METH_VARARGS | METH_KEYWORDS,
   "_is_matching(matching, types=None)\n\n"
   "Internal function, undocumented.\n\n"
  },
  {"_is_maximal_matching", (PyCFunction)igraphmodule_Graph_is_maximal_matching,
   METH_VARARGS | METH_KEYWORDS,
   "_is_maximal_matching(matching, types=None)\n\n"
   "Internal function, undocumented.\n\n"
   "Use L{Matching.is_maximal} instead.\n"
  },
  {"_maximum_bipartite_matching", (PyCFunction)igraphmodule_Graph_maximum_bipartite_matching,
   METH_VARARGS | METH_KEYWORDS,
   "_maximum_bipartite_matching(types, weights=None)\n\n"
   "Internal function, undocumented.\n\n"
   "@see: L{Graph.maximum_bipartite_matching}\n"
  },

  /**********************/
  /* INTERNAL FUNCTIONS */
  /**********************/
#ifdef IGRAPH_PYTHON3
  {"__graph_as_capsule",
   (PyCFunction) igraphmodule_Graph___graph_as_capsule__,
   METH_VARARGS | METH_KEYWORDS,
   "__graph_as_capsule()\n\n"
   "Returns the igraph graph encapsulated by the Python object as\n"
   "a PyCapsule\n\n."
   "A PyCapsule is practically a regular C pointer, wrapped in a\n"
   "Python object. This function should not be used directly by igraph\n"
   "users, it is useful only in the case when the underlying igraph object\n"
   "must be passed to other C code through Python.\n\n"},
#else
  {"__graph_as_cobject",
   (PyCFunction) igraphmodule_Graph___graph_as_cobject__,
   METH_VARARGS | METH_KEYWORDS,
   "__graph_as_cobject()\n\n"
   "Returns the igraph graph encapsulated by the Python object as\n"
   "a PyCObject\n\n."
   "A PyCObject is practically a regular C pointer, wrapped in a\n"
   "Python object. This function should not be used directly by igraph\n"
   "users, it is useful only in the case when the underlying igraph object\n"
   "must be passed to other C code through Python.\n\n"},
#endif

  {"__register_destructor",
   (PyCFunction) igraphmodule_Graph___register_destructor__,
   METH_VARARGS | METH_KEYWORDS,
   "__register_destructor(destructor)\n\n"
   "Registers a destructor to be called when the object is freed by\n"
   "Python. This function should not be used directly by igraph users."},
  {NULL}
};

/** \ingroup python_interface_graph
 * This structure is the collection of functions necessary to implement
 * the graph as a mapping (i.e. to allow the retrieval and setting of
 * igraph attributes in Python as if it were of a Python mapping type)
 */
PyMappingMethods igraphmodule_Graph_as_mapping = {
  /* __len__ function intentionally left unimplemented */
  0,
  /* returns an attribute by name or returns part of the adjacency matrix */
  (binaryfunc) igraphmodule_Graph_mp_subscript,
  /* sets an attribute by name or sets part of the adjacency matrix */
  (objobjargproc) igraphmodule_Graph_mp_assign_subscript
};

/** \ingroup python_interface
 * \brief Collection of methods to allow numeric operators to be used on the graph
 */
PyNumberMethods igraphmodule_Graph_as_number = {
  0,                            /* nb_add */
  0,                            /*nb_subtract */
  0,                            /*nb_multiply */
#ifndef IGRAPH_PYTHON3
  0,                            /*nb_divide */
#endif
  0,                            /*nb_remainder */
  0,                            /*nb_divmod */
  0,                            /*nb_power */
  0,                            /*nb_negative */
  0,                            /*nb_positive */
  0,                            /*nb_absolute */
  0,                            /*nb_nonzero (2.x) / nb_bool (3.x) */
  (unaryfunc) igraphmodule_Graph_complementer_op, /*nb_invert */
  0,                            /*nb_lshift */
  0,                            /*nb_rshift */
  (binaryfunc) igraphmodule_Graph_intersection, /*nb_and */
  0,                            /*nb_xor */
  (binaryfunc) igraphmodule_Graph_union,  /*nb_or */
#ifndef IGRAPH_PYTHON3
  0,                            /*nb_coerce */
#endif
  0,                            /*nb_int */
  0,                            /*nb_long (2.x) / nb_reserved (3.x)*/
  0,                            /*nb_float */
#ifndef IGRAPH_PYTHON3
  0,                            /*nb_oct */
  0,                            /*nb_hex */
#endif
  0,                            /*nb_inplace_add */
  0,                            /*nb_inplace_subtract */
  0,                            /*nb_inplace_multiply */
#ifndef IGRAPH_PYTHON3
  0,                            /*nb_inplace_divide */
#endif
  0,                            /*nb_inplace_remainder */
  0,                            /*nb_inplace_power */
  0,                            /*nb_inplace_lshift */
  0,                            /*nb_inplace_rshift */
  0,                            /*nb_inplace_and */
  0,                            /*nb_inplace_xor */
  0,                            /*nb_inplace_or */

#ifdef IGRAPH_PYTHON3
  0,                            /*nb_floor_divide */
  0,                            /*nb_true_divide */
  0,                            /*nb_inplace_floor_divide */
  0,                            /*nb_inplace_true_divide */
  0,                            /*nb_index */
#endif
};

/** \ingroup python_interface_graph
 * Python type object referencing the methods Python calls when it performs various operations on an igraph (creating, printing and so on)
 */
PyTypeObject igraphmodule_GraphType = {
  PyVarObject_HEAD_INIT(0, 0)
  "igraph.Graph",               /* tp_name */
  sizeof(igraphmodule_GraphObject), /* tp_basicsize */
  0,                            /* tp_itemsize */
  (destructor) igraphmodule_Graph_dealloc,  /* tp_dealloc */
  0,                            /* tp_print */
  0,                            /* tp_getattr */
  0,                            /* tp_setattr */
  0,                            /* tp_compare (2.x) / tp_reserved (3.x) */
  0,                            /* tp_repr */
  &igraphmodule_Graph_as_number,  /* tp_as_number */
  0,                            /* tp_as_sequence */
  &igraphmodule_Graph_as_mapping, /* tp_as_mapping */
  0,                            /* tp_hash */
  0,                            /* tp_call */
  (reprfunc) igraphmodule_Graph_str,  /* tp_str */
  0,                            /* tp_getattro */
  0,                            /* tp_setattro */
  0,                            /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,  /* tp_flags */
  "Low-level representation of a graph.\n\n"
  "Don't use it directly, use L{igraph.Graph} instead.\n\n"
  "@undocumented: _Bipartite, _Full_Bipartite, _GRG, _Incidence, _is_matching,\n"
  "  _is_maximal_matching, _layout_sugiyama, _maximum_bipartite_matching,\n"
  "  _spanning_tree\n"
  "@deffield ref: Reference",  /* tp_doc */
  (traverseproc) igraphmodule_Graph_traverse, /* tp_traverse */
  (inquiry) igraphmodule_Graph_clear, /* tp_clear */
  0,                            /* tp_richcompare */
  offsetof(igraphmodule_GraphObject, weakreflist),  /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  igraphmodule_Graph_methods,   /* tp_methods */
  0,                            /* tp_members */
  0,                            /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc) igraphmodule_Graph_init, /* tp_init */
  0,                            /* tp_alloc */
  igraphmodule_Graph_new,       /* tp_new */
  0,                            /* tp_free */
};

#undef CREATE_GRAPH

