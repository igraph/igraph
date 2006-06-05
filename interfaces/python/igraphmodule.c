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

#include <Python.h>
#include <pythonrun.h>
#include "igraph.h"
#include "attributes.h"
#include "common.h"
#include "error.h"
#include "convert.h"
#include "graphobject.h"
#include "vertexseqobject.h"
#include "vertexobject.h"
#include "edgeseqobject.h"
#include "edgeobject.h"
#include "bfsiter.h"

/**
 * \defgroup python_interface Python module implementation
 * \brief Functions implementing a Python interface to \a igraph
 * 
 * These functions provide a way to access \a igraph functions from Python.
 * It should be of interest of \a igraph developers only. Classes, functions
 * and methods exposed to Python are still to be documented. Until it is done,
 * just type the following to get help about \a igraph functions in Python
 * (assuming you have \c igraph.so somewhere in your Python library path):
 * 
 * \verbatim
import igraph
help(igraph)
help(igraph.Graph)
\endverbatim
 * 
 * Most of the functions provided here share the same calling conventions
 * (which are determined by the Python/C API). Since the role of the
 * arguments are the same across many functions, I won't explain them
 * everywhere, just give a quick overview of the common argument names here.
 * 
 * \param self the Python igraph.Graph object the method is working on
 * \param args pointer to the Python tuple containing the arguments
 * \param kwds pointer to the Python hash containing the keyword parameters
 * \param type the type object of a Python igraph.Graph object. Used usually
 * in constructors and class methods.
 * 
 * Any arguments not documented here should be mentioned at the documentation
 * of the appropriate method.
 * 
 * The functions which implement a Python method always return a pointer to
 * a \c PyObject. According to Python conventions, this is \c NULL if and
 * only if an exception was thrown by the method (or any of the functions
 * it has called). When I explain the return value of a function which
 * provides interface to an \a igraph function, I won't cover the case of
 * returning a \c NULL value, because this is the same for every such method.
 * The conclusion is that a method can return \c NULL even if I don't state
 * it explicitly.
 * 
 * Also please take into consideration that I'm documenting the C calls
 * with the abovementioned parameters here, and \em not the Python methods
 * which are presented to the user using the Python interface of \a igraph.
 * If you are looking for the documentation of the classes, methods and
 * functions exposed to Python, please use the \c help calls from Python
 * or use \c pydoc to generate a formatted version.
 *
 * \section weakrefs The usage of weak references in the Python interface
 * 
 * Many classes implemented in the Python interface (e.g. VertexSeq, Vertex...)
 * use weak references to keep track of the graph they are referencing to.
 * The use of weak references is twofold:
 * 
 * -# If we assign a VertexSeq or a Vertex of a given graph to a local
 *    variable and then destroy the graph, real references keep the graph
 *    alive and do not return the memory back to Python.
 * -# If we use real references, a Graph object will hold a reference
 *    to its VertexSeq (because we don't want to allocate a new VertexSeq
 *    object for the same graph every time it is requested), and the
 *    VertexSeq will also hold a reference to the Graph. This is a circular
 *    reference. Python does not reclaim the memory occupied by the Graph
 *    back when the Graph is destroyed, because the VertexSeq is holding a
 *    reference to it. Similarly, VertexSeq doesn't get freed because the
 *    Graph is holding a reference to it. These situations can only be
 *    resolved by the Python garbage collector which is invoked at regular
 *    intervals. Unfortunately, the garbage collector refuses to break
 *    circular references and free the objects participating in the circle
 *    when any of the objects has a \c __del__ method. In this case,
 *    \c igraph.Graph has one (which frees the underlying \c igraph_t
 *    graph), therefore our graphs never get freed when we use real
 *    references.
 */

static PyObject* igraphmodule_progress_handler=NULL;

static int igraphmodule_igraph_interrupt_hook(void* data) {
  if (PyErr_CheckSignals()) {
    return IGRAPH_FAILURE;
  }
  return IGRAPH_SUCCESS;
}

int igraphmodule_igraph_progress_hook(const char* message, igraph_real_t percent,
				       void* data) {
  if (igraphmodule_progress_handler) {
    PyObject *result;
    if (PyCallable_Check(igraphmodule_progress_handler)) {
      result=PyObject_CallFunction(igraphmodule_progress_handler,
				   "sd", message, (double)percent);
      Py_DECREF(result);
    }
  }
  
  return 0;
}


PyObject* igraphmodule_set_progress_handler(PyObject* self, PyObject* args) {
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  if (!PyCallable_Check(o)) {
    PyErr_SetString(PyExc_TypeError, "Progress handler must be callable.");
    return NULL;
  }
  igraphmodule_progress_handler=o;
  Py_RETURN_NONE;
}


PyObject* igraphmodule_convex_hull(PyObject* self, PyObject* args, PyObject* kwds) {
  char* kwlist[] = {"vs", "coords", NULL};
  PyObject *vs, *o, *o1=0, *o2=0, *coords = Py_False;
  igraph_matrix_t mtrx;
  igraph_vector_t result;
  igraph_matrix_t resmat;
  long no_of_nodes, i;
  
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O", kwlist, &PyList_Type, &vs, &coords))
    return NULL;
  
  no_of_nodes=PyList_Size(vs);
  if (igraph_matrix_init(&mtrx, no_of_nodes, 2)) {
    igraphmodule_handle_igraph_error();
    return NULL;
  }
  for (i=0; i<no_of_nodes; i++) {
    o=PyList_GetItem(vs, i);
    if (PyList_Check(o)) {
      if (PyList_Size(o) >= 2) {
	o1=PyList_GetItem(o, 0);
	o2=PyList_GetItem(o, 1);
	if (PyList_Size(o) > 2)
	  PyErr_Warn(PyExc_Warning, "vertex with more than 2 coordinates found, considering only the first 2");
      } else {
	PyErr_SetString(PyExc_TypeError, "vertex with less than 2 coordinates found");
	igraph_matrix_destroy(&mtrx);
	return NULL;
      }
    } else if (PyTuple_Check(o)) {
      if (PyTuple_Size(o) >= 2) {
	o1=PyTuple_GetItem(o, 0);
	o2=PyTuple_GetItem(o, 1);
	if (PyTuple_Size(o) > 2)
	  PyErr_Warn(PyExc_Warning, "vertex with more than 2 coordinates found, considering only the first 2");
      } else {
	PyErr_SetString(PyExc_TypeError, "vertex with less than 2 coordinates found");
	igraph_matrix_destroy(&mtrx);
	return NULL;
      }
    }
    
    if (!PyNumber_Check(o1) || !PyNumber_Check(o2)) {
      PyErr_SetString(PyExc_TypeError, "vertex coordinates must be numeric");
      igraph_matrix_destroy(&mtrx);
      return NULL;
    }
    /* o, o1 and o2 were borrowed, but now o1 and o2 are actual references! */
    o1=PyNumber_Float(o1); o2=PyNumber_Float(o2);
    if (!o1 || !o2) {
      PyErr_SetString(PyExc_TypeError, "vertex coordinate conversion to float failed");
      Py_XDECREF(o1);
      Py_XDECREF(o2);
      igraph_matrix_destroy(&mtrx);
      return NULL;
    }
    MATRIX(mtrx, i, 0)=(igraph_real_t)PyFloat_AsDouble(o1);
    MATRIX(mtrx, i, 1)=(igraph_real_t)PyFloat_AsDouble(o2);
    Py_DECREF(o1);
    Py_DECREF(o2);
  }

  if (!PyObject_IsTrue(coords)) {
    if (igraph_vector_init(&result, 0)) {
      igraphmodule_handle_igraph_error();
      igraph_matrix_destroy(&mtrx);
      return NULL;
    }
    if (igraph_convex_hull(&mtrx, &result, 0)) {
      igraphmodule_handle_igraph_error();
      igraph_matrix_destroy(&mtrx);
      igraph_vector_destroy(&result);
      return NULL;
    }    
    o=igraphmodule_vector_t_to_PyList(&result);
    igraph_vector_destroy(&result);
  } else {
    if (igraph_matrix_init(&resmat, 0, 0)) {
      igraphmodule_handle_igraph_error();
      igraph_matrix_destroy(&mtrx);
      return NULL;
    }
    if (igraph_convex_hull(&mtrx, 0, &resmat)) {
      igraphmodule_handle_igraph_error();
      igraph_matrix_destroy(&mtrx);
      igraph_matrix_destroy(&resmat);
      return NULL;
    }        
    o=igraphmodule_matrix_t_to_PyList(&resmat, IGRAPHMODULE_TYPE_FLOAT);
    igraph_matrix_destroy(&resmat);
  }
  
  igraph_matrix_destroy(&mtrx);

  return o;
}


/* Attribute handlers for the Python interface */

/* Initialization */ 
static int igraphmodule_i_attribute_init(igraph_t *graph) {
  PyObject** attrs;
  int i;
  
  attrs=(PyObject**)calloc(3, sizeof(PyObject*));
  /* printf("Created attribute table at %p\n", attrs); */
  if (!attrs)
    IGRAPH_ERROR("not enough memory to allocate attribute hashes", IGRAPH_ENOMEM);
  
  for (i=0; i<3; i++) {
    attrs[i] = PyDict_New();
  }
  graph->attr=(void*)attrs;
  
  return IGRAPH_SUCCESS;
}

/* Destruction */
static void igraphmodule_i_attribute_destroy(igraph_t *graph) {
  PyObject** attrs;
  int i;
 
  /* printf("Destroying attribute table\n"); */
  if (graph->attr) {
    attrs=(PyObject**)graph->attr;
    for (i=0; i<3; i++) {
      Py_DECREF(attrs[i]);
    }
    free(attrs);
  }
}

/* Copying */
static int igraphmodule_i_attribute_copy(igraph_t *to, const igraph_t *from) {
  PyObject **fromattrs, **toattrs, *key, *value, *newval, *o;
  int i, j, pos;
 
  /* printf("Copying attribute table\n"); */
  if (from->attr) {
    fromattrs=(PyObject**)from->attr;
    /* what to do with the original value of toattrs? */
    toattrs=to->attr=(PyObject**)calloc(3, sizeof(PyObject*));
    for (i=0; i<3; i++) {
      if (!PyDict_Check(fromattrs[i])) {
	toattrs[i]=fromattrs[i];
	Py_XINCREF(o);
	continue;
      }
      
      toattrs[i]=PyDict_New();
      
      pos=0;
      while (PyDict_Next(fromattrs[i], &pos, &key, &value)) {
	/* value is only borrowed, so copy it */
	if (i>0) {
	  newval=PyList_New(PyList_GET_SIZE(value));
	  for (j=0; j<PyList_GET_SIZE(value); j++) {
	    o=PyList_GetItem(value, j);
	    Py_INCREF(o);
	    PyList_SetItem(newval, j, o);
	  }
	} else {
	  newval=value;
	  Py_INCREF(newval);
	}
	PyDict_SetItem(toattrs[i], key, newval);
      }
    }
  }
  return IGRAPH_SUCCESS;
}

/* Adding vertices */
static int igraphmodule_i_attribute_add_vertices(igraph_t *graph, long int nv, igraph_vector_ptr_t *attr) {
  /* Extend the end of every value in the vertex hash with nv pieces of None */
  PyObject *key, *value, *dict, *l, *r;
  int pos=0;
  long int i;
  
  if (!graph->attr) return IGRAPH_SUCCESS;
  if (nv<=0) return IGRAPH_SUCCESS;
  
  dict=((PyObject**)graph->attr)[1];
  if (!PyDict_Check(dict)) 
    IGRAPH_ERROR("vertex attribute hash type mismatch", IGRAPH_EINVAL);
  while (PyDict_Next(dict, &pos, &key, &value)) {
    /* Create an object with nv pieces of None */
    if (!PyList_Check(value))
      IGRAPH_ERROR("vertex attribute hash member is not a list", IGRAPH_EINVAL);

    for (i=0; i<nv; i++) {
      Py_INCREF(Py_None); /* is it correct here? */
      if (PyList_Append(value, Py_None) == -1) {
	Py_DECREF(Py_None);
	IGRAPH_ERROR("can't extend a vertex attribute hash member", IGRAPH_FAILURE);
      }
    }
  }
  
  return IGRAPH_SUCCESS;
}

static void igraphmodule_i_attribute_delete_edges(igraph_t *graph, const igraph_vector_t *idx);

/* Deleting vertices */
static void igraphmodule_i_attribute_delete_vertices(igraph_t *graph,
						     const igraph_vector_t *eidx,
						     const igraph_vector_t *vidx) {
  long int n, i, ndeleted=0;
  PyObject *key, *value, *dict, *o;
  int pos=0;
  
  /* Reindexing vertices */
  dict=((PyObject**)graph->attr)[1];
  if (!PyDict_Check(dict)) return;

  n=igraph_vector_size(vidx);
  for (i=0; i<n; i++) {
    /*printf("%ld:%f ", i, VECTOR(*idx)[i]);*/
    if (!VECTOR(*vidx)[i]) {
      ndeleted++;
      continue;
    }

    pos=0;
    /* TODO: maybe it would be more efficient to get the values from the
     * hash in advance? */
    while (PyDict_Next(dict, &pos, &key, &value)) {
      /* Move the element from index i to VECTOR(*idx)[i]-1 */
      o=PyList_GetItem(value, i);
      if (!o) {
	/* IndexError is already set, clear it and return */
	PyErr_Clear();
	return;
      }
      Py_INCREF(o);
      PyList_SetItem(value, VECTOR(*vidx)[i]-1, o);
    }
  }
  /*printf("\n");*/
  
  /* Clear the remaining parts of the lists that aren't needed anymore */
  pos=0;
  while (PyDict_Next(dict, &pos, &key, &value)) {
    n=PySequence_Size(value);
    if (PySequence_DelSlice(value, n-ndeleted, n) == -1) return;
    /*printf("key: "); PyObject_Print(key, stdout, Py_PRINT_RAW); printf("\n");
    printf("value: "); PyObject_Print(value, stdout, Py_PRINT_RAW); printf("\n");*/
  }
  
  igraphmodule_i_attribute_delete_edges(graph, eidx);

  return;
}

/* Adding edges */
static int igraphmodule_i_attribute_add_edges(igraph_t *graph, const igraph_vector_t *edges, igraph_vector_ptr_t *attr) {
  /* Extend the end of every value in the edge hash with ne pieces of None */
  PyObject *key, *value, *dict, *l, *r;
  int pos=0;
  long int i, ne;

  ne=igraph_vector_size(edges);
  if (!graph->attr) return IGRAPH_SUCCESS;
  if (ne<=0) return IGRAPH_SUCCESS;
  
  dict=((PyObject**)graph->attr)[2];
  if (!PyDict_Check(dict)) 
    IGRAPH_ERROR("edge attribute hash type mismatch", IGRAPH_EINVAL);
  while (PyDict_Next(dict, &pos, &key, &value)) {
    if (!PyList_Check(value))
      IGRAPH_ERROR("edge attribute hash member is not a list", IGRAPH_EINVAL);

    for (i=0; i<ne; i++) {
      Py_INCREF(Py_None); /* is it correct here? */
      if (PyList_Append(value, Py_None) == -1) {
	Py_DECREF(Py_None);
	IGRAPH_ERROR("can't extend a vertex attribute hash member", IGRAPH_FAILURE);
      }
    }
  }
  
  /*pos=0;
  while (PyDict_Next(dict, &pos, &key, &value)) {
    printf("key: "); PyObject_Print(key, stdout, Py_PRINT_RAW); printf("\n");
    printf("value: "); PyObject_Print(value, stdout, Py_PRINT_RAW); printf("\n");
  }*/
  
  return IGRAPH_SUCCESS;
}

/* Deleting edges */
static void igraphmodule_i_attribute_delete_edges(igraph_t *graph, const igraph_vector_t *idx) {
  long int n, i, ndeleted=0;
  PyObject *key, *value, *dict, *o;
  int pos=0;
  
  dict=((PyObject**)graph->attr)[2];
  if (!PyDict_Check(dict)) return;

  n=igraph_vector_size(idx);
  for (i=0; i<n; i++) {
    /*printf("%ld:%f ", i, VECTOR(*idx)[i]);*/
    if (!VECTOR(*idx)[i]) {
      ndeleted++;
      continue;
    }

    pos=0;
    /* TODO: maybe it would be more efficient to get the values from the
     * hash in advance? */
    while (PyDict_Next(dict, &pos, &key, &value)) {
      /* Move the element from index i to VECTOR(*idx)[i]-1 */
      o=PyList_GetItem(value, i);
      if (!o) {
	/* IndexError is already set, clear it and return */
	PyErr_Clear();
	return;
      }
      Py_INCREF(o);
      PyList_SetItem(value, VECTOR(*idx)[i]-1, o);
    }
  }
  /*printf("\n");*/
  
  /* Clear the remaining parts of the lists that aren't needed anymore */
  pos=0;
  while (PyDict_Next(dict, &pos, &key, &value)) {
    n=PySequence_Size(value);
    if (PySequence_DelSlice(value, n-ndeleted, n) == -1) return;
    /*printf("key: "); PyObject_Print(key, stdout, Py_PRINT_RAW); printf("\n");
    printf("value: "); PyObject_Print(value, stdout, Py_PRINT_RAW); printf("\n");*/
  }
  
  return;
}

/* Getting attribute names and types */
static int igraphmodule_i_attribute_get_info(const igraph_t *graph,
					     igraph_strvector_t *gnames,
					     igraph_vector_t *gtypes,
					     igraph_strvector_t *vnames,
					     igraph_vector_t *vtypes,
					     igraph_strvector_t *enames,
					     igraph_vector_t *etypes) {
  igraph_strvector_t *names[3] = { gnames, vnames, enames };
  igraph_vector_t *types[3] = { gtypes, vtypes, etypes };
  long int i, j, k;
  
  for (i=0; i<3; i++) {
    igraph_strvector_t *n = names[i];
    igraph_vector_t *t = types[i];
    PyObject *dict = ((PyObject**)graph->attr)[i];
    PyObject *keys;
    keys=PyDict_Keys(dict);
    if (!keys) IGRAPH_ERROR("Internal error in PyDict_Keys", IGRAPH_FAILURE);
    
    if (n) {
      j=igraphmodule_PyList_to_strvector_t(keys, n);
      if (j) return j;
    }
    if (t) {
      k=PyList_Size(keys);
      igraph_vector_init(t, k);
      for (j=0; j<k; j++) VECTOR(*t)[j]=IGRAPH_ATTRIBUTE_PY_OBJECT;
    }
    
    Py_DECREF(keys);
  }
  
  return 0;
}

/* Checks whether the graph has a graph/vertex/edge attribute with the given name */
igraph_bool_t igraphmodule_i_attribute_has_attr(const igraph_t *graph,
						igraph_attribute_elemtype_t type,
						const char* name) {
  long int attrnum;
  PyObject *o, *dict;
  switch (type) {
  case IGRAPH_ATTRIBUTE_GRAPH: attrnum=0; break;
  case IGRAPH_ATTRIBUTE_VERTEX: attrnum=1; break;
  case IGRAPH_ATTRIBUTE_EDGE: attrnum=2; break;
  default: return 0; break;
  }
  dict = ((PyObject**)graph->attr)[attrnum];
  o = PyDict_GetItemString(dict, name);
  return o != 0;
}

/* Returns the type of a given attribute */
int igraphmodule_i_attribute_get_type(const igraph_t *graph,
				      igraph_attribute_type_t *type,
				      igraph_attribute_elemtype_t elemtype,
				      const char *name) {
  long int attrnum;
  PyObject *o, *dict;
  switch (elemtype) {
  case IGRAPH_ATTRIBUTE_GRAPH: attrnum=0; break;
  case IGRAPH_ATTRIBUTE_VERTEX: attrnum=1; break;
  case IGRAPH_ATTRIBUTE_EDGE: attrnum=2; break;
  default: IGRAPH_ERROR("No such attribute type", IGRAPH_EINVAL); break;
  }
  dict = ((PyObject**)graph->attr)[attrnum];
  o = PyDict_GetItemString(dict, name);
  if (o != 0) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  *type = IGRAPH_ATTRIBUTE_PY_OBJECT;
  return 0;
}

static igraph_attribute_table_t igraphmodule_i_attribute_table = {
  igraphmodule_i_attribute_init,
  igraphmodule_i_attribute_destroy,
  igraphmodule_i_attribute_copy,
  igraphmodule_i_attribute_add_vertices,
  igraphmodule_i_attribute_delete_vertices,
  igraphmodule_i_attribute_add_edges,
  igraphmodule_i_attribute_delete_edges,
  igraphmodule_i_attribute_get_info,
  igraphmodule_i_attribute_has_attr,
  igraphmodule_i_attribute_get_type,
  0,
  0,
  0,
  0,
  0,
  0,
};

/** \ingroup python_interface
 * \brief Method table for the igraph Python module
 */
static PyMethodDef igraphmodule_methods[] = 
{
  {"convex_hull", (PyCFunction)igraphmodule_convex_hull, METH_VARARGS,
      "Calculates the convex hull of a given point set\n\n"
      "Keyword arguments:\n"
      "vs     -- the point set as a list of lists\n"
      "coords -- if True, the function returns the coordinates of the\n"
      "          corners of the convex hull polygon, otherwise returns\n"
      "          the corner indices. Optional, defaults to False."
  },
  {"set_progress_handler", igraphmodule_set_progress_handler, METH_VARARGS,
      "Sets the handler to be called when igraph is performing a long operation."
  },
  {NULL, NULL, 0, NULL}
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

extern PyObject* igraphmodule_InternalError;

PyMODINIT_FUNC
initigraph(void)
{
  PyObject* m;  ///< igraph module object
  
  igraphmodule_VertexSeqType.tp_traverse = (traverseproc)igraphmodule_VertexSeq_traverse;
  igraphmodule_VertexSeqType.tp_clear = (inquiry)igraphmodule_VertexSeq_clear;
  if (PyType_Ready(&igraphmodule_VertexSeqType) < 0) return;
  
  igraphmodule_VertexType.tp_traverse = (traverseproc)igraphmodule_Vertex_traverse;
  igraphmodule_VertexType.tp_clear = (inquiry)igraphmodule_Vertex_clear;
  if (PyType_Ready(&igraphmodule_VertexType) < 0) return;
  
  igraphmodule_EdgeSeqType.tp_traverse = (traverseproc)igraphmodule_EdgeSeq_traverse;
  igraphmodule_EdgeSeqType.tp_clear = (inquiry)igraphmodule_EdgeSeq_clear;
  if (PyType_Ready(&igraphmodule_EdgeSeqType) < 0) return;
  
  igraphmodule_EdgeType.tp_traverse = (traverseproc)igraphmodule_Edge_traverse;
  igraphmodule_EdgeType.tp_clear = (inquiry)igraphmodule_Edge_clear;
  if (PyType_Ready(&igraphmodule_EdgeType) < 0) return;
  
  igraphmodule_GraphType.tp_methods = igraphmodule_Graph_methods;
  igraphmodule_GraphType.tp_getset = igraphmodule_Graph_getseters;
  if (PyType_Ready(&igraphmodule_GraphType) < 0) return;
  
  if (PyType_Ready(&igraphmodule_BFSIterType) < 0) return;
  
  igraphmodule_InternalError =
    PyErr_NewException("igraph.InternalError", NULL, NULL);
  
  Py_INCREF(igraphmodule_InternalError);
  
  m = Py_InitModule3("igraph", igraphmodule_methods,
		     "Python interface for the igraph library");
  
  Py_INCREF(&igraphmodule_GraphType);
  
  PyModule_AddObject(m, "Graph", (PyObject*)&igraphmodule_GraphType);
  PyModule_AddObject(m, "BFSIter", (PyObject*)&igraphmodule_BFSIterType);
  
  PyModule_AddObject(m, "InternalError", igraphmodule_InternalError);
  PyModule_AddIntConstant(m, "OUT", IGRAPH_OUT);
  PyModule_AddIntConstant(m, "IN", IGRAPH_IN);
  PyModule_AddIntConstant(m, "ALL", IGRAPH_ALL);
  PyModule_AddIntConstant(m, "STAR_OUT", IGRAPH_STAR_OUT);
  PyModule_AddIntConstant(m, "STAR_IN", IGRAPH_STAR_IN);
  PyModule_AddIntConstant(m, "STAR_UNDIRECTED", IGRAPH_STAR_UNDIRECTED);
  PyModule_AddIntConstant(m, "TREE_OUT", IGRAPH_TREE_OUT);
  PyModule_AddIntConstant(m, "TREE_IN", IGRAPH_TREE_IN);
  PyModule_AddIntConstant(m, "TREE_UNDIRECTED", IGRAPH_TREE_UNDIRECTED);
  PyModule_AddIntConstant(m, "STRONG", IGRAPH_STRONG);
  PyModule_AddIntConstant(m, "WEAK", IGRAPH_WEAK);
  PyModule_AddIntConstant(m, "GET_ADJACENCY_UPPER", IGRAPH_GET_ADJACENCY_UPPER);
  PyModule_AddIntConstant(m, "GET_ADJACENCY_LOWER", IGRAPH_GET_ADJACENCY_LOWER);
  PyModule_AddIntConstant(m, "GET_ADJACENCY_BOTH", IGRAPH_GET_ADJACENCY_BOTH);
  PyModule_AddIntConstant(m, "REWIRING_SIMPLE", IGRAPH_REWIRING_SIMPLE);
  PyModule_AddIntConstant(m, "ADJ_DIRECTED", IGRAPH_ADJ_DIRECTED);
  PyModule_AddIntConstant(m, "ADJ_UNDIRECTED", IGRAPH_ADJ_UNDIRECTED);
  PyModule_AddIntConstant(m, "ADJ_MAX", IGRAPH_ADJ_MAX);
  PyModule_AddIntConstant(m, "ADJ_MIN", IGRAPH_ADJ_MIN);
  PyModule_AddIntConstant(m, "ADJ_PLUS", IGRAPH_ADJ_PLUS);
  PyModule_AddIntConstant(m, "ADJ_UPPER", IGRAPH_ADJ_UPPER);
  PyModule_AddIntConstant(m, "ADJ_LOWER", IGRAPH_ADJ_LOWER);
  
  /* initialize error, progress and interruption handler */
  igraph_set_error_handler(igraphmodule_igraph_error_hook);
  igraph_set_progress_handler(igraphmodule_igraph_progress_hook);
  igraph_set_interruption_handler(igraphmodule_igraph_interrupt_hook);
  
  /* initialize attribute handlers */
  igraph_i_set_attribute_table(&igraphmodule_i_attribute_table);
}
