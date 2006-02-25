#include <Python.h>
#include "igraph.h"
#include "common.h"
#include "error.h"
#include "convert.h"
#include "graphobject.h"
#include "vertexseqobject.h"
#include "vertexobject.h"
#include "edgeseqobject.h"
#include "edgeobject.h"

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

int igraphmodule_igraph_progress_hook(const char* message, real_t percent,
				       void* data) {
  /* Look up _progress_handler in module namespace */
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
  const char* kwlist[] = {"vs", "coords", NULL};
  PyObject *vs, *o, *o1, *o2, *coords = Py_False;
  igraph_matrix_t mtrx;
  igraph_vector_t result;
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
    MATRIX(mtrx, i, 0)=(real_t)PyFloat_AsDouble(o1);
    MATRIX(mtrx, i, 1)=(real_t)PyFloat_AsDouble(o2);
    Py_DECREF(o1);
    Py_DECREF(o2);
  }

  if (igraph_vector_init(&result, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&mtrx);
    return NULL;
  }
  
  if (igraph_convex_hull(&mtrx, &result, PyObject_IsTrue(coords))) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&mtrx);
    igraph_vector_destroy(&result);
    return NULL;
  }
  
  if (PyObject_IsTrue(coords))
    o=igraphmodule_vector_t_to_PyList_pairs(&result);
  else
    o=igraphmodule_vector_t_to_PyList(&result);
  
  igraph_matrix_destroy(&mtrx);
  igraph_vector_destroy(&result);

  return o;
}

/** \ingroup python_interface
 * \brief Method table for the igraph Python module
 */
static PyMethodDef igraphmodule_methods[] = 
{
  {"convex_hull", igraphmodule_convex_hull, METH_VARARGS,
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
  
  igraphmodule_GraphType.tp_new = igraphmodule_Graph_new;
  igraphmodule_GraphType.tp_init = (initproc)igraphmodule_Graph_init;
  igraphmodule_GraphType.tp_methods = igraphmodule_Graph_methods;
  igraphmodule_GraphType.tp_getset = igraphmodule_Graph_getseters;
  
  if (PyType_Ready(&igraphmodule_GraphType) < 0) return;
  
  igraphmodule_InternalError =
    PyErr_NewException("igraph.InternalError", NULL, NULL);
  
  Py_INCREF(igraphmodule_InternalError);
  
  m = Py_InitModule3("igraph", igraphmodule_methods,
		     "Python interface for the igraph library");
  
  Py_INCREF(&igraphmodule_GraphType);
  // Maybe the next line is unnecessary?
  // Py_INCREF(&igraphmodule_VertexSeqType);
  
  PyModule_AddObject(m, "Graph", (PyObject*)&igraphmodule_GraphType);
  // Maybe the next line is unnecessary?
  // PyModule_AddObject(m, "VertexSeq", (PyObject*)&igraphmodule_VertexSeqType);
  
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
  
  /* initialize error and progress handler */
  igraph_set_error_handler(igraphmodule_igraph_error_hook);
  igraph_set_progress_handler(igraphmodule_igraph_progress_hook);
}
