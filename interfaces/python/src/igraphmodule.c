/* vim:set ts=2 sw=2 sts=2 et: */
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
#include <igraph/igraph.h>
#include "common.h"
#include "error.h"
#include "convert.h"
#include "graphobject.h"
#include "vertexseqobject.h"
#include "vertexobject.h"
#include "edgeseqobject.h"
#include "edgeobject.h"
#include "arpackobject.h"
#include "bfsiter.h"
//#include "config.h"

extern double igraph_i_fdiv(double, double);

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
    IGRAPH_FINALLY_FREE();
    return IGRAPH_INTERRUPTED;
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
      if (result) Py_DECREF(result);
      else return IGRAPH_INTERRUPTED;
    }
  }
  
  return IGRAPH_SUCCESS;
}

PyObject* igraphmodule_set_progress_handler(PyObject* self, PyObject* args) {
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  if (!PyCallable_Check(o) && o != Py_None) {
    PyErr_SetString(PyExc_TypeError, "Progress handler must be callable.");
    return NULL;
  }
  Py_XDECREF(igraphmodule_progress_handler);
  if (o == Py_None)
     igraphmodule_progress_handler=NULL;
  else {
    Py_INCREF(o);
    igraphmodule_progress_handler=o;
  }
  Py_RETURN_NONE;
}


PyObject* igraphmodule_convex_hull(PyObject* self, PyObject* args, PyObject* kwds) {
  static char* kwlist[] = {"vs", "coords", NULL};
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
    o=igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
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


PyObject* igraphmodule_community_to_membership(PyObject *self,
  PyObject *args, PyObject *kwds) {
  static char* kwlist[] = { "merges", "nodes", "steps", "return_csize", NULL };
  PyObject *merges_o, *return_csize = Py_False, *result_o;
  igraph_matrix_t merges;
  igraph_vector_t result, csize, *csize_p = 0;
  long int nodes, steps;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!ll|O", kwlist,
      &PyList_Type, &merges_o, &nodes, &steps, &return_csize)) return NULL;

  if (igraphmodule_PyList_to_matrix_t(merges_o, &merges)) return NULL;

  if (igraph_vector_init(&result, nodes)) {
    igraphmodule_handle_igraph_error();
    igraph_matrix_destroy(&merges);
    return NULL;
  }

  if (PyObject_IsTrue(return_csize)) {
	igraph_vector_init(&csize, 0);
	csize_p = &csize;
  }

  if (igraph_community_to_membership(&merges, nodes, steps, &result, csize_p)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&result);
    if (csize_p) igraph_vector_destroy(csize_p);
    igraph_matrix_destroy(&merges);
    return NULL;
  }
  igraph_matrix_destroy(&merges);

  result_o = igraphmodule_vector_t_to_PyList(&result, IGRAPHMODULE_TYPE_INT);
  igraph_vector_destroy(&result);

  if (csize_p) {
	PyObject* csize_o = igraphmodule_vector_t_to_PyList(csize_p, IGRAPHMODULE_TYPE_INT);
	igraph_vector_destroy(csize_p);
	if (csize_o) return Py_BuildValue("NN", result_o, csize_o);
	Py_DECREF(result_o);
	return NULL;
  }

  return result_o;
}

/* Attribute handlers for the Python interface */

/* Initialization */ 
static int igraphmodule_i_attribute_init(igraph_t *graph, igraph_vector_ptr_t *attr) {
  PyObject** attrs;
  long int i, n;
  
  attrs=(PyObject**)calloc(3, sizeof(PyObject*));
  if (!attrs)
    IGRAPH_ERROR("not enough memory to allocate attribute hashes", IGRAPH_ENOMEM);
  
  for (i=0; i<3; i++) {
    attrs[i] = PyDict_New();
    RC_ALLOC("dict", attrs[i]);
  }
  graph->attr=(void*)attrs;

  /* See if we have graph attributes */
  if (attr) {
    PyObject *dict=attrs[ATTRHASH_IDX_GRAPH], *value;
    char *s;
    n = igraph_vector_ptr_size(attr);
    for (i=0; i<n; i++) {
      igraph_i_attribute_record_t *attr_rec;
      attr_rec = VECTOR(*attr)[i];
      switch (attr_rec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
        value=PyFloat_FromDouble((double)VECTOR(*(igraph_vector_t*)attr_rec->value)[0]);
        break;
      case IGRAPH_ATTRIBUTE_STRING:
        igraph_strvector_get((igraph_strvector_t*)attr_rec->value, 0, &s);
        if (s == 0) value=PyString_FromString("");
        else value=PyString_FromString(s);
        break;
      default:
        IGRAPH_WARNING("unsupported attribute type (not string and not numeric)");
        value=0;
        break;
      }
      if (value) {
        if (PyDict_SetItemString(dict, attr_rec->name, value)) {
          Py_DECREF(value);
          Py_DECREF(attrs[0]);
          Py_DECREF(attrs[1]);
          Py_DECREF(attrs[2]);
          free(graph->attr); graph->attr = 0;
          IGRAPH_ERROR("failed to add attributes to graph attribute hash",
                       IGRAPH_FAILURE);
        }
        Py_DECREF(value);
        value=0;
      }
    }
  }

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
static int igraphmodule_i_attribute_copy(igraph_t *to, const igraph_t *from,
  igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea) {
  PyObject **fromattrs, **toattrs, *key, *value, *newval, *o=NULL;
  igraph_bool_t copy_attrs[3] = { ga, va, ea };
  int i, j;
  Py_ssize_t pos = 0;
 
  if (from->attr) {
    fromattrs=(PyObject**)from->attr;
    /* what to do with the original value of toattrs? */
    toattrs=to->attr=(PyObject**)calloc(3, sizeof(PyObject*));
    for (i=0; i<3; i++) {
      if (!copy_attrs[i]) {
        toattrs[i] = PyDict_New();
        RC_ALLOC("dict (copying, empty)", toattrs[i]);
        continue;
      }

      if (!PyDict_Check(fromattrs[i])) {
        toattrs[i]=fromattrs[i];
        Py_XINCREF(o);
        continue;
      }
      
      toattrs[i]=PyDict_New();
      RC_ALLOC("dict (copying)", toattrs[i]);
      
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
        Py_DECREF(newval); /* compensate for PyDict_SetItem */
      }
    }
  }
  return IGRAPH_SUCCESS;
}

/* Adding vertices */
static int igraphmodule_i_attribute_add_vertices(igraph_t *graph, long int nv, igraph_vector_ptr_t *attr) {
  /* Extend the end of every value in the vertex hash with nv pieces of None */
  PyObject *key, *value, *dict;
  long int i, j, k, l;
  igraph_i_attribute_record_t *attr_rec;
  igraph_bool_t *added_attrs=0;
  Py_ssize_t pos = 0;

  if (!graph->attr) return IGRAPH_SUCCESS;
  if (nv<=0) return IGRAPH_SUCCESS;

  if (attr) {
    added_attrs = (igraph_bool_t*)calloc((size_t)igraph_vector_ptr_size(attr),
                                         sizeof(igraph_bool_t));
    if (!added_attrs)
      IGRAPH_ERROR("can't add vertex attributes", IGRAPH_ENOMEM);
    IGRAPH_FINALLY(free, added_attrs);
  }

  dict=((PyObject**)graph->attr)[1];
  if (!PyDict_Check(dict)) 
    IGRAPH_ERROR("vertex attribute hash type mismatch", IGRAPH_EINVAL);

  while (PyDict_Next(dict, &pos, &key, &value)) {
    if (!PyString_Check(key))
      IGRAPH_ERROR("vertex attribute hash key is not a string", IGRAPH_EINVAL);
    if (!PyList_Check(value))
      IGRAPH_ERROR("vertex attribute hash member is not a list", IGRAPH_EINVAL);
    /* Check if we have specific values for the given attribute */
    attr_rec=0;
    if (attr) {
      j=igraph_vector_ptr_size(attr);
      for (i=0; i<j; i++) {
        attr_rec=VECTOR(*attr)[i];
        if (!strcmp(attr_rec->name, PyString_AS_STRING(key))) {
          added_attrs[i]=1;
          break;
        }
        attr_rec=0;
      }
    }
    /* If we have specific values for the given attribute, attr_rec contains
     * the appropriate vector. If not, it is null. */
    if (attr_rec) {
      for (i=0; i<nv; i++) {
        char *s;
        PyObject *o;
        switch (attr_rec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
          o=PyFloat_FromDouble((double)VECTOR(*(igraph_vector_t*)attr_rec->value)[i]);
          break;
        case IGRAPH_ATTRIBUTE_STRING:
          igraph_strvector_get((igraph_strvector_t*)attr_rec->value, i, &s);
          o=PyString_FromString(s);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string and not numeric)");
          o=0;
          break;
        }
        if (o) {
          if (PyList_Append(value, o) == -1)
            IGRAPH_ERROR("can't extend a vertex attribute hash member", IGRAPH_FAILURE);
          else Py_DECREF(o);
        }
      }
    } else {
      for (i=0; i<nv; i++) {
        if (PyList_Append(value, Py_None) == -1) {
          IGRAPH_ERROR("can't extend a vertex attribute hash member", IGRAPH_FAILURE);
        }
      }
    }
  }

  /* Okay, now we added the new attribute values for the already existing
   * attribute keys. Let's see if we have something left */
  if (attr) {
    l=igraph_vector_ptr_size(attr);
    j=igraph_vcount(graph)-nv;
    /* j contains the number of vertices EXCLUDING the new ones! */
    for (k=0; k<l; k++) {
      if (added_attrs[k]) continue;
      attr_rec=(igraph_i_attribute_record_t*)VECTOR(*attr)[k];

      value=PyList_New(j + nv);
      if (!value) {
        IGRAPH_ERROR("can't add attributes", IGRAPH_ENOMEM);
      }

      for (i=0; i<j; i++) {
        Py_INCREF(Py_None);
        PyList_SET_ITEM(value, i, Py_None);
      }

      for (i=0; i<nv; i++) {
        char *s;
        PyObject *o;
        switch (attr_rec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
          o=PyFloat_FromDouble((double)VECTOR(*(igraph_vector_t*)attr_rec->value)[i]);
          break;
        case IGRAPH_ATTRIBUTE_STRING:
          igraph_strvector_get((igraph_strvector_t*)attr_rec->value, i, &s);
          o=PyString_FromString(s);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string and not numeric)");
          o=0;
          break;
        }
        if (o) PyList_SET_ITEM(value, i+j, o);
      }

      PyDict_SetItemString(dict, attr_rec->name, value);
      Py_DECREF(value);   /* compensate for PyDict_SetItemString */
    }
    free(added_attrs);
    IGRAPH_FINALLY_CLEAN(1);
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
  Py_ssize_t pos=0;
  
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
      Py_INCREF(o);   /* take ownership, since PyList_SetItem will steal it */
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
  PyObject *key, *value, *dict;
  Py_ssize_t pos=0;
  long int i, j, k, l, ne;
  igraph_bool_t *added_attrs=0;
  igraph_i_attribute_record_t *attr_rec;

  ne=igraph_vector_size(edges)/2;
  if (!graph->attr) return IGRAPH_SUCCESS;
  if (ne<=0) return IGRAPH_SUCCESS;
  
  if (attr) {
    added_attrs = (igraph_bool_t*)calloc((size_t)igraph_vector_ptr_size(attr),
                                         sizeof(igraph_bool_t));
    if (!added_attrs)
      IGRAPH_ERROR("can't add vertex attributes", IGRAPH_ENOMEM);
    IGRAPH_FINALLY(free, added_attrs);
  }

  dict=((PyObject**)graph->attr)[2];
  if (!PyDict_Check(dict)) 
    IGRAPH_ERROR("edge attribute hash type mismatch", IGRAPH_EINVAL);
  while (PyDict_Next(dict, &pos, &key, &value)) {
    if (!PyString_Check(key))
      IGRAPH_ERROR("edge attribute hash key is not a string", IGRAPH_EINVAL);
    if (!PyList_Check(value))
      IGRAPH_ERROR("edge attribute hash member is not a list", IGRAPH_EINVAL);

    /* Check if we have specific values for the given attribute */
    attr_rec=0;
    if (attr) {
      j=igraph_vector_ptr_size(attr);
      for (i=0; i<j; i++) {
        attr_rec=VECTOR(*attr)[i];
        if (!strcmp(attr_rec->name, PyString_AS_STRING(key))) {
          added_attrs[i]=1;
          break;
        }
        attr_rec=0;
      }
    }
    /* If we have specific values for the given attribute, attr_rec contains
     * the appropriate vector. If not, it is null. */
    if (attr_rec) {
      for (i=0; i<ne; i++) {
        char *s;
        PyObject *o;
        switch (attr_rec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
          o=PyFloat_FromDouble((double)VECTOR(*(igraph_vector_t*)attr_rec->value)[i]);
          break;
        case IGRAPH_ATTRIBUTE_STRING:
          igraph_strvector_get((igraph_strvector_t*)attr_rec->value, i, &s);
          o=PyString_FromString(s);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string and not numeric)");
          o=0;
          break;
        }
        if (o) {
          if (PyList_Append(value, o) == -1)
            IGRAPH_ERROR("can't extend an edge attribute hash member", IGRAPH_FAILURE);
          else Py_DECREF(o);
        }
      }
    } else {
      for (i=0; i<ne; i++) {
        if (PyList_Append(value, Py_None) == -1) {
          IGRAPH_ERROR("can't extend an edge attribute hash member", IGRAPH_FAILURE);
        }
      }
    }
  }
  
  /*pos=0;
  while (PyDict_Next(dict, &pos, &key, &value)) {
    printf("key: "); PyObject_Print(key, stdout, Py_PRINT_RAW); printf("\n");
    printf("value: "); PyObject_Print(value, stdout, Py_PRINT_RAW); printf("\n");
  }*/
  
  /* Okay, now we added the new attribute values for the already existing
   * attribute keys. Let's see if we have something left */
  if (attr) {
    l=igraph_vector_ptr_size(attr);
    j=igraph_ecount(graph)-ne;
    /* j contains the number of edges EXCLUDING the new ones! */
    for (k=0; k<l; k++) {
      if (added_attrs[k]) continue;
      attr_rec=(igraph_i_attribute_record_t*)VECTOR(*attr)[k];

      value=PyList_New(j+ne);
      if (!value) {
        IGRAPH_ERROR("can't add attributes", IGRAPH_ENOMEM);
      }

      for (i=0; i<j; i++) {
        Py_INCREF(Py_None);
        PyList_SET_ITEM(value, i, Py_None);
      }

      for (i=0; i<ne; i++) {
        char *s;
        PyObject *o;
        switch (attr_rec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
          o=PyFloat_FromDouble((double)VECTOR(*(igraph_vector_t*)attr_rec->value)[i]);
          break;
        case IGRAPH_ATTRIBUTE_STRING:
          igraph_strvector_get((igraph_strvector_t*)attr_rec->value, i, &s);
          o=PyString_FromString(s);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string and not numeric)");
          o=0;
          break;
        }
        if (o) PyList_SET_ITEM(value, i+j, o);
      }

      PyDict_SetItemString(dict, attr_rec->name, value);
      Py_DECREF(value);   /* compensate for PyDict_SetItemString */
    }
    free(added_attrs);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return IGRAPH_SUCCESS;
}

/* Deleting edges */
static void igraphmodule_i_attribute_delete_edges(igraph_t *graph, const igraph_vector_t *idx) {
  long int n, i, ndeleted=0;
  PyObject *key, *value, *dict, *o;
  Py_ssize_t pos=0;
  
  dict=((PyObject**)graph->attr)[2];
  if (!PyDict_Check(dict)) return;

  n=igraph_vector_size(idx);
  for (i=0; i<n; i++) {
    /* printf("%ld:%f ", i, VECTOR(*idx)[i]); */
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

/* Permuting edges */
static int igraphmodule_i_attribute_permute_edges(igraph_t *graph,
						  const igraph_vector_t *idx) { 
  long int n, i;
  PyObject *key, *value, *dict, *newdict, *newlist, *o;
  Py_ssize_t pos=0;

  dict=((PyObject**)graph->attr)[2];
  if (!PyDict_Check(dict)) return 1;

  newdict=PyDict_New();
  if (!newdict) return 1;

  n=igraph_vector_size(idx);
  pos=0;

  while (PyDict_Next(dict, &pos, &key, &value)) {
    newlist=PyList_New(n);
    for (i=0; i<n; i++) {
      o=PyList_GetItem(value, VECTOR(*idx)[i]-1);
      if (!o) {
	PyErr_Clear();
	return 1;
      }
      Py_INCREF(o);
      PyList_SET_ITEM(newlist, i, o);
    }
    PyDict_SetItem(newdict, key, newlist);
    Py_DECREF(newlist);
  }

  ((PyObject**)graph->attr)[2]=newdict;
  Py_DECREF(dict);

  return 0;
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
  long int i, j, k, l, m;
  
  for (i=0; i<3; i++) {
    igraph_strvector_t *n = names[i];
    igraph_vector_t *t = types[i];
    PyObject *dict = ((PyObject**)graph->attr)[i];
    PyObject *keys;
    PyObject *values;
    PyObject *o=0;
    keys=PyDict_Keys(dict);
    if (!keys) IGRAPH_ERROR("Internal error in PyDict_Keys", IGRAPH_FAILURE);
 
    if (n) {
      j=igraphmodule_PyList_to_strvector_t(keys, n);
      if (j) return j;
    }
    if (t) {
      k=PyList_Size(keys);
      igraph_vector_init(t, k);
      for (j=0; j<k; j++) {
	int is_numeric = 1; 
	values=PyDict_GetItem(dict, PyList_GetItem(keys, j));
	if (PyList_Check(values)) {
	  m=PyList_Size(values);
	  for (l=0; l<m && is_numeric; l++) {
	    o=PyList_GetItem(values, l);
	    if (o != Py_None && !PyNumber_Check(o)) is_numeric=0;
	  }
	} else if (o != Py_None && !PyNumber_Check(values)) is_numeric=0;
      
	VECTOR(*t)[j]=is_numeric ? IGRAPH_ATTRIBUTE_NUMERIC : IGRAPH_ATTRIBUTE_STRING;
      }
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
  long int attrnum, i, j;
  int is_numeric;
  PyObject *o, *dict;
  switch (elemtype) {
  case IGRAPH_ATTRIBUTE_GRAPH: attrnum=0; break;
  case IGRAPH_ATTRIBUTE_VERTEX: attrnum=1; break;
  case IGRAPH_ATTRIBUTE_EDGE: attrnum=2; break;
  default: IGRAPH_ERROR("No such attribute type", IGRAPH_EINVAL); break;
  }
  dict = ((PyObject**)graph->attr)[attrnum];
  o = PyDict_GetItemString(dict, name);
  if (o == 0) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  is_numeric = 1;
  if (attrnum>0) {
    if (!PyList_Check(o)) IGRAPH_ERROR("attribute hash type mismatch", IGRAPH_EINVAL);
    if (!PyList_Size(o))  IGRAPH_ERROR("attribute hash type mismatch", IGRAPH_EINVAL);
    j = PyList_Size(o);
    for (i=0; i<j && is_numeric; i++) {
      PyObject *item = PyList_GET_ITEM(o, i);
      if (item != Py_None && !PyNumber_Check(item)) is_numeric=0;
    }
  } else if (o != Py_None && !PyNumber_Check(o)) is_numeric=0;
  if (is_numeric)
    *type = IGRAPH_ATTRIBUTE_NUMERIC;
  else
    *type = IGRAPH_ATTRIBUTE_STRING;
  return 0;
}

/* Getting numeric graph attributes */
int igraphmodule_i_get_numeric_graph_attr(const igraph_t *graph,
					  const char *name, igraph_vector_t *value) {
  PyObject *dict, *o, *result;
  dict = ((PyObject**)graph->attr)[0];
  /* No error checking, if we get here, the type has already been checked by previous
     attribute handler calls... hopefully :) Same applies for the other handlers. */
  o = PyDict_GetItemString(dict, name);
  if (!o) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  IGRAPH_CHECK(igraph_vector_resize(value, 1));
  if (o == Py_None) {
    VECTOR(*value)[0] = IGRAPH_NAN;
    return 0;
  }
  result = PyNumber_Float(o);
  if (result) {
    VECTOR(*value)[0] = PyFloat_AsDouble(o);
    Py_DECREF(result);
  } else IGRAPH_ERROR("Internal error in PyFloat_AsDouble", IGRAPH_EINVAL); 

  return 0;
}

/* Getting string graph attributes */
int igraphmodule_i_get_string_graph_attr(const igraph_t *graph,
					 const char *name, igraph_strvector_t *value) {
  PyObject *dict, *o, *result;
  dict = ((PyObject**)graph->attr)[0];
  o = PyDict_GetItemString(dict, name);
  if (!o) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  IGRAPH_CHECK(igraph_strvector_resize(value, 1));
  result = PyObject_Str(o);
  if (result) {
    IGRAPH_CHECK(igraph_strvector_set(value, 0, PyString_AsString(result)));
    Py_DECREF(result);
  } else IGRAPH_ERROR("Internal error in PyObject_Str", IGRAPH_EINVAL); 

  return 0;
}

/* Getting numeric vertex attributes */
int igraphmodule_i_get_numeric_vertex_attr(const igraph_t *graph,
					   const char *name,
					   igraph_vs_t vs,
					   igraph_vector_t *value) {
  PyObject *dict, *list, *result, *o;
  igraph_vector_t newvalue;

  dict = ((PyObject**)graph->attr)[1];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_vs_is_all(&vs)) {
    if (igraphmodule_PyObject_float_to_vector_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_vector_copy(value, &newvalue);
    igraph_vector_destroy(&newvalue);
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_VIT_SIZE(it)));
    while (!IGRAPH_VIT_END(it)) {
      long int v=IGRAPH_VIT_GET(it);
      o = PyList_GetItem(list, v);
      if (o != Py_None) {
	result = PyNumber_Float(o);
	VECTOR(*value)[i] = PyFloat_AsDouble(result);
	Py_XDECREF(result);
      } else VECTOR(*value)[i] = IGRAPH_NAN;
      IGRAPH_VIT_NEXT(it);
      i++;
    }
    igraph_vit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

/* Getting string vertex attributes */
int igraphmodule_i_get_string_vertex_attr(const igraph_t *graph,
					  const char *name,
					  igraph_vs_t vs,
					  igraph_strvector_t *value) {
  PyObject *dict, *list, *result;
  igraph_strvector_t newvalue;

  dict = ((PyObject**)graph->attr)[1];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_vs_is_all(&vs)) {
    if (igraphmodule_PyList_to_strvector_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_strvector_destroy(value);
    *value=newvalue;
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_VIT_SIZE(it)));
    while (!IGRAPH_VIT_END(it)) {
      long int v=IGRAPH_VIT_GET(it);
      result = PyObject_Str(PyList_GetItem(list, v));
      igraph_strvector_set(value, i, PyString_AsString(result));
      Py_XDECREF(result);
      IGRAPH_VIT_NEXT(it);
      i++;
    }
    igraph_vit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

/* Getting numeric edge attributes */
int igraphmodule_i_get_numeric_edge_attr(const igraph_t *graph,
					 const char *name,
					 igraph_es_t es,
					 igraph_vector_t *value) {
  PyObject *dict, *list, *result, *o;
  igraph_vector_t newvalue;

  dict = ((PyObject**)graph->attr)[2];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_es_is_all(&es)) {
    if (igraphmodule_PyObject_float_to_vector_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_vector_copy(value, &newvalue);
    igraph_vector_destroy(&newvalue);
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_EIT_SIZE(it)));
    while (!IGRAPH_EIT_END(it)) {
      long int v=IGRAPH_EIT_GET(it);
      o = PyList_GetItem(list, v);
      if (o != Py_None) {
        result = PyNumber_Float(o);
        VECTOR(*value)[i] = PyFloat_AsDouble(result);
        Py_XDECREF(result);
      } else VECTOR(*value)[i] = IGRAPH_NAN;
      IGRAPH_EIT_NEXT(it);
      i++;
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

/* Getting string edge attributes */
int igraphmodule_i_get_string_edge_attr(const igraph_t *graph,
					const char *name,
					igraph_es_t es,
					igraph_strvector_t *value) {
  PyObject *dict, *list, *result;
  igraph_strvector_t newvalue;

  dict = ((PyObject**)graph->attr)[2];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_es_is_all(&es)) {
    if (igraphmodule_PyList_to_strvector_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_strvector_destroy(value);
    *value=newvalue;
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_EIT_SIZE(it)));
    while (!IGRAPH_EIT_END(it)) {
      long int v=IGRAPH_EIT_GET(it);
      result = PyObject_Str(PyList_GetItem(list, v));
      igraph_strvector_set(value, i, PyString_AsString(result));
      Py_XDECREF(result);
      IGRAPH_EIT_NEXT(it);
      i++;
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

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
  igraphmodule_i_attribute_permute_edges,
  igraphmodule_i_attribute_get_info,
  igraphmodule_i_attribute_has_attr,
  igraphmodule_i_attribute_get_type,
  igraphmodule_i_get_numeric_graph_attr,
  igraphmodule_i_get_string_graph_attr,
  igraphmodule_i_get_numeric_vertex_attr,
  igraphmodule_i_get_string_vertex_attr,
  igraphmodule_i_get_numeric_edge_attr,
  igraphmodule_i_get_string_edge_attr,
};

/** \ingroup python_interface
 * \brief Method table for the igraph Python module
 */
static PyMethodDef igraphmodule_methods[] = 
{
  {"community_to_membership", (PyCFunction)igraphmodule_community_to_membership,
    METH_VARARGS | METH_KEYWORDS,
    "community_to_membership(merges, nodes, steps, return_csize=False)\n\n"
  },
  {"convex_hull", (PyCFunction)igraphmodule_convex_hull, METH_VARARGS,
      "convex_hull(vs, coords=False)\n\n"
      "Calculates the convex hull of a given point set.\n\n"
      "@param vs: the point set as a list of lists\n"
      "@param coords: if C{True}, the function returns the\n"
      "  coordinates of the corners of the convex hull polygon,\n"
      "  otherwise returns the corner indices.\n"
      "@return: either the hull's corner coordinates or the point\n"
      "  indices corresponding to them, depending on the C{coords}\n"
      "  parameter."
  },
  {"set_progress_handler", igraphmodule_set_progress_handler, METH_VARARGS,
      "set_progress_handler(handler)\n\n"
      "Sets the handler to be called when igraph is performing a long operation.\n"
      "@param handler: the progress handler function. It must accept two\n"
      "  arguments, the first is the message informing the user about\n"
      "  what igraph is doing right now, the second is the actual\n"
      "  progress information (a percentage).\n"
  },
  {NULL, NULL, 0, NULL}
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

extern PyObject* igraphmodule_InternalError;
extern PyObject* igraphmodule_arpack_options_default;

PyMODINIT_FUNC initcore(void) {
  PyObject* m;
  
  if (PyType_Ready(&igraphmodule_VertexSeqType) < 0) return;
  if (PyType_Ready(&igraphmodule_EdgeSeqType) < 0) return;
  
  igraphmodule_VertexType.tp_clear = (inquiry)igraphmodule_Vertex_clear;
  if (PyType_Ready(&igraphmodule_VertexType) < 0) return;
  
  igraphmodule_EdgeType.tp_clear = (inquiry)igraphmodule_Edge_clear;
  if (PyType_Ready(&igraphmodule_EdgeType) < 0) return;
  
  if (PyType_Ready(&igraphmodule_GraphType) < 0) return;
  if (PyType_Ready(&igraphmodule_BFSIterType) < 0) return;
  if (PyType_Ready(&igraphmodule_ARPACKOptionsType) < 0) return;

  m = Py_InitModule3("igraph.core", igraphmodule_methods,
		     "Low-level Python interface for the igraph library. "
		     "Should not be used directly.");
  
  PyModule_AddObject(m, "GraphBase", (PyObject*)&igraphmodule_GraphType);
  PyModule_AddObject(m, "BFSIter", (PyObject*)&igraphmodule_BFSIterType);
  PyModule_AddObject(m, "ARPACKOptions", (PyObject*)&igraphmodule_ARPACKOptionsType);
  PyModule_AddObject(m, "Edge", (PyObject*)&igraphmodule_EdgeType);
  PyModule_AddObject(m, "EdgeSeq", (PyObject*)&igraphmodule_EdgeSeqType);
  PyModule_AddObject(m, "Vertex", (PyObject*)&igraphmodule_VertexType);
  PyModule_AddObject(m, "VertexSeq", (PyObject*)&igraphmodule_VertexSeqType);
 
  /* Internal error exception type */
  igraphmodule_InternalError =
    PyErr_NewException("igraph.core.InternalError", PyExc_Exception, NULL);
  PyModule_AddObject(m, "InternalError", igraphmodule_InternalError);

  /* ARPACK default options variable */
  igraphmodule_arpack_options_default = igraphmodule_ARPACKOptions_new();
  PyModule_AddObject(m, "arpack_options", igraphmodule_arpack_options_default);

  /* Useful constants */
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

  PyModule_AddIntConstant(m, "BLISS_F", IGRAPH_BLISS_F);
  PyModule_AddIntConstant(m, "BLISS_FL", IGRAPH_BLISS_FL);
  PyModule_AddIntConstant(m, "BLISS_FS", IGRAPH_BLISS_FS);
  PyModule_AddIntConstant(m, "BLISS_FM", IGRAPH_BLISS_FM);
  PyModule_AddIntConstant(m, "BLISS_FLM", IGRAPH_BLISS_FLM);
  PyModule_AddIntConstant(m, "BLISS_FSM", IGRAPH_BLISS_FSM);

  /* More useful constants */
  PyModule_AddStringConstant(m, "__version__", "0.5.2");
  PyModule_AddStringConstant(m, "__build_date__", __DATE__);

  /* initialize error, progress, warning and interruption handler */
  igraph_set_error_handler(igraphmodule_igraph_error_hook);
  igraph_set_progress_handler(igraphmodule_igraph_progress_hook);
  igraph_set_warning_handler(igraphmodule_igraph_warning_hook);
  igraph_set_interruption_handler(igraphmodule_igraph_interrupt_hook);
  
  /* initialize attribute handlers */
  igraph_i_set_attribute_table(&igraphmodule_i_attribute_table);
}
