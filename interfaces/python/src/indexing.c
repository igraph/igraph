/* vim:set ts=4 sw=2 sts=2 et:  */
/* 
   IGraph library - Python interface.
   Copyright (C) 2006-2011  Tamas Nepusz <ntamas@gmail.com>
   5 Avenue Road, Staines, Middlesex, TW18 3AW, United Kingdom
   
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
#include "indexing.h"
#include "platform.h"
#include "py2compat.h"
#include "pyhelpers.h"

/***************************************************************************/

static PyObject* igraphmodule_i_Graph_adjmatrix_indexing_get_value_for_vertex_pair(
    igraph_t* graph, igraph_integer_t from, igraph_integer_t to, PyObject* values) {
  igraph_integer_t eid;
  PyObject* result;

  /* Retrieving a single edge */
  igraph_get_eid(graph, &eid, from, to, /* directed = */1, /* error = */0);
  if (eid >= 0) {
    /* Edge found, get the value of the attribute */
    if (values == 0) {
      return PyInt_FromLong(1L);
    } else {
      result = PyList_GetItem(values, eid);
      Py_XINCREF(result);
      return result;
    }
  } else {
    /* No such edge, return zero */
    return PyInt_FromLong(0L);
  }
}

static PyObject* igraphmodule_i_Graph_adjmatrix_get_index_row(igraph_t* graph, 
    igraph_integer_t from, igraph_vs_t* to, igraph_neimode_t neimode,
    PyObject* values);

PyObject* igraphmodule_Graph_adjmatrix_get_index(igraph_t* graph,
        PyObject* row_index, PyObject* column_index, PyObject* attr_name) {
    PyObject *result = 0, *values;
    igraph_vs_t vs1, vs2;
    igraph_integer_t vid1 = -1, vid2 = -1;
    char* attr;

    if (igraphmodule_PyObject_to_vs_t(row_index, &vs1, graph, 0, &vid1))
      return NULL;
    if (igraphmodule_PyObject_to_vs_t(column_index, &vs2, graph, 0, &vid2))
      return NULL;

    if (attr_name == 0) {
      /* Using the "weight" attribute by default */
      values = igraphmodule_get_edge_attribute_values(graph, "weight");
    } else {
      /* Specifying the name of the attribute */
      attr = PyObject_ConvertToCString(attr_name);
      values = igraphmodule_get_edge_attribute_values(graph, attr);
      free(attr);
    }

    if (vid1 >= 0 && vid2 >= 0) {
      /* Retrieving an edge between vid1 and vid2 */
      result = igraphmodule_i_Graph_adjmatrix_indexing_get_value_for_vertex_pair(
          graph, vid1, vid2, values);
    } else if (vid1 >= 0) {
      /* Retrieving the successors of vid1 */
      result = igraphmodule_i_Graph_adjmatrix_get_index_row(
          graph, vid1, &vs2, IGRAPH_OUT, values);
    } else if (vid2 >= 0) {
      /* Retrieving the predecessors of vid2 */
      result = igraphmodule_i_Graph_adjmatrix_get_index_row(
          graph, vid2, &vs1, IGRAPH_IN, values);
    } else {
      /* Retrieving a submatrix */
      igraph_vit_t vit;
      PyObject *item;

      if (igraph_vit_create(graph, vs1, &vit)) {
        igraphmodule_handle_igraph_error();
        result = 0;
      } else {
        result = PyList_New(0);
        if (result != 0) {
          while (!IGRAPH_VIT_END(vit)) {
            vid1 = IGRAPH_VIT_GET(vit);
            item = igraphmodule_i_Graph_adjmatrix_get_index_row(graph, vid1, &vs2, IGRAPH_OUT, values);
            if (item == 0) {
              Py_DECREF(result);
              result = 0;
              break;
            }
            if (PyList_Append(result, item)) {
              /* error while appending */
              Py_DECREF(item);
              Py_DECREF(result);
              result = 0;
              break;
            }
            Py_DECREF(item);
            IGRAPH_VIT_NEXT(vit);
          }
        }
        igraph_vit_destroy(&vit);
      }
    }

    igraph_vs_destroy(&vs1);
    igraph_vs_destroy(&vs2);

    return result;
}

static PyObject* igraphmodule_i_Graph_adjmatrix_get_index_row(igraph_t* graph, 
    igraph_integer_t from, igraph_vs_t* to, igraph_neimode_t neimode,
    PyObject* values) {
  igraph_vector_t eids;
  igraph_integer_t eid;
  igraph_vit_t vit;
  PyObject *result = 0, *item;
  long int i, n;
  igraph_integer_t v;

  if (igraph_vs_is_all(to)) {
    /* Simple case: all edges */
    IGRAPH_PYCHECK(igraph_vector_init(&eids, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &eids);
    IGRAPH_PYCHECK(igraph_incident(graph, &eids, from, neimode));
    
    n = igraph_vector_size(&eids);
    result = PyList_Zeroes(igraph_vcount(graph));
    if (result == 0) {
      IGRAPH_FINALLY_FREE();
      return 0;
    }

    for (i = 0; i < n; i++) {
      eid = (igraph_integer_t)VECTOR(eids)[i];
      v = IGRAPH_OTHER(graph, eid, from);
      if (values)
        item = PyList_GetItem(values, eid);
      else
        item = PyInt_FromLong(1);
      Py_INCREF(item);
      PyList_SetItem(result, v, item);   /* reference stolen here */
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_vector_destroy(&eids);

    return result;
  }

  /* More complicated case: only some vertices */
  IGRAPH_PYCHECK(igraph_vit_create(graph, *to, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  result = PyList_New(0);
  if (result == 0) {
    IGRAPH_FINALLY_FREE();
    return 0;
  }

  while (!IGRAPH_VIT_END(vit)) {
    v = IGRAPH_VIT_GET(vit);
    if (neimode == IGRAPH_OUT) {
      item = igraphmodule_i_Graph_adjmatrix_indexing_get_value_for_vertex_pair(
          graph, from, v, values);
    } else {
      item = igraphmodule_i_Graph_adjmatrix_indexing_get_value_for_vertex_pair(
          graph, v, from, values);
    }
    if (item == 0) {
      IGRAPH_FINALLY_FREE();
      Py_DECREF(result);
      return 0;
    }
    if (PyList_Append(result, item)) {
      /* error while appending */
      Py_DECREF(item);
      Py_DECREF(result);
      result = 0;
      break;
    }
    Py_DECREF(item);
    IGRAPH_VIT_NEXT(vit);
  }

  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(1);

  return result;
}

/***************************************************************************/

/**
 * Determines whether the given Python value means that the user would like
 * to delete the edge the value is being assigned to in the adjacency matrix
 * assignment syntax.
 */
static INLINE igraph_bool_t deleting_edge(PyObject* value) {
  return value == Py_None || value == Py_False ||
      (PyInt_Check(value) && PyInt_AsLong(value) == 0);
}

/**
 * Structure to hold data related to newly added/removed edges during an
 * adjacency matrix assignment.
 */
typedef struct {
  igraph_vector_t to_add;
  PyObject* to_add_values;
  igraph_vector_t to_delete;
} igraphmodule_i_Graph_adjmatrix_set_index_data_t;

int igraphmodule_i_Graph_adjmatrix_set_index_data_init(
    igraphmodule_i_Graph_adjmatrix_set_index_data_t* data) {
  if (igraph_vector_init(&data->to_add, 0)) {
    igraphmodule_handle_igraph_error();
    return -1;
  }

  if (igraph_vector_init(&data->to_delete, 0)) {
    igraphmodule_handle_igraph_error();
    igraph_vector_destroy(&data->to_delete);
    return -1;
  }

  data->to_add_values = PyList_New(0);
  if (data->to_add_values == 0) {
    igraph_vector_destroy(&data->to_add);
    igraph_vector_destroy(&data->to_delete);
    return -1;
  }

  return 0;
}

void igraphmodule_i_Graph_adjmatrix_set_index_data_destroy(
    igraphmodule_i_Graph_adjmatrix_set_index_data_t* data) {
  igraph_vector_destroy(&data->to_add);
  igraph_vector_destroy(&data->to_delete);
  Py_DECREF(data->to_add_values);
}

static int igraphmodule_i_Graph_adjmatrix_set_index_row(igraph_t* graph, 
    igraph_integer_t from, igraph_vs_t* to, igraph_neimode_t neimode,
    PyObject* values, PyObject* new_value,
    igraphmodule_i_Graph_adjmatrix_set_index_data_t* data) {
  PyObject *iter = 0, *item;
  igraph_vit_t vit;
  igraph_integer_t v, v1, v2, eid;
  igraph_bool_t deleting, ok = 1;

  /* Check whether new_value is an iterable (and not a string). If not,
   * every assignment will use the same value (that is, new_value) */
  if (!PyBaseString_Check(new_value)) {
    iter = PyObject_GetIter(new_value);
    if (PyErr_Occurred()) {
      /* Object is not an iterable. Clear the exception */
      iter = 0;
      PyErr_Clear();
    }
  }

  if (igraph_vit_create(graph, *to, &vit)) {
    Py_XDECREF(iter);
    igraphmodule_handle_igraph_error();
    return -1;
  }

  v1 = from; v2 = from;

  /* The two branches of the next `if' are almost the same; make sure
   * you make changes to both branches if appropriate! */

  if (iter != 0) {
    /* The new value is an iterable, so it must have exactly as many elements
     * as the number of vertices in the graph. If it has less, we simply
     * skip the rest (with a warning) */
    while (!IGRAPH_VIT_END(vit) && (item = PyIter_Next(iter)) != 0) {
      v = IGRAPH_VIT_GET(vit);

      /* Get the ID of the edge between from and v */
      if (neimode == IGRAPH_OUT) {
        v2 = v;
      } else {
        v1 = v;
      }
      igraph_get_eid(graph, &eid, v1, v2, /* directed = */1, /* error = */0);
      if (deleting_edge(item)) {
        /* Deleting edges if eid != -1 */
        if (eid != -1) {
          if (igraph_vector_push_back(&data->to_delete, eid)) {
            igraphmodule_handle_igraph_error();
            igraph_vector_clear(&data->to_delete);
            ok = 0;
            break;
          }
        }
      } else {
        if (eid == -1) {
          /* Adding edges */
          if (igraph_vector_push_back(&data->to_add, v1) ||
              igraph_vector_push_back(&data->to_add, v2)) {
            igraphmodule_handle_igraph_error();
            igraph_vector_clear(&data->to_add);
            ok = 0;
            break;
          }
          if (values != 0) {
            Py_INCREF(new_value);
            if (PyList_Append(data->to_add_values, new_value)) {
              Py_DECREF(new_value);
              igraph_vector_clear(&data->to_add);
              ok = 0;
              break;
            }
          }
        } else if (values != 0) {
          /* Setting attribute */
          Py_INCREF(item);
          if (PyList_SetItem(values, eid, item)) {
            Py_DECREF(item);
            igraph_vector_clear(&data->to_add);
          }
        }
      }
      Py_DECREF(item);
      IGRAPH_VIT_NEXT(vit);
    }
    if (!IGRAPH_VIT_END(vit)) {
      PyErr_WarnEx(PyExc_RuntimeWarning,
          "iterable was shorter than the number of vertices in the vertex "
          "sequence", 1);
    }
  } else {
    /* The new value is not an iterable; setting the same value for
     * more than one edge */
    deleting = deleting_edge(new_value);
    while (!IGRAPH_VIT_END(vit)) {
      v = IGRAPH_VIT_GET(vit);

      /* Get the ID of the edge between from and v */
      if (neimode == IGRAPH_OUT) {
        v2 = v;
      } else {
        v1 = v;
      }
      igraph_get_eid(graph, &eid, v1, v2, /* directed = */1, /* error = */0);

      if (deleting) {
        /* Deleting edges if eid != -1 */
        if (eid != -1) {
          if (igraph_vector_push_back(&data->to_delete, eid)) {
            igraphmodule_handle_igraph_error();
            igraph_vector_clear(&data->to_delete);
            ok = 0;
            break;
          }
        }
      } else {
        if (eid == -1) {
          /* Adding edges */
          if (igraph_vector_push_back(&data->to_add, v1) ||
              igraph_vector_push_back(&data->to_add, v2)) {
            igraphmodule_handle_igraph_error();
            igraph_vector_clear(&data->to_add);
            ok = 0;
            break;
          }
          if (values != 0) {
            Py_INCREF(new_value);
            if (PyList_Append(data->to_add_values, new_value)) {
              Py_DECREF(new_value);
              igraph_vector_clear(&data->to_add);
              ok = 0;
              break;
            }
          }
        } else if (values != 0) {
          /* Setting attribute */
          Py_INCREF(new_value);
          if (PyList_SetItem(values, eid, new_value)) {
            Py_DECREF(new_value);
            igraph_vector_clear(&data->to_add);
          }
        }
      }
      IGRAPH_VIT_NEXT(vit);
    }
  }

  Py_XDECREF(iter);
  igraph_vit_destroy(&vit);

  return ok ? 0 : -1;
}

int igraphmodule_Graph_adjmatrix_set_index(igraph_t* graph,
        PyObject* row_index, PyObject* column_index, PyObject* attr_name,
        PyObject* new_value) {
  PyObject *values;
  igraph_vs_t vs1, vs2;
  igraph_vit_t vit;
  igraph_integer_t vid1 = -1, vid2 = -1, eid = -1;
  igraph_bool_t ok = 1;
  igraphmodule_i_Graph_adjmatrix_set_index_data_t data;
  char* attr;

  if (igraphmodule_PyObject_to_vs_t(row_index, &vs1, graph, 0, &vid1))
    return -1;
  if (igraphmodule_PyObject_to_vs_t(column_index, &vs2, graph, 0, &vid2))
    return -1;

  if (attr_name == 0) {
    /* Using the "weight" attribute by default */
    values = igraphmodule_get_edge_attribute_values(graph, "weight");
  } else {
    /* Specifying the name of the attribute */
    attr = PyObject_ConvertToCString(attr_name);
    values = igraphmodule_create_or_get_edge_attribute_values(graph, attr);
    free(attr);
  }

  if (vid1 >= 0 && vid2 >= 0) {
    /* Setting an edge between vid1 and vid2 */
    igraph_get_eid(graph, &eid, vid1, vid2, /* directed = */1, /* error = */0);
    if (deleting_edge(new_value)) {
      if (eid != -1) {
        /* Deleting the edge between vid1 and vid2 if it is there */
        if (igraph_delete_edges(graph, igraph_ess_1(eid))) {
          igraphmodule_handle_igraph_error();
          ok = 0;
        }
      }
    } else {
      /* Adding the edge between vid1 and vid2 if it is not there */
      if (eid == -1) {
        eid = igraph_ecount(graph);
        if (igraph_add_edge(graph, vid1, vid2)) {
          igraphmodule_handle_igraph_error();
          ok = 0;
        }
      }
      if (ok && values != 0) {
        /* Set the attribute value */
        Py_INCREF(new_value);
        PyList_SetItem(values, eid, new_value); /* reference stolen here */
      }
    }
  } else {
    /* In all the non-trivial cases, we do the modifications in three phases;
     * in the first phase, we modify the attribute values of edges that are to
     * stay (but possibly with a different attribute value) and collect the
     * list of edges to be added (and their attribute values) and the list of
     * edge to be deleted. In the second phase, we do the deletions in one
     * batch. Finally, we add the edges to be added.
     */
    igraphmodule_i_Graph_adjmatrix_set_index_data_init(&data);

    /* First phase */
    if (vid1 >= 0) {
      /* vs1 is a single vertex, vs2 is not */
      ok = (igraphmodule_i_Graph_adjmatrix_set_index_row(
              graph, vid1, &vs2, IGRAPH_OUT, values, new_value, &data) == 0);
    } else if (vid2 >= 0) {
      /* vs2 is a single vertex, vs1 is not */
      ok = (igraphmodule_i_Graph_adjmatrix_set_index_row(
              graph, vid2, &vs1, IGRAPH_IN, values, new_value, &data) == 0);
    } else {
      /* Complete submatrix */
      if (igraph_vit_create(graph, vs1, &vit)) {
        igraphmodule_handle_igraph_error();
        ok = 0;
      } else {
        while (!IGRAPH_VIT_END(vit)) {
          vid1 = IGRAPH_VIT_GET(vit);
          if (igraphmodule_i_Graph_adjmatrix_set_index_row(
                graph, vid1, &vs2, IGRAPH_OUT, values, new_value, &data) == 0) {
            ok = 0;
            break;
          }
          IGRAPH_VIT_NEXT(vit);
        }
        igraph_vit_destroy(&vit);
      }
    }

    if (ok) {
      /* Second phase: do the deletions in one batch */
      if (igraph_delete_edges(graph, igraph_ess_vector(&data.to_delete))) {
        igraphmodule_handle_igraph_error();
        ok = 0;
      }
    }

    if (ok) {
      /* Third phase: add the new edges in one batch */
      if (!igraph_vector_empty(&data.to_add)) {
        eid = igraph_ecount(graph);
        igraph_add_edges(graph, &data.to_add, 0);
        if (values != 0) {
          PyList_SetSlice(values, eid, eid+PyList_Size(data.to_add_values),
              data.to_add_values);
          if (PyList_Size(values) != igraph_ecount(graph)) {
            PyErr_SetString(PyExc_ValueError, "hmmm, attribute value list "
                "length mismatch, this is most likely a bug.");
            ok = 0;
          }
        }
      }
    }

    igraphmodule_i_Graph_adjmatrix_set_index_data_destroy(&data);
  }

  igraph_vs_destroy(&vs1);
  igraph_vs_destroy(&vs2);

  return ok ? 0 : -1;
}

