/* vim:set ts=2 sw=2 sts=2 et: */
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

#include <Python.h>
#include "attributes.h"
#include "common.h"
#include "convert.h"
#include "py2compat.h"
#include "pyhelpers.h"

int igraphmodule_i_attribute_struct_init(igraphmodule_i_attribute_struct *attrs) {
  int i;
  for (i=0; i<3; i++) {
    attrs->attrs[i] = PyDict_New();
    if (PyErr_Occurred())
      return 1;
    RC_ALLOC("dict", attrs->attrs[i]);
  }
  attrs->vertex_name_index = 0;
  return 0;
}

void igraphmodule_i_attribute_struct_destroy(igraphmodule_i_attribute_struct *attrs) {
  int i;
  for (i=0; i<3; i++) {
    if (attrs->attrs[i]) {
      RC_DEALLOC("dict", attrs->attrs[i]);
      Py_DECREF(attrs->attrs[i]);
    }
  }
  if (attrs->vertex_name_index) {
    RC_DEALLOC("dict", attrs->vertex_name_index);
    Py_DECREF(attrs->vertex_name_index);
  }
}

int igraphmodule_i_attribute_struct_index_vertex_names(
    igraphmodule_i_attribute_struct *attrs, igraph_bool_t force) {
  Py_ssize_t n = 0;
  PyObject *name_list, *key, *value;

  if (attrs->vertex_name_index && !force)
    return 0;

  if (attrs->vertex_name_index == 0) {
    attrs->vertex_name_index = PyDict_New();
    if (attrs->vertex_name_index == 0) {
      return 1;
    }
  } else
    PyDict_Clear(attrs->vertex_name_index);

  name_list = PyDict_GetItemString(attrs->attrs[1], "name");
  if (name_list == 0)
    return 0;    /* no name attribute */

  n = PyList_Size(name_list) - 1;
  while (n >= 0) {
    key = PyList_GET_ITEM(name_list, n);    /* we don't own a reference to key */
    value = PyInt_FromLong(n);              /* we do own a reference to value */
    if (value == 0)
      return 1;
    PyDict_SetItem(attrs->vertex_name_index, key, value);
    /* PyDict_SetItem did an INCREF for both the key and a value, therefore we
     * have to drop our reference on value */
    Py_DECREF(value);

    n--;
  }

  return 0;
}

void igraphmodule_i_attribute_struct_invalidate_vertex_name_index(
    igraphmodule_i_attribute_struct *attrs) {
  if (attrs->vertex_name_index == 0)
    return;

  Py_DECREF(attrs->vertex_name_index);
  attrs->vertex_name_index = 0;
}

void igraphmodule_invalidate_vertex_name_index(igraph_t *graph) {
  igraphmodule_i_attribute_struct_invalidate_vertex_name_index(ATTR_STRUCT(graph));
}

void igraphmodule_index_vertex_names(igraph_t *graph, igraph_bool_t force) {
  igraphmodule_i_attribute_struct_index_vertex_names(ATTR_STRUCT(graph), force);
}

int igraphmodule_get_vertex_id_by_name(igraph_t *graph, PyObject* o, igraph_integer_t* vid) {
  igraphmodule_i_attribute_struct* attrs = ATTR_STRUCT(graph);
  PyObject* o_vid = NULL;
  int tmp;

  if (graph) {
    attrs = ATTR_STRUCT(graph);
    if (igraphmodule_i_attribute_struct_index_vertex_names(attrs, 0))
      return 1;
    o_vid = PyDict_GetItem(attrs->vertex_name_index, o);
  }

  if (o_vid == NULL) {
#ifdef IGRAPH_PYTHON3
    PyErr_Format(PyExc_ValueError, "no such vertex: %R", o);
#else
    PyObject* s = PyObject_Repr(o);
    if (s) {
      PyErr_Format(PyExc_ValueError, "no such vertex: %s", PyString_AS_STRING(s));
      Py_DECREF(s);
    } else {
      PyErr_Format(PyExc_ValueError, "no such vertex: %p", o);
    }
#endif
    return 1;
  }

  if (!PyInt_Check(o_vid)) {
    PyErr_SetString(PyExc_ValueError, "non-numeric vertex ID assigned to vertex name. This is most likely a bug.");
    return 1;
  }
  
  if (PyInt_AsInt(o_vid, &tmp))
    return 1;
  
  *vid = tmp;

  return 0;
}

/**
 * \brief Checks whether the given graph has the given graph attribute.
 *
 * \param  graph  the graph
 * \param  name   the name of the attribute being searched for
 */
igraph_bool_t igraphmodule_has_graph_attribute(const igraph_t *graph, const char* name) {
  PyObject *dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_GRAPH];
  return name != 0 && dict != 0 && PyDict_GetItemString(dict, name) != 0;
}

/**
 * \brief Checks whether the given graph has the given vertex attribute.
 *
 * \param  graph  the graph
 * \param  name   the name of the attribute being searched for
 */
igraph_bool_t igraphmodule_has_vertex_attribute(const igraph_t *graph, const char* name) {
  PyObject *dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_VERTEX];
  return name != 0 && dict != 0 && PyDict_GetItemString(dict, name) != 0;
}

/**
 * \brief Checks whether the given graph has the given edge attribute.
 *
 * \param  graph  the graph
 * \param  name   the name of the attribute being searched for
 */
igraph_bool_t igraphmodule_has_edge_attribute(const igraph_t *graph, const char* name) {
  PyObject *dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  return name != 0 && dict != 0 && PyDict_GetItemString(dict, name) != 0;
}

/**
 * \brief Creates a new edge attribute and sets the values to None.
 *
 * This returns the actual list that we use to store the edge attributes, so
 * be careful when modifying it - any modification will propagate back to the
 * graph itself. You have been warned.
 *
 * \param  graph  the graph
 * \param  name   the name of the attribute being created
 * \returns  a Python list of the values or \c NULL if the given
 *           attribute exists already (no exception set). The returned
 *           reference is borrowed.
 */
PyObject* igraphmodule_create_edge_attribute(const igraph_t* graph,
    const char* name) {
  PyObject *dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  PyObject *values;
  Py_ssize_t i, n;

  if (dict == 0) {
    dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE] = PyDict_New();
  }
  if (PyDict_GetItemString(dict, name))
    return 0;

  n = igraph_ecount(graph);
  values = PyList_New(n);
  if (values == 0)
    return 0;

  for (i = 0; i < n; i++) {
    Py_INCREF(Py_None);
    PyList_SET_ITEM(values, i, Py_None);   /* reference stolen */
  }

  if (PyDict_SetItemString(dict, name, values)) {
    Py_DECREF(values);
    return 0;
  }

  Py_DECREF(values);
  return values;
}

/**
 * \brief Returns the values of the given edge attribute for all edges in the
 *        given graph.
 *
 * This returns the actual list that we use to store the edge attributes, so
 * be careful when modifying it - any modification will propagate back to the
 * graph itself. You have been warned.
 *
 * \param  graph  the graph
 * \param  name   the name of the attribute being searched for
 * \returns  a Python list or \c NULL if there is no such attribute
 *           (no exception set). The returned reference is borrowed.
 */
PyObject* igraphmodule_get_edge_attribute_values(const igraph_t* graph,
    const char* name) {
  PyObject *dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  if (dict == 0)
    return 0;
  return PyDict_GetItemString(dict, name);
}

/**
 * \brief Returns the values of the given edge attribute for all edges in the
 *        given graph, optionally creating it if it does not exist.
 *
 * This returns the actual list that we use to store the edge attributes, so
 * be careful when modifying it - any modification will propagate back to the
 * graph itself. You have been warned.
 *
 * \param  graph  the graph
 * \param  name   the name of the attribute being searched for
 * \returns  a Python list (borrowed reference)
 */
PyObject* igraphmodule_create_or_get_edge_attribute_values(const igraph_t* graph,
    const char* name) {
  PyObject *dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE], *result;
  if (dict == 0)
    return 0;
  result = PyDict_GetItemString(dict, name);
  if (result != 0)
    return result;
  return igraphmodule_create_edge_attribute(graph, name);
}

/* Attribute handlers for the Python interface */

/* Initialization */ 
static int igraphmodule_i_attribute_init(igraph_t *graph, igraph_vector_ptr_t *attr) {
  igraphmodule_i_attribute_struct* attrs;
  long int i, n;
  
  attrs=(igraphmodule_i_attribute_struct*)calloc(1, sizeof(igraphmodule_i_attribute_struct));
  if (!attrs)
    IGRAPH_ERROR("not enough memory to allocate attribute hashes", IGRAPH_ENOMEM);
  if (igraphmodule_i_attribute_struct_init(attrs)) {
    PyErr_Clear();
    free(attrs);
    IGRAPH_ERROR("not enough memory to allocate attribute hashes", IGRAPH_ENOMEM);
  }
  graph->attr=(void*)attrs;

  /* See if we have graph attributes */
  if (attr) {
    PyObject *dict=attrs->attrs[0], *value;
    char *s;
    n = igraph_vector_ptr_size(attr);
    for (i=0; i<n; i++) {
      igraph_attribute_record_t *attr_rec;
      attr_rec = VECTOR(*attr)[i];
      switch (attr_rec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
        value=PyFloat_FromDouble((double)VECTOR(*(igraph_vector_t*)attr_rec->value)[0]);
        break;
      case IGRAPH_ATTRIBUTE_STRING:
        igraph_strvector_get((igraph_strvector_t*)attr_rec->value, 0, &s);
        if (s == 0)
          value=PyString_FromString("");
        else
          value=PyString_FromString(s);
        break;
      case IGRAPH_ATTRIBUTE_BOOLEAN:
        value=VECTOR(*(igraph_vector_bool_t*)attr_rec->value)[0] ? Py_True : Py_False;
        Py_INCREF(value);
        break;
      default:
        IGRAPH_WARNING("unsupported attribute type (not string, numeric or Boolean)");
        value=0;
        break;
      }
      if (value) {
        if (PyDict_SetItemString(dict, attr_rec->name, value)) {
          Py_DECREF(value);
          igraphmodule_i_attribute_struct_destroy(attrs);
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
  igraphmodule_i_attribute_struct* attrs;
 
  /* printf("Destroying attribute table\n"); */
  if (graph->attr) {
    attrs=(igraphmodule_i_attribute_struct*)graph->attr;
    igraphmodule_i_attribute_struct_destroy(attrs);
    free(attrs);
  }
}

/* Copying */
static int igraphmodule_i_attribute_copy(igraph_t *to, const igraph_t *from,
  igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea) {
  igraphmodule_i_attribute_struct *fromattrs, *toattrs;
  PyObject *key, *value, *newval, *o=NULL;
  igraph_bool_t copy_attrs[3] = { ga, va, ea };
  int i, j;
  Py_ssize_t pos = 0;
 
  if (from->attr) {
    fromattrs=ATTR_STRUCT(from);
    /* what to do with the original value of toattrs? */
    toattrs=(igraphmodule_i_attribute_struct*)calloc(1, sizeof(igraphmodule_i_attribute_struct));
    if (!toattrs)
      IGRAPH_ERROR("not enough memory to allocate attribute hashes", IGRAPH_ENOMEM);
    if (igraphmodule_i_attribute_struct_init(toattrs)) {
      PyErr_Clear();
      free(toattrs);
      IGRAPH_ERROR("not enough memory to allocate attribute hashes", IGRAPH_ENOMEM);
    }
    to->attr=toattrs;

    for (i=0; i<3; i++) {
      if (!copy_attrs[i])
        continue;

      if (!PyDict_Check(fromattrs->attrs[i])) {
        toattrs->attrs[i]=fromattrs->attrs[i];
        Py_XINCREF(fromattrs->attrs[i]);
        continue;
      }
      
      pos = 0;
      while (PyDict_Next(fromattrs->attrs[i], &pos, &key, &value)) {
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
        PyDict_SetItem(toattrs->attrs[i], key, newval);
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
  igraph_attribute_record_t *attr_rec;
  igraph_bool_t *added_attrs=0;
  Py_ssize_t pos = 0;

  if (!graph->attr) return IGRAPH_SUCCESS;
  if (nv<0) return IGRAPH_SUCCESS;

  if (attr) {
    added_attrs = (igraph_bool_t*)calloc((size_t)igraph_vector_ptr_size(attr),
                                         sizeof(igraph_bool_t));
    if (!added_attrs)
      IGRAPH_ERROR("can't add vertex attributes", IGRAPH_ENOMEM);
    IGRAPH_FINALLY(free, added_attrs);
  }

  dict=ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_VERTEX];
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
        if (PyString_IsEqualToASCIIString(key, attr_rec->name)) {
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
        case IGRAPH_ATTRIBUTE_BOOLEAN:
          o=VECTOR(*(igraph_vector_bool_t*)attr_rec->value)[i] ? Py_True : Py_False;
          Py_INCREF(o);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string, numeric or Boolean)");
          o=0;
          break;
        }
        if (o) {
          if (PyList_Append(value, o) == -1)
            IGRAPH_ERROR("can't extend a vertex attribute hash member", IGRAPH_FAILURE);
          else Py_DECREF(o);
        }
      }

      /* Invalidate the vertex name index if needed */
      if (!strcmp(attr_rec->name, "name"))
        igraphmodule_i_attribute_struct_invalidate_vertex_name_index(ATTR_STRUCT(graph));
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
      attr_rec=(igraph_attribute_record_t*)VECTOR(*attr)[k];

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
        case IGRAPH_ATTRIBUTE_BOOLEAN:
          o=VECTOR(*(igraph_vector_bool_t*)attr_rec->value)[i] ? Py_True : Py_False;
          Py_INCREF(o);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string, numeric or Boolean)");
          o=0;
          break;
        }
        if (o) PyList_SET_ITEM(value, i+j, o);
      }

      /* Invalidate the vertex name index if needed */
      if (!strcmp(attr_rec->name, "name"))
        igraphmodule_i_attribute_struct_invalidate_vertex_name_index(ATTR_STRUCT(graph));

      PyDict_SetItemString(dict, attr_rec->name, value);
      Py_DECREF(value);   /* compensate for PyDict_SetItemString */
    }
    free(added_attrs);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return IGRAPH_SUCCESS;
}

/* Permuting vertices */
static int igraphmodule_i_attribute_permute_vertices(const igraph_t *graph,
    igraph_t *newgraph, const igraph_vector_t *idx) {
  long int n, i;
  PyObject *key, *value, *dict, *newdict, *newlist, *o;
  Py_ssize_t pos=0;
  
  dict=ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_VERTEX];
  if (!PyDict_Check(dict)) return 1;

  newdict=PyDict_New();
  if (!newdict) return 1;

  n=igraph_vector_size(idx);
  pos=0;

  while (PyDict_Next(dict, &pos, &key, &value)) {
    newlist=PyList_New(n);
    for (i=0; i<n; i++) {
      o = PyList_GetItem(value, (Py_ssize_t)VECTOR(*idx)[i]);
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

  dict = ATTR_STRUCT_DICT(newgraph)[ATTRHASH_IDX_VERTEX];
  ATTR_STRUCT_DICT(newgraph)[ATTRHASH_IDX_VERTEX]=newdict;
  Py_DECREF(dict);

  /* Invalidate the vertex name index */
  igraphmodule_i_attribute_struct_invalidate_vertex_name_index(ATTR_STRUCT(newgraph));

  return 0;
}

/* Adding edges */
static int igraphmodule_i_attribute_add_edges(igraph_t *graph, const igraph_vector_t *edges, igraph_vector_ptr_t *attr) {
  /* Extend the end of every value in the edge hash with ne pieces of None */
  PyObject *key, *value, *dict;
  Py_ssize_t pos=0;
  long int i, j, k, l, ne;
  igraph_bool_t *added_attrs=0;
  igraph_attribute_record_t *attr_rec;

  ne=igraph_vector_size(edges)/2;
  if (!graph->attr) return IGRAPH_SUCCESS;
  if (ne<0) return IGRAPH_SUCCESS;
  
  if (attr) {
    added_attrs = (igraph_bool_t*)calloc((size_t)igraph_vector_ptr_size(attr),
                                         sizeof(igraph_bool_t));
    if (!added_attrs)
      IGRAPH_ERROR("can't add vertex attributes", IGRAPH_ENOMEM);
    IGRAPH_FINALLY(free, added_attrs);
  }

  dict=ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
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
        if (PyString_IsEqualToASCIIString(key, attr_rec->name)) {
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
        case IGRAPH_ATTRIBUTE_BOOLEAN:
          o=VECTOR(*(igraph_vector_bool_t*)attr_rec->value)[i] ? Py_True : Py_False;
          Py_INCREF(o);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string, numeric or Boolean)");
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
      attr_rec=(igraph_attribute_record_t*)VECTOR(*attr)[k];

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
        case IGRAPH_ATTRIBUTE_BOOLEAN:
          o=VECTOR(*(igraph_vector_bool_t*)attr_rec->value)[i] ? Py_True : Py_False;
          Py_INCREF(o);
          break;
        default:
          IGRAPH_WARNING("unsupported attribute type (not string, numeric or Boolean)");
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

/* Deleting edges, currently unused */
/*
static void igraphmodule_i_attribute_delete_edges(igraph_t *graph, const igraph_vector_t *idx) {
  long int n, i, ndeleted=0;
  PyObject *key, *value, *dict, *o;
  Py_ssize_t pos=0;
  
  dict=ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  if (!PyDict_Check(dict)) return;

  n=igraph_vector_size(idx);
  for (i=0; i<n; i++) {
    if (!VECTOR(*idx)[i]) {
      ndeleted++;
      continue;
    }

    pos=0;
    while (PyDict_Next(dict, &pos, &key, &value)) {
      o=PyList_GetItem(value, i);
      if (!o) {
        PyErr_Clear();
        return;
      }
      Py_INCREF(o);
      PyList_SetItem(value, VECTOR(*idx)[i]-1, o);
    }
  }
  
  pos=0;
  while (PyDict_Next(dict, &pos, &key, &value)) {
    n=PySequence_Size(value);
    if (PySequence_DelSlice(value, n-ndeleted, n) == -1) return;
  }
  
  return;
}
*/

/* Permuting edges */
static int igraphmodule_i_attribute_permute_edges(const igraph_t *graph,
    igraph_t *newgraph, const igraph_vector_t *idx) { 
  long int n, i;
  PyObject *key, *value, *dict, *newdict, *newlist, *o;
  Py_ssize_t pos=0;

  dict=ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  if (!PyDict_Check(dict)) return 1;

  newdict=PyDict_New();
  if (!newdict) return 1;

  n=igraph_vector_size(idx);
  pos=0;

  while (PyDict_Next(dict, &pos, &key, &value)) {
    newlist=PyList_New(n);
    for (i=0; i<n; i++) {
      o=PyList_GetItem(value, (Py_ssize_t)VECTOR(*idx)[i]);
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

  dict = ATTR_STRUCT_DICT(newgraph)[ATTRHASH_IDX_EDGE];
  ATTR_STRUCT_DICT(newgraph)[ATTRHASH_IDX_EDGE]=newdict;
  Py_DECREF(dict);

  return 0;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged, the source
 * attribute values and a Python callable to be called for every merge,
 * returns a new list with the new attribute values. Each new attribute
 * is derived by calling func with the old attributes of the set of
 * merged vertices/edges as the first argument.
 */
static PyObject* igraphmodule_i_ac_func(PyObject* values,
    const igraph_vector_ptr_t *merges, PyObject* func) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *list, *item;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    long int j, n = igraph_vector_size(v);
    list = PyList_New(n);
    for (j = 0; j < n; j++) {
      item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[j]);
      Py_INCREF(item);
      PyList_SET_ITEM(list, j, item);   /* reference to item stolen */
    }
    item = PyObject_CallFunctionObjArgs(func, list, 0);
    Py_DECREF(list);
    if (item == 0) {
      Py_DECREF(res);
      return 0;
    }
    PyList_SET_ITEM(res, i, item);   /* reference to item stolen */
  }

  return res;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged, the source
 * attribute values and a name of a Python builtin function,
 * returns a new list with the new attribute values. Each new attribute
 * is derived by calling the given builtin function with the old
 * attributes of the set of merged vertices/edges as the first argument.
 */
static PyObject* igraphmodule_i_ac_builtin_func(PyObject* values,
    const igraph_vector_ptr_t *merges, const char* func_name) {
  static PyObject* builtin_module_dict = 0;
  PyObject* func = 0;

  if (builtin_module_dict == 0) {
#ifdef IGRAPH_PYTHON3
    PyObject* builtin_module = PyImport_ImportModule("builtins");
#else
    PyObject* builtin_module = PyImport_ImportModule("__builtin__");
#endif
    if (builtin_module == 0)
      return 0;
    builtin_module_dict = PyModule_GetDict(builtin_module);
    Py_DECREF(builtin_module);
    if (builtin_module_dict == 0)
      return 0;
  }

  func = PyDict_GetItemString(builtin_module_dict, func_name);
  if (func == 0) {
    PyErr_Format(PyExc_NameError, "no such builtin function; %s", func_name);
    return 0;
  }

  return igraphmodule_i_ac_func(values, merges, func);
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged and the source
 * attribute values, returns a new list with the new attribute values.
 * Each new attribute is derived from the sum of the attributes of
 * the merged vertices/edges.
 */
static PyObject* igraphmodule_i_ac_sum(PyObject* values,
    const igraph_vector_ptr_t *merges) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *item;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    igraph_real_t num = 0.0, sum = 0.0;
    long int j, n = igraph_vector_size(v);

    for (j = 0; j < n; j++) {
      item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[j]);
      if (igraphmodule_PyObject_to_real_t(item, &num)) {
        PyErr_SetString(PyExc_TypeError, "product can only be invoked on numeric attributes");
        Py_DECREF(res);
        return 0;
      }
      sum += num;
    }

    /* reference to new float stolen */
    PyList_SET_ITEM(res, i, PyFloat_FromDouble((double)sum));
  }

  return res;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged and the source
 * attribute values, returns a new list with the new attribute values.
 * Each new attribute is derived from the product of the attributes of
 * the merged vertices/edges.
 */
static PyObject* igraphmodule_i_ac_prod(PyObject* values,
    const igraph_vector_ptr_t *merges) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *item;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    igraph_real_t num = 1.0, prod = 1.0;
    long int j, n = igraph_vector_size(v);

    for (j = 0; j < n; j++) {
      item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[j]);
      if (igraphmodule_PyObject_to_real_t(item, &num)) {
        PyErr_SetString(PyExc_TypeError, "product can only be invoked on numeric attributes");
        Py_DECREF(res);
        return 0;
      }
      prod *= num;
    }

    /* reference to new float stolen */
    PyList_SET_ITEM(res, i, PyFloat_FromDouble((double)prod));
  }

  return res;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged and the source
 * attribute values, returns a new list with the new attribute values.
 * Each new attribute is derived from the first entry of the set of merged
 * vertices/edges.
 */
static PyObject* igraphmodule_i_ac_first(PyObject* values,
    const igraph_vector_ptr_t *merges) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *item;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    long int n = igraph_vector_size(v);

    if (n == 0)
      continue;

    item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[0]);
    Py_INCREF(item);
    PyList_SET_ITEM(res, i, item);   /* reference to item stolen */
  }

  return res;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged and the source
 * attribute values, returns a new list with the new attribute values.
 * Each new attribute is derived from a randomly selected entry of the set of
 * merged vertices/edges.
 */
static PyObject* igraphmodule_i_ac_random(PyObject* values,
    const igraph_vector_ptr_t *merges) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *item, *num;
  PyObject *random_module = PyImport_ImportModule("random");
  PyObject *random_func;

  if (random_module == 0)
    return 0;

  random_func = PyObject_GetAttrString(random_module, "random");
  Py_DECREF(random_module);

  if (random_func == 0)
    return 0;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    long int n = igraph_vector_size(v);

    if (n == 0)
      continue;

    num = PyObject_CallObject(random_func, 0);
    if (num == 0) {
      Py_DECREF(random_func);
      Py_DECREF(res);
      return 0;
    }

    item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[(long int)(n*PyFloat_AsDouble(num))]);
    Py_INCREF(item);
    Py_DECREF(num);
    PyList_SET_ITEM(res, i, item);   /* reference to item stolen */
  }

  Py_DECREF(random_func);

  return res;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged and the source
 * attribute values, returns a new list with the new attribute values.
 * Each new attribute is derived from the last entry of the set of merged
 * vertices/edges.
 */
static PyObject* igraphmodule_i_ac_last(PyObject* values,
    const igraph_vector_ptr_t *merges) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *item;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    long int n = igraph_vector_size(v);

    if (n == 0)
      continue;

    item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[n-1]);
    Py_INCREF(item);
    PyList_SET_ITEM(res, i, item);   /* reference to item stolen */
  }

  return res;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged and the source
 * attribute values, returns a new list with the new attribute values.
 * Each new attribute is derived from the mean of the attributes of
 * the merged vertices/edges.
 */
static PyObject* igraphmodule_i_ac_mean(PyObject* values,
    const igraph_vector_ptr_t *merges) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *item;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    igraph_real_t num = 0.0, mean = 0.0;
    long int j, n = igraph_vector_size(v);

    for (j = 0; j < n; ) {
      item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[j]);
      if (igraphmodule_PyObject_to_real_t(item, &num)) {
        PyErr_SetString(PyExc_TypeError, "mean can only be invoked on numeric attributes");
        Py_DECREF(res);
        return 0;
      }
      j++;
      num -= mean;
      mean += num / j;
    }

    /* reference to new float stolen */
    PyList_SET_ITEM(res, i, PyFloat_FromDouble((double)mean));
  }

  return res;
}

/* Auxiliary function for combining vertices/edges. Given a merge list
 * (which specifies the vertex/edge IDs that were merged and the source
 * attribute values, returns a new list with the new attribute values.
 * Each new attribute is derived from the median of the attributes of
 * the merged vertices/edges.
 */
static PyObject* igraphmodule_i_ac_median(PyObject* values,
    const igraph_vector_ptr_t *merges) {
  long int i, len = igraph_vector_ptr_size(merges);
  PyObject *res, *list, *item;

  res = PyList_New(len);
  for (i = 0; i < len; i++) {
    igraph_vector_t *v = (igraph_vector_t*)VECTOR(*merges)[i];
    long int j, n = igraph_vector_size(v);
    list = PyList_New(n);
    for (j = 0; j < n; j++) {
      item = PyList_GET_ITEM(values, (Py_ssize_t)VECTOR(*v)[j]);
      Py_INCREF(item);
      PyList_SET_ITEM(list, j, item);   /* reference to item stolen */
    }
    /* sort the list */
    if (PyList_Sort(list)) {
      Py_DECREF(list);
      Py_DECREF(res);
      return 0;
    }
    if (n % 2 == 1) {
      item = PyList_GET_ITEM(list, n / 2);
    } else {
      igraph_real_t num1, num2;
      item = PyList_GET_ITEM(list, n / 2 - 1);
      if (igraphmodule_PyObject_to_real_t(item, &num1)) {
        Py_DECREF(list);
        Py_DECREF(res);
        return 0;
      }
      item = PyList_GET_ITEM(list, n / 2);
      if (igraphmodule_PyObject_to_real_t(item, &num2)) {
        Py_DECREF(list);
        Py_DECREF(res);
        return 0;
      }
      item = PyFloat_FromDouble((num1 + num2) / 2);
    }
    /* reference to item stolen */
    PyList_SET_ITEM(res, i, item);
  }

  return res;
}

static void igraphmodule_i_free_attribute_combination_records(
    igraph_attribute_combination_record_t* records) {
  igraph_attribute_combination_record_t* ptr = records;
  while (ptr->name != 0) {
    free((char*)ptr->name);
    ptr++;
  }
  free(records);
}

/* Auxiliary function for the common parts of
 * igraphmodule_i_attribute_combine_vertices and
 * igraphmodule_i_attribute_combine_edges
 */
static int igraphmodule_i_attribute_combine_dicts(PyObject *dict,
    PyObject *newdict, const igraph_vector_ptr_t *merges,
    const igraph_attribute_combination_t *comb) {
  PyObject *key, *value;
  Py_ssize_t pos;
  igraph_attribute_combination_record_t* todo;
  Py_ssize_t i, n;
  if (!PyDict_Check(dict) || !PyDict_Check(newdict)) return 1;

  /* Allocate memory for the attribute_combination_records */
  n = PyDict_Size(dict);
  todo = (igraph_attribute_combination_record_t*)calloc(
    n+1, sizeof(igraph_attribute_combination_record_t)
  );
  if (todo == 0) {
    IGRAPH_ERROR("cannot allocate memory for attribute combination", IGRAPH_ENOMEM);
  }
  for (i = 0; i < n+1; i++)
    todo[i].name = 0;       /* sentinel elements */
  IGRAPH_FINALLY(igraphmodule_i_free_attribute_combination_records, todo);

  /* Collect what to do for each attribute in the source dict */
  pos = 0; i = 0;
  while (PyDict_Next(dict, &pos, &key, &value)) {
    todo[i].name = PyString_CopyAsString(key);
    if (todo[i].name == 0)
      IGRAPH_ERROR("PyString_CopyAsString failed", IGRAPH_FAILURE);
    igraph_attribute_combination_query(comb, todo[i].name,
        &todo[i].type, &todo[i].func);
    i++;
  }

  /* Combine the attributes. Here we make use of the fact that PyDict_Next
   * will iterate over the dict in the same order */
  pos = 0; i = 0;
  while (PyDict_Next(dict, &pos, &key, &value)) {
    PyObject *empty_str;
    PyObject *func;
    PyObject *newvalue;

    /* Safety check */
    if (!PyString_IsEqualToASCIIString(key, todo[i].name)) {
      IGRAPH_ERROR("PyDict_Next iteration order not consistent. "
          "This should never happen. Please report the bug to the igraph "
          "developers!", IGRAPH_FAILURE);
    }

    newvalue = 0;
    switch (todo[i].type) {
      case IGRAPH_ATTRIBUTE_COMBINE_DEFAULT:
      case IGRAPH_ATTRIBUTE_COMBINE_IGNORE:
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
        func = (PyObject*)todo[i].func;
        newvalue = igraphmodule_i_ac_func(value, merges, func);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_SUM:
        newvalue = igraphmodule_i_ac_sum(value, merges);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_PROD:
        newvalue = igraphmodule_i_ac_prod(value, merges);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_MIN:
        newvalue = igraphmodule_i_ac_builtin_func(value, merges, "min");
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_MAX:
        newvalue = igraphmodule_i_ac_builtin_func(value, merges, "max");
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
        newvalue = igraphmodule_i_ac_random(value, merges);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
        newvalue = igraphmodule_i_ac_first(value, merges);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_LAST:
        newvalue = igraphmodule_i_ac_last(value, merges);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
        newvalue = igraphmodule_i_ac_mean(value, merges);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
        newvalue = igraphmodule_i_ac_median(value, merges);
        break;

      case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
        empty_str = PyString_FromString("");
        func = PyObject_GetAttrString(empty_str, "join");
        newvalue = igraphmodule_i_ac_func(value, merges, func);
        Py_DECREF(func);
        Py_DECREF(empty_str);
        break;

      default:
        IGRAPH_ERROR("Unsupported combination type. "
            "This should never happen. Please report the bug to the igraph "
            "developers!", IGRAPH_FAILURE);
    }

    if (newvalue) {
      if (PyDict_SetItem(newdict, key, newvalue)) {
        Py_DECREF(newvalue);  /* PyDict_SetItem does not steal reference */
        IGRAPH_ERROR("PyDict_SetItem failed when combining attributes.", IGRAPH_FAILURE);
      }
      Py_DECREF(newvalue);    /* PyDict_SetItem does not steal reference */
    } else {
      /* We can arrive here for two reasons: first, if the attribute is to
       * be ignored explicitly; second, if there was an error. */
      if (PyErr_Occurred()) {
        IGRAPH_ERROR("Unexpected failure when combining attributes", IGRAPH_FAILURE);
      }
    }

    i++;
  }

  igraphmodule_i_free_attribute_combination_records(todo);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/* Combining vertices */
static int igraphmodule_i_attribute_combine_vertices(const igraph_t *graph,
    igraph_t *newgraph, const igraph_vector_ptr_t *merges,
    const igraph_attribute_combination_t *comb) {
  PyObject *dict, *newdict;
  int result;

  /* Get the attribute dicts */
  dict=ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_VERTEX];
  newdict=ATTR_STRUCT_DICT(newgraph)[ATTRHASH_IDX_VERTEX];

  /* Combine the attribute dicts */
  result = igraphmodule_i_attribute_combine_dicts(dict, newdict,
      merges, comb);

  /* Invalidate vertex name index */
  igraphmodule_i_attribute_struct_invalidate_vertex_name_index(ATTR_STRUCT(graph));

  return result;
}

/* Combining edges */
static int igraphmodule_i_attribute_combine_edges(const igraph_t *graph,
    igraph_t *newgraph, const igraph_vector_ptr_t *merges,
    const igraph_attribute_combination_t *comb) {
  PyObject *dict, *newdict;

  /* Get the attribute dicts */
  dict=ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  newdict=ATTR_STRUCT_DICT(newgraph)[ATTRHASH_IDX_EDGE];

  return igraphmodule_i_attribute_combine_dicts(dict, newdict,
      merges, comb);
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
  int retval;
  long int i, j, k, l, m;
  
  for (i=0; i<3; i++) {
    igraph_strvector_t *n = names[i];
    igraph_vector_t *t = types[i];
    PyObject *dict = ATTR_STRUCT_DICT(graph)[i];
    PyObject *keys;
    PyObject *values;
    PyObject *o=0;
    keys=PyDict_Keys(dict);
    if (!keys) IGRAPH_ERROR("Internal error in PyDict_Keys", IGRAPH_FAILURE);
 
    if (n) {
      retval = igraphmodule_PyList_to_strvector_t(keys, n);
      if (retval)
        return retval;
    }
    if (t) {
      k=PyList_Size(keys);
      igraph_vector_resize(t, k);
      for (j=0; j<k; j++) {
        int is_numeric = 1;
        int is_string = 1;
        int is_boolean = 1;
        values=PyDict_GetItem(dict, PyList_GetItem(keys, j));
        if (PyList_Check(values)) {
          m=PyList_Size(values);
          for (l=0; l<m && is_numeric; l++) {
            o=PyList_GetItem(values, l);
            if (o != Py_None && !PyNumber_Check(o))
              is_numeric=0;
          }
          for (l=0; l<m && is_string; l++) {
            o=PyList_GetItem(values, l);
            if (o != Py_None && !PyBaseString_Check(o))
              is_string=0;
          }
          for (l=0; l<m && is_boolean; l++) {
            o=PyList_GetItem(values, l);
            if (o != Py_None && o != Py_False && o != Py_True)
              is_boolean=0;
          }
        } else {
          if (values != Py_None && !PyNumber_Check(values))
            is_numeric=0;
          if (values != Py_None && !PyBaseString_Check(values))
            is_string=0;
          if (values != Py_None && values != Py_False && values != Py_True)
            is_boolean=0;
        }
        if (is_boolean)
          VECTOR(*t)[j] = IGRAPH_ATTRIBUTE_BOOLEAN;
        else if (is_numeric)
          VECTOR(*t)[j] = IGRAPH_ATTRIBUTE_NUMERIC;
        else if (is_string)
          VECTOR(*t)[j] = IGRAPH_ATTRIBUTE_STRING;
        else
          VECTOR(*t)[j] = IGRAPH_ATTRIBUTE_PY_OBJECT;
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
  switch (type) {
  case IGRAPH_ATTRIBUTE_GRAPH:
    return igraphmodule_has_graph_attribute(graph, name);
  case IGRAPH_ATTRIBUTE_VERTEX:
    return igraphmodule_has_vertex_attribute(graph, name);
  case IGRAPH_ATTRIBUTE_EDGE:
    return igraphmodule_has_edge_attribute(graph, name);
  default:
    return 0;
  }
}

/* Returns the type of a given attribute */
int igraphmodule_i_attribute_get_type(const igraph_t *graph,
				      igraph_attribute_type_t *type,
				      igraph_attribute_elemtype_t elemtype,
				      const char *name) {
  long int attrnum, i, j;
  int is_numeric, is_string, is_boolean;
  PyObject *o, *dict;

  switch (elemtype) {
  case IGRAPH_ATTRIBUTE_GRAPH:  attrnum=ATTRHASH_IDX_GRAPH;  break;
  case IGRAPH_ATTRIBUTE_VERTEX: attrnum=ATTRHASH_IDX_VERTEX; break;
  case IGRAPH_ATTRIBUTE_EDGE:   attrnum=ATTRHASH_IDX_EDGE;   break;
  default: IGRAPH_ERROR("No such attribute type", IGRAPH_EINVAL); break;
  }

  /* Get the attribute dict */
  dict = ATTR_STRUCT_DICT(graph)[attrnum];

  /* Check whether the attribute exists */
  o = PyDict_GetItemString(dict, name);
  if (o == 0) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  /* Basic type check */
  if (!PyList_Check(o)) IGRAPH_ERROR("attribute hash type mismatch", IGRAPH_EINVAL);
  j = PyList_Size(o);
  if (j == 0) {
    *type = IGRAPH_ATTRIBUTE_NUMERIC;
    return 0;
  }

  /* Go on with the checks */
  is_numeric = is_string = is_boolean = 1;
  if (attrnum>0) {

    for (i=0; i<j && is_numeric; i++) {
      PyObject *item = PyList_GET_ITEM(o, i);
      if (item != Py_None && !PyNumber_Check(item)) is_numeric=0;
    }
    for (i=0; i<j && is_string; i++) {
      PyObject *item = PyList_GET_ITEM(o, i);
      if (item != Py_None && !PyBaseString_Check(item))
        is_string=0;
    }
    for (i=0; i<j && is_boolean; i++) {
      PyObject *item = PyList_GET_ITEM(o, i);
      if (item != Py_None && item != Py_True && item != Py_False)
        is_boolean=0;
    }
  } else {
    if (o != Py_None && !PyNumber_Check(o))
      is_numeric=0;
    if (o != Py_None && !PyBaseString_Check(o))
      is_string=0;
    if (o != Py_None && o != Py_True && o != Py_False)
      is_boolean=0;
  }
  if (is_boolean)
    *type = IGRAPH_ATTRIBUTE_BOOLEAN;
  else if (is_numeric)
    *type = IGRAPH_ATTRIBUTE_NUMERIC;
  else if (is_string)
    *type = IGRAPH_ATTRIBUTE_STRING;
  else
    *type = IGRAPH_ATTRIBUTE_PY_OBJECT;
  return 0;
}

/* Getting Boolean graph attributes */
int igraphmodule_i_get_boolean_graph_attr(const igraph_t *graph,
					  const char *name, igraph_vector_bool_t *value) {
  PyObject *dict, *o;
  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_GRAPH];
  /* No error checking, if we get here, the type has already been checked by previous
     attribute handler calls... hopefully :) Same applies for the other handlers. */
  o = PyDict_GetItemString(dict, name);
  if (!o)
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  IGRAPH_CHECK(igraph_vector_bool_resize(value, 1));
  VECTOR(*value)[0] = PyObject_IsTrue(o);
  return 0;
}

/* Getting numeric graph attributes */
int igraphmodule_i_get_numeric_graph_attr(const igraph_t *graph,
					  const char *name, igraph_vector_t *value) {
  PyObject *dict, *o, *result;
  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_GRAPH];
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
  PyObject *dict, *o, *str = 0;
  const char* c_str;

  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_GRAPH];
  o = PyDict_GetItemString(dict, name);
  if (!o)
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  IGRAPH_CHECK(igraph_strvector_resize(value, 1));

#ifdef IGRAPH_PYTHON3
  /* For Python 3.x, we simply call PyObject_Str, which produces a
   * Unicode string, then encode it into UTF-8, except when we
   * already have a PyBytes object -- this is assumed to be in
   * UTF-8.
   */
  if (PyBytes_Check(o)) {
    str = o;
    Py_INCREF(str);
  } else {
    PyObject* unicode = PyObject_Str(o);
    if (unicode == 0)
      IGRAPH_ERROR("Internal error in PyObject_Str", IGRAPH_EINVAL);
    str = PyUnicode_AsEncodedString(unicode, "utf-8", "xmlcharrefreplace");
    Py_DECREF(unicode);
  }

#else
  /* For Python 2.x, we check whether we have received a string or a
   * Unicode string. Unicode strings are encoded into UTF-8, strings
   * are used intact.
   */
  if (PyUnicode_Check(o)) {
    str = PyUnicode_AsEncodedString(o, "utf-8", "xmlcharrefreplace");
  } else {
    str = PyObject_Str(o);
  }
#endif

  if (str == 0)
    IGRAPH_ERROR("Internal error in PyObject_Str", IGRAPH_EINVAL);
#ifdef IGRAPH_PYTHON3
  c_str = PyBytes_AS_STRING(str);
#else
  c_str = PyString_AS_STRING(str);
#endif
  IGRAPH_CHECK(igraph_strvector_set(value, 0, c_str));
  Py_XDECREF(str);

  return 0;
}

/* Getting numeric vertex attributes */
int igraphmodule_i_get_numeric_vertex_attr(const igraph_t *graph,
					   const char *name,
					   igraph_vs_t vs,
					   igraph_vector_t *value) {
  PyObject *dict, *list, *result, *o;
  igraph_vector_t newvalue;

  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_VERTEX];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_vs_is_all(&vs)) {
    if (igraphmodule_PyObject_float_to_vector_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_vector_update(value, &newvalue);
    igraph_vector_destroy(&newvalue);
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_VIT_SIZE(it)));
    while (!IGRAPH_VIT_END(it)) {
      o = PyList_GetItem(list, (Py_ssize_t)IGRAPH_VIT_GET(it));
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

  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_VERTEX];
  list = PyDict_GetItemString(dict, name);
  if (!list)
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

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
      int v=(int)IGRAPH_VIT_GET(it);
      char* str;

      result = PyList_GetItem(list, v);
      if (result == 0)
        IGRAPH_ERROR("null element in PyList", IGRAPH_EINVAL);

      str = PyObject_ConvertToCString(result);
      if (str == 0)
        IGRAPH_ERROR("error while calling PyObject_ConvertToCString", IGRAPH_EINVAL);

      /* Note: this is a bit inefficient here, PyObject_ConvertToCString
       * allocates a new string which could be copied into the string
       * vector straight away. Instead of that, the string vector makes
       * another copy. Probably the performance hit is not too severe.
       */
      igraph_strvector_set(value, i, str);
      free(str);

      IGRAPH_VIT_NEXT(it);
      i++;
    }
    igraph_vit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

/* Getting boolean vertex attributes */
int igraphmodule_i_get_boolean_vertex_attr(const igraph_t *graph,
					   const char *name,
					   igraph_vs_t vs,
					   igraph_vector_bool_t *value) {
  PyObject *dict, *list, *o;
  igraph_vector_bool_t newvalue;

  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_VERTEX];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_vs_is_all(&vs)) {
    if (igraphmodule_PyObject_to_vector_bool_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_vector_bool_update(value, &newvalue);
    igraph_vector_bool_destroy(&newvalue);
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_bool_resize(value, IGRAPH_VIT_SIZE(it)));
    while (!IGRAPH_VIT_END(it)) {
      o = PyList_GetItem(list, (Py_ssize_t)IGRAPH_VIT_GET(it));
      VECTOR(*value)[i] = PyObject_IsTrue(o);
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

  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_es_is_all(&es)) {
    if (igraphmodule_PyObject_float_to_vector_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_vector_update(value, &newvalue);
    igraph_vector_destroy(&newvalue);
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_EIT_SIZE(it)));
    while (!IGRAPH_EIT_END(it)) {
      o = PyList_GetItem(list, (Py_ssize_t)IGRAPH_EIT_GET(it));
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

  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
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
      char* str;

      result = PyList_GetItem(list, (Py_ssize_t)IGRAPH_EIT_GET(it));
      if (result == 0)
        IGRAPH_ERROR("null element in PyList", IGRAPH_EINVAL);

      str = PyObject_ConvertToCString(result);
      if (str == 0)
        IGRAPH_ERROR("error while calling PyObject_ConvertToCString", IGRAPH_EINVAL);

      /* Note: this is a bit inefficient here, PyObject_ConvertToCString
       * allocates a new string which could be copied into the string
       * vector straight away. Instead of that, the string vector makes
       * another copy. Probably the performance hit is not too severe.
       */
      igraph_strvector_set(value, i, str);
      free(str);

      IGRAPH_EIT_NEXT(it);
      i++;
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

/* Getting boolean edge attributes */
int igraphmodule_i_get_boolean_edge_attr(const igraph_t *graph,
					 const char *name,
					 igraph_es_t es,
					 igraph_vector_bool_t *value) {
  PyObject *dict, *list, *o;
  igraph_vector_bool_t newvalue;

  dict = ATTR_STRUCT_DICT(graph)[ATTRHASH_IDX_EDGE];
  list = PyDict_GetItemString(dict, name);
  if (!list) IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);

  if (igraph_es_is_all(&es)) {
    if (igraphmodule_PyObject_to_vector_bool_t(list, &newvalue))
      IGRAPH_ERROR("Internal error", IGRAPH_EINVAL);
    igraph_vector_bool_update(value, &newvalue);
    igraph_vector_bool_destroy(&newvalue);
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_bool_resize(value, IGRAPH_EIT_SIZE(it)));
    while (!IGRAPH_EIT_END(it)) {
      o = PyList_GetItem(list, (Py_ssize_t)IGRAPH_EIT_GET(it));
      VECTOR(*value)[i] = PyObject_IsTrue(o);
      IGRAPH_EIT_NEXT(it);
      i++;
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

static igraph_attribute_table_t igraphmodule_attribute_table = {
  igraphmodule_i_attribute_init,
  igraphmodule_i_attribute_destroy,
  igraphmodule_i_attribute_copy,
  igraphmodule_i_attribute_add_vertices,
  igraphmodule_i_attribute_permute_vertices,
  igraphmodule_i_attribute_combine_vertices,
  igraphmodule_i_attribute_add_edges,
  igraphmodule_i_attribute_permute_edges,
  igraphmodule_i_attribute_combine_edges,
  igraphmodule_i_attribute_get_info,
  igraphmodule_i_attribute_has_attr,
  igraphmodule_i_attribute_get_type,
  igraphmodule_i_get_numeric_graph_attr,
  igraphmodule_i_get_string_graph_attr,
  igraphmodule_i_get_boolean_graph_attr,
  igraphmodule_i_get_numeric_vertex_attr,
  igraphmodule_i_get_string_vertex_attr,
  igraphmodule_i_get_boolean_vertex_attr,
  igraphmodule_i_get_numeric_edge_attr,
  igraphmodule_i_get_string_edge_attr,
  igraphmodule_i_get_boolean_edge_attr,
};

void igraphmodule_initialize_attribute_handler(void) {
  igraph_i_set_attribute_table(&igraphmodule_attribute_table);
}

/**
 * Checks whether the given Python object can be a valid attribute name or not.
 * Returns 1 if the object could be used as an attribute name, 0 otherwise.
 * Also raises a suitable Python exception if needed.
 */
int igraphmodule_attribute_name_check(PyObject* obj) {
  PyObject* type_str;

  if (obj != 0 && PyBaseString_Check(obj))
    return 1;

  type_str = obj ? PyObject_Str((PyObject*)obj->ob_type) : 0;
  if (type_str != 0) {
    PyErr_Format(PyExc_TypeError, "igraph supports string attribute names only, got %s",
        PyString_AS_STRING(type_str));
    Py_DECREF(type_str);
  } else {
    PyErr_Format(PyExc_TypeError, "igraph supports string attribute names only");
  }

  return 0;
}

