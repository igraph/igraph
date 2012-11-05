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

/************************ Miscellaneous functions *************************/

#include <Python.h>
#include <limits.h>
#include "attributes.h"
#include "graphobject.h"
#include "vertexseqobject.h"
#include "vertexobject.h"
#include "edgeseqobject.h"
#include "edgeobject.h"
#include "convert.h"
#include "error.h"
#include "memory.h"
#include "py2compat.h"

/**
 * \brief Converts a Python integer to a C int
 *
 * This is similar to PyInt_AsLong, but it checks for overflow first and throws
 * an exception if necessary.
 *
 * Returns -1 if there was an error, 0 otherwise.
 */
int PyInt_AsInt(PyObject* obj, int* result) {
  long dummy = PyInt_AsLong(obj);
  if (dummy < INT_MIN) {
    PyErr_SetString(PyExc_OverflowError,
        "integer too small for conversion to C int");
    return -1;
  }
  if (dummy > INT_MAX) {
    PyErr_SetString(PyExc_OverflowError,
        "integer too large for conversion to C int");
    return -1;
  }
  *result = (int)dummy;
  return 0;
}

/**
 * \brief Converts a Python long to a C int
 *
 * This is similar to PyLong_AsLong, but it checks for overflow first and
 * throws an exception if necessary.
 *
 * Returns -1 if there was an error, 0 otherwise.
 */
int PyLong_AsInt(PyObject* obj, int* result) {
  long dummy = PyLong_AsLong(obj);
  if (dummy < INT_MIN) {
    PyErr_SetString(PyExc_OverflowError,
        "long integer too small for conversion to C int");
    return -1;
  }
  if (dummy > INT_MAX) {
    PyErr_SetString(PyExc_OverflowError,
        "long integer too large for conversion to C int");
    return -1;
  }
  *result = (int)dummy;
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to a corresponding igraph enum.
 *
 * The numeric value is returned as an integer that must be converted
 * explicitly to the corresponding igraph enum type. This is to allow one
 * to use the same common conversion routine for multiple enum types.
 *
 * \param o a Python object to be converted
 * \param translation the translation table between strings and the
 *   enum values. Strings are treated as case-insensitive, but it is
 *   assumed that the translation table keys are lowercase. The last
 *   entry of the table must contain NULL values.
 * \param result the result is returned here. The default value must be
 *   passed in before calling this function, since this value is
 *   returned untouched if the given Python object is Py_None.
 * \return 0 if everything is OK, 1 otherwise. An appropriate exception
 *   is raised in this case.
 */
int igraphmodule_PyObject_to_enum(PyObject *o,
  igraphmodule_enum_translation_table_entry_t* table,
  int *result) {
    char *s, *s2;
    int i, best, best_result, best_unique;
    
    if (o == 0 || o == Py_None)
      return 0;
    if (PyInt_Check(o))
      return PyInt_AsInt(o, result);
    if (PyLong_Check(o))
      return PyLong_AsInt(o, result);
    s = PyString_CopyAsString(o);
    if (s == 0) {
        PyErr_SetString(PyExc_TypeError, "int, long or string expected");
        return -1;
    }
    /* Convert string to lowercase */
    for (s2=s; *s2; s2++)
      *s2 = tolower(*s2);
    best = 0; best_unique = 0; best_result = -1;
    /* Search for matches */
    while (table->name != 0) {
        if (strcmp(s, table->name) == 0) {
          *result = table->value;
          free(s);
          return 0;
        }
        for (i=0; s[i] == table->name[i]; i++);
        if (i > best) {
            best = i; best_unique = 1; best_result = table->value;
        } else if (i == best) best_unique = 0;
        table++;
    }
    free(s);
    if (best_unique) { *result = best_result; return 0; }
    PyErr_SetObject(PyExc_ValueError, o);
    return -1;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_neimode_t
 */
int igraphmodule_PyObject_to_neimode_t(PyObject *o,
  igraph_neimode_t *result) {
  static igraphmodule_enum_translation_table_entry_t neimode_tt[] = {
        {"in", IGRAPH_IN},
        {"out", IGRAPH_OUT},
        {"all", IGRAPH_ALL},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, neimode_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_add_weights_t
 */
int igraphmodule_PyObject_to_add_weights_t(PyObject *o,
  igraph_add_weights_t *result) {
  static igraphmodule_enum_translation_table_entry_t add_weights_tt[] = {
        {"true", IGRAPH_ADD_WEIGHTS_YES},
        {"yes", IGRAPH_ADD_WEIGHTS_YES},
        {"false", IGRAPH_ADD_WEIGHTS_NO},
        {"no", IGRAPH_ADD_WEIGHTS_NO},
        {"auto", IGRAPH_ADD_WEIGHTS_IF_PRESENT},
        {"if_present", IGRAPH_ADD_WEIGHTS_IF_PRESENT},
        {0,0}
    };

  if (o == Py_True) {
    *result = IGRAPH_ADD_WEIGHTS_YES;
    return 0;
  }

  if (o == Py_False) {
    *result = IGRAPH_ADD_WEIGHTS_NO;
    return 0;
  }

  return igraphmodule_PyObject_to_enum(o, add_weights_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_adjacency_t
 */
int igraphmodule_PyObject_to_adjacency_t(PyObject *o,
  igraph_adjacency_t *result) {
  static igraphmodule_enum_translation_table_entry_t adjacency_tt[] = {
        {"directed", IGRAPH_ADJ_DIRECTED},
        {"undirected", IGRAPH_ADJ_UNDIRECTED},
        {"upper", IGRAPH_ADJ_UPPER},
        {"lower", IGRAPH_ADJ_LOWER},
        {"minimum", IGRAPH_ADJ_MIN},
        {"maximum", IGRAPH_ADJ_MAX},
        {"plus", IGRAPH_ADJ_PLUS},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, adjacency_tt, (int*)result);
}

int igraphmodule_PyObject_to_attribute_combination_type_t(PyObject* o,
    igraph_attribute_combination_type_t *result) {
  static igraphmodule_enum_translation_table_entry_t attribute_combination_type_tt[] = {
        {"ignore", IGRAPH_ATTRIBUTE_COMBINE_IGNORE},
        {"sum", IGRAPH_ATTRIBUTE_COMBINE_SUM},
        {"product", IGRAPH_ATTRIBUTE_COMBINE_PROD},
        {"min", IGRAPH_ATTRIBUTE_COMBINE_MIN},
        {"max", IGRAPH_ATTRIBUTE_COMBINE_MAX},
        {"random", IGRAPH_ATTRIBUTE_COMBINE_RANDOM},
        {"first", IGRAPH_ATTRIBUTE_COMBINE_FIRST},
        {"last", IGRAPH_ATTRIBUTE_COMBINE_LAST},
        {"mean", IGRAPH_ATTRIBUTE_COMBINE_MEAN},
        {"median", IGRAPH_ATTRIBUTE_COMBINE_MEDIAN},
        {"concatenate", IGRAPH_ATTRIBUTE_COMBINE_CONCAT},
        {0, 0}
  };

  if (o == Py_None) {
    *result = IGRAPH_ATTRIBUTE_COMBINE_IGNORE;
    return 0;
  }

  if (PyCallable_Check(o)) {
    *result = IGRAPH_ATTRIBUTE_COMBINE_FUNCTION;
    return 0;
  }

  return igraphmodule_PyObject_to_enum(o, attribute_combination_type_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_barabasi_algorithm_t
 */
int igraphmodule_PyObject_to_barabasi_algorithm_t(PyObject *o,
  igraph_barabasi_algorithm_t *result) {
  static igraphmodule_enum_translation_table_entry_t barabasi_algorithm_tt[] = {
        {"bag", IGRAPH_BARABASI_BAG},
        {"psumtree", IGRAPH_BARABASI_PSUMTREE},
        {"psumtree_multiple", IGRAPH_BARABASI_PSUMTREE_MULTIPLE},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, barabasi_algorithm_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_connectedness_t
 */
int igraphmodule_PyObject_to_connectedness_t(PyObject *o,
  igraph_connectedness_t *result) {
  static igraphmodule_enum_translation_table_entry_t connectedness_tt[] = {
        {"weak", IGRAPH_WEAK},
        {"string", IGRAPH_STRONG},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, connectedness_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_vconn_nei_t
 */
int igraphmodule_PyObject_to_vconn_nei_t(PyObject *o,
  igraph_vconn_nei_t *result) {
  static igraphmodule_enum_translation_table_entry_t vconn_nei_tt[] = {
        {"error", IGRAPH_VCONN_NEI_ERROR},
        {"negative", IGRAPH_VCONN_NEI_NEGATIVE},
        {"number_of_nodes", IGRAPH_VCONN_NEI_NUMBER_OF_NODES},
        {"nodes", IGRAPH_VCONN_NEI_NUMBER_OF_NODES},
        {"ignore", IGRAPH_VCONN_NEI_IGNORE},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, vconn_nei_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_bliss_sh_t
 */
int igraphmodule_PyObject_to_bliss_sh_t(PyObject *o,
  igraph_bliss_sh_t *result) {
  static igraphmodule_enum_translation_table_entry_t bliss_sh_tt[] = {
        {"f", IGRAPH_BLISS_F},
        {"fl", IGRAPH_BLISS_FL},
        {"fs", IGRAPH_BLISS_FS},
        {"fm", IGRAPH_BLISS_FM},
        {"flm", IGRAPH_BLISS_FLM},
        {"fsm", IGRAPH_BLISS_FSM},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, bliss_sh_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_community_comparison_t
 */
int igraphmodule_PyObject_to_community_comparison_t(PyObject *o,
                  igraph_community_comparison_t *result) {
  static igraphmodule_enum_translation_table_entry_t commcmp_tt[] = {
        {"vi", IGRAPH_COMMCMP_VI},
        {"meila", IGRAPH_COMMCMP_VI},
        {"nmi", IGRAPH_COMMCMP_NMI},
        {"danon", IGRAPH_COMMCMP_NMI},
        {"split-join", IGRAPH_COMMCMP_SPLIT_JOIN},
        {"rand", IGRAPH_COMMCMP_RAND},
        {"adjusted_rand", IGRAPH_COMMCMP_ADJUSTED_RAND},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, commcmp_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_degseq_t
 */
int igraphmodule_PyObject_to_degseq_t(PyObject *o,
  igraph_degseq_t *result) {
  static igraphmodule_enum_translation_table_entry_t degseq_tt[] = {
        {"simple", IGRAPH_DEGSEQ_SIMPLE},
        {"no_multiple", IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE},
        {"vl", IGRAPH_DEGSEQ_VL},
        {"viger-latapy", IGRAPH_DEGSEQ_VL},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, degseq_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_fas_algorithm_t
 */
int igraphmodule_PyObject_to_fas_algorithm_t(PyObject *o,
  igraph_fas_algorithm_t *result) {
  static igraphmodule_enum_translation_table_entry_t fas_algorithm_tt[] = {
        {"approx_eades", IGRAPH_FAS_APPROX_EADES},
        {"eades", IGRAPH_FAS_APPROX_EADES},
        {"exact", IGRAPH_FAS_EXACT_IP},
        {"exact_ip", IGRAPH_FAS_EXACT_IP},
        {"ip", IGRAPH_FAS_EXACT_IP},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, fas_algorithm_tt, (int*)result);
}

/**
 * \brief Converts a Python object to an igraph \c igraph_reciprocity_t
 */
int igraphmodule_PyObject_to_reciprocity_t(PyObject *o, igraph_reciprocity_t *result) {
  static igraphmodule_enum_translation_table_entry_t reciprocity_tt[] = {
    {"default", IGRAPH_RECIPROCITY_DEFAULT},
    {"ratio", IGRAPH_RECIPROCITY_RATIO},
    {0,0}
  };

  return igraphmodule_PyObject_to_enum(o, reciprocity_tt, (int*)result);
}

/**
 * \brief Converts a Python object to an igraph \c igraph_rewiring_t
 */
int igraphmodule_PyObject_to_rewiring_t(PyObject *o, igraph_rewiring_t *result) {
  static igraphmodule_enum_translation_table_entry_t rewiring_tt[] = {
    {"simple", IGRAPH_REWIRING_SIMPLE},
    {"simple_loops", IGRAPH_REWIRING_SIMPLE_LOOPS},
    {"loops", IGRAPH_REWIRING_SIMPLE_LOOPS},
    {0,0}
  };

  return igraphmodule_PyObject_to_enum(o, rewiring_tt, (int*)result);
}

/**
 * \brief Converts a Python object to an igraph \c igraph_spinglass_implementation_t
 */
int igraphmodule_PyObject_to_spinglass_implementation_t(PyObject *o, igraph_spinglass_implementation_t *result) {
  static igraphmodule_enum_translation_table_entry_t spinglass_implementation_tt[] = {
    {"original", IGRAPH_SPINCOMM_IMP_ORIG},
    {"negative", IGRAPH_SPINCOMM_IMP_NEG},
    {0,0}
  };

  return igraphmodule_PyObject_to_enum(o, spinglass_implementation_tt, (int*)result);
}

/**
 * \brief Converts a Python object to an igraph \c igraph_spincomm_update_t
 */
int igraphmodule_PyObject_to_spincomm_update_t(PyObject *o, igraph_spincomm_update_t *result) {
  static igraphmodule_enum_translation_table_entry_t spincomm_update_tt[] = {
    {"simple", IGRAPH_SPINCOMM_UPDATE_SIMPLE},
    {"config", IGRAPH_SPINCOMM_UPDATE_CONFIG},
    {0,0}
  };

  return igraphmodule_PyObject_to_enum(o, spincomm_update_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_star_mode_t
 */
int igraphmodule_PyObject_to_star_mode_t(PyObject *o,
  igraph_star_mode_t *result) {
  static igraphmodule_enum_translation_table_entry_t star_mode_tt[] = {
        {"in", IGRAPH_STAR_IN},
        {"out", IGRAPH_STAR_OUT},
        {"mutual", IGRAPH_STAR_MUTUAL},
        {"undirected", IGRAPH_STAR_UNDIRECTED},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, star_mode_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_subgraph_implementation_t
 */
int igraphmodule_PyObject_to_subgraph_implementation_t(PyObject *o,
  igraph_subgraph_implementation_t *result) {
  static igraphmodule_enum_translation_table_entry_t subgraph_impl_tt[] = {
        {"auto", IGRAPH_SUBGRAPH_AUTO},
        {"copy_and_delete", IGRAPH_SUBGRAPH_COPY_AND_DELETE},
        {"old", IGRAPH_SUBGRAPH_COPY_AND_DELETE},
        {"create_from_scratch", IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH},
        {"new", IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, subgraph_impl_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_to_undirected_t
 */
int igraphmodule_PyObject_to_to_undirected_t(PyObject *o,
  igraph_to_undirected_t *result) {
  static igraphmodule_enum_translation_table_entry_t to_undirected_tt[] = {
        {"each",     IGRAPH_TO_UNDIRECTED_EACH},
        {"collapse", IGRAPH_TO_UNDIRECTED_COLLAPSE},
        {"mutual",   IGRAPH_TO_UNDIRECTED_MUTUAL},
        {0,0}
  };

  if (o == Py_True) {
    *result = IGRAPH_TO_UNDIRECTED_COLLAPSE;
    return 0;
  } else if (o == Py_False) {
    *result = IGRAPH_TO_UNDIRECTED_EACH;
    return 0;
  }

  return igraphmodule_PyObject_to_enum(o, to_undirected_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an \c igraph_transitivity_mode_t
 */
int igraphmodule_PyObject_to_transitivity_mode_t(PyObject *o,
  igraph_transitivity_mode_t *result) {
  static igraphmodule_enum_translation_table_entry_t transitivity_mode_tt[] = {
        {"zero", IGRAPH_TRANSITIVITY_ZERO},
        {"0", IGRAPH_TRANSITIVITY_ZERO},
        {"nan", IGRAPH_TRANSITIVITY_NAN},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, transitivity_mode_tt, (int*)result);
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_tree_mode_t
 */
int igraphmodule_PyObject_to_tree_mode_t(PyObject *o,
  igraph_tree_mode_t *result) {
  static igraphmodule_enum_translation_table_entry_t tree_mode_tt[] = {
        {"in", IGRAPH_TREE_IN},
        {"out", IGRAPH_TREE_OUT},
        {"all", IGRAPH_TREE_UNDIRECTED},
        {"undirected", IGRAPH_TREE_UNDIRECTED},
        {"tree_in", IGRAPH_TREE_IN},
        {"tree_out", IGRAPH_TREE_OUT},
        {"tree_all", IGRAPH_TREE_UNDIRECTED},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, tree_mode_tt, (int*)result);
}

/**
 * \brief Extracts a pointer to the internal \c igraph_t from a graph object
 *
 * Raises suitable Python exceptions when needed.
 *
 * \param object the Python object to be converted. If it is Py_None, the
 *               result pointer is untouched (so it should be null by default).
 * \param result the pointer is stored here
 *
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_igraph_t(PyObject *o, igraph_t **result) {
  if (o == Py_None)
    return 0;

  if (!PyObject_TypeCheck(o, &igraphmodule_GraphType)) {
    PyErr_Format(PyExc_TypeError,
        "expected graph object, got %s", o->ob_type->tp_name);
    return 1;
  }

  *result = &((igraphmodule_GraphObject*)o)->g;
  return 0;
}

/**
 * \brief Converts a Python object to an igraph \c igraph_integer_t
 *
 * Raises suitable Python exceptions when needed.
 *
 * \param object the Python object to be converted
 * \param v the result is returned here
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_integer_t(PyObject *object, igraph_integer_t *v) {
  int retval, num;

  if (object == NULL) {
  } else if (PyLong_Check(object)) {
    retval = PyLong_AsInt(object, &num);
    if (retval)
      return retval;
    *v = num;
    return 0;
#ifdef IGRAPH_PYTHON3
  } else if (PyNumber_Check(object)) {
    PyObject *i = PyNumber_Int(object);
    if (i == NULL)
      return 1;
    retval = PyInt_AsInt(i, &num);
    Py_DECREF(i);
    if (retval)
      return retval;
    *v = num;
    return 0;
  }
#else
  } else if (PyInt_Check(object)) {
    retval = PyInt_AsInt(object, &num);
    if (retval)
      return retval;
    *v = num;
    return 0;
  } else if (PyNumber_Check(object)) {
    PyObject *i = PyNumber_Int(object);
    if (i == NULL)
      return 1;
    retval = PyInt_AsInt(i, &num);
    Py_DECREF(i);
    if (retval)
      return retval;
    *v = num;
    return 0;
  }
#endif
  PyErr_BadArgument();
  return 1;
}

/**
 * \brief Converts a Python object to an igraph \c igraph_real_t
 *
 * Raises suitable Python exceptions when needed.
 *
 * \param object the Python object to be converted
 * \param v the result is returned here
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_real_t(PyObject *object, igraph_real_t *v) {
  if (object == NULL) {
  } else if (PyLong_Check(object)) {
    double d = PyLong_AsDouble(object);
    *v=(igraph_real_t)d;
    return 0;
#ifndef IGRAPH_PYTHON3
  } else if (PyInt_Check(object)) {
    long l = PyInt_AS_LONG((PyIntObject*)object);
    *v=(igraph_real_t)l;
    return 0;
#endif
  } else if (PyFloat_Check(object)) {
    double d = PyFloat_AS_DOUBLE((PyFloatObject*)object);
    *v=(igraph_real_t)d;
    return 0;
  } else if (PyNumber_Check(object)) {
    PyObject *i = PyNumber_Float(object);
    double d;
    if (i == NULL) return 1;
    d = PyFloat_AS_DOUBLE((PyFloatObject*)i);
    Py_DECREF(i);
    *v = (igraph_real_t)d;
    return 0;
  }
  PyErr_BadArgument();
  return 1;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python object to an igraph \c igraph_vector_t
 * The incoming \c igraph_vector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * 
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_t containing the result
 * \param need_non_negative if true, checks whether all elements are non-negative
 * \param pairs if true, assumes that every list element is a pair of integers
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vector_t(PyObject *list, igraph_vector_t *v, igraph_bool_t need_non_negative, igraph_bool_t pairs) {
  PyObject *item, *i1, *i2;
  Py_ssize_t i, j, k;
  int ok;
  long int idx=0, idx2=0;

  if (PyBaseString_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    if (pairs)
      PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing integer pairs");
    else
      PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing integers");
    return 1;
  }

  if (PySequence_Check(list)) {
    ok=1;
    
    if (pairs && PyTuple_Check(list) && PyTuple_Size(list)==2 &&
        PyInt_Check(PyTuple_GetItem(list,0)) &&
        PyInt_Check(PyTuple_GetItem(list,1))) {
      /* a pair was given instead of a list */
      /* Assume that the user meant a list consisting of this single pair */
      i1=i2=NULL;
      i1=PyTuple_GetItem(list, 0);
      if (i1) i2=PyTuple_GetItem(list, 1);
      if (i1 && i2) {
        idx=PyInt_AsLong(i1); idx2=PyInt_AsLong(i2);
        if (need_non_negative && (idx<0 || idx2<0)) ok=0;
      } else ok=0;
      if (ok) {
        igraph_vector_init(v, 2);
        VECTOR(*v)[0]=(igraph_real_t)idx;
        VECTOR(*v)[1]=(igraph_real_t)idx2;
        return 0;
      } else if (!ok && need_non_negative) {
        PyErr_SetString(PyExc_TypeError, "sequence elements must be non-negative integers");
        return 1;
      } else {
        PyErr_SetString(PyExc_TypeError, "sequence elements must be integers");
        return 1;
      }
    }
  } else if (!pairs && PyInt_Check(list)) {
    /* a single integer was given instead of a list */
    /* Let's assume that the user meant a list consisting of this single item */
    igraph_vector_init(v, 1);
    VECTOR(*v)[0]=(igraph_real_t)PyInt_AsLong(list);
    return 0;
  } else if (!pairs && PyLong_Check(list)) {
    /* a single long was given instead of a list */
    /* Let's assume that the user meant a list consisting of this single item */
    igraph_vector_init(v, 1);
    VECTOR(*v)[0]=(igraph_real_t)PyLong_AsDouble(list);
    return 0;
  }

  if (!PySequence_Check(list)) {
    /* try to use an iterator */
    PyObject *it = PyObject_GetIter(list);
    if (it) {
      PyObject *item;
      igraph_vector_init(v, 0);
      while ((item = PyIter_Next(it)) != 0) {
        ok = 1;
        if (pairs) {
          if (!PySequence_Check(item) || PySequence_Size(item) != 2) {
            PyErr_SetString(PyExc_TypeError, "iterable must return pairs of integers");
            ok=0;
          } else {
            i1=i2=NULL;
            i1=PySequence_GetItem(item, 0);
            if (i1) i2=PySequence_GetItem(item, 1);
            if (i1 && i2 && PyInt_Check(i1) && PyInt_Check(i2)) {
              idx=PyInt_AsLong(i1); idx2=PyInt_AsLong(i2);
              if (need_non_negative && (idx<0 || idx2<0)) {
                PyErr_SetString(PyExc_ValueError, "iterable must return non-negative integer pairs");
                ok=0;
              }
            } else {
              PyErr_SetString(PyExc_ValueError, "iterable must return pairs of integers");
              ok=0;
            }
            Py_XDECREF(i1); Py_XDECREF(i2); /* PySeq_GetItem returned new ref */
          }
        } else {
          if (!PyInt_Check(item)) {
            PyErr_SetString(PyExc_ValueError, "iterable must return integers");
            ok=0;
          } else {
            idx=PyInt_AsLong(item);
            if (need_non_negative && idx<0) {
              PyErr_SetString(PyExc_ValueError, "iterable must return non-negative integers");
              ok=0;
            }
          }
        }
       
        if (ok == 0) {
          igraph_vector_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        if (igraph_vector_push_back(v, idx)) {
          igraphmodule_handle_igraph_error();
          igraph_vector_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        if (pairs) {
          if (igraph_vector_push_back(v, idx2)) {
            igraphmodule_handle_igraph_error();
            igraph_vector_destroy(v);
            Py_DECREF(item);
            Py_DECREF(it);
            return 1;
          }
        }
        Py_DECREF(item);
      }
      Py_DECREF(it);
      return 0;
    } else {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      return 1;
    }
    return 0;
  }

  j=PySequence_Size(list);
  if (pairs)
    igraph_vector_init(v, 2*j);
  else
    igraph_vector_init(v, j);
  for (i=0, k=0; i<j; i++) {
    item=PySequence_GetItem(list, i);
    if (item) {
      ok=1;
      if (pairs) {
        if (PySequence_Check(item) && PySequence_Size(item)==2) {
          i1=NULL; i2=NULL;
          i1=PySequence_GetItem(item, 0);
          if (i1) i2=PySequence_GetItem(item, 1);
          if (i1 && i2) {
            if (PyInt_Check(i1) && PyInt_Check(i2)) {
              idx=PyInt_AsLong(i1);
              idx2=PyInt_AsLong(i2);
              if (need_non_negative && (idx<0 || idx2<0)) {
                PyErr_SetString(PyExc_TypeError, "sequence elements must be non-negative integer pairs");
                ok=0;
              }
            } else {
              PyErr_SetString(PyExc_TypeError, "sequence elements must be integer pairs");
              ok=0;
            }
            Py_XDECREF(i1); Py_XDECREF(i2);
          } else {
            /* this should not happen, but we return anyway.
             * An IndexError exception was set by PySequence_GetItem
             * at this point */
            igraph_vector_destroy(v);
            Py_XDECREF(i1); Py_XDECREF(i2);
            Py_XDECREF(item);
            return 1;
          }
        } else {
          PyErr_SetString(PyExc_TypeError, "sequence elements must be integer pairs");
          ok=0;
        }
      } else {
        if (PyInt_Check(item)) {
          idx=PyInt_AsLong(item);
          if (need_non_negative && idx<0) {
            PyErr_SetString(PyExc_TypeError, "sequence elements must be non-negative integers");
            ok=0;
          }
        } else {
          PyErr_SetString(PyExc_TypeError, "sequence elements must be integers");
          ok=0;
        }
        Py_XDECREF(item);
      }
         
      if (!ok) {
        igraph_vector_destroy(v);
        return 1;
      }
          
      /* add idx into index vector */
      VECTOR(*v)[k]=(igraph_real_t)idx;
      k++;
      if (pairs) {
        VECTOR(*v)[k]=(igraph_real_t)idx2;
        k++;
      }
    } else {
      /* this should not happen, but we return anyway.
       * an IndexError exception was set by PyList_GetItem
       * at this point */
      igraph_vector_destroy(v);
      return 1;
    }
  }
   
  return 0;
}


/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of floats to an igraph \c igraph_vector_t
 * The incoming \c igraph_vector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * 
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_float_to_vector_t(PyObject *list, igraph_vector_t *v) {
  PyObject *item;
  igraph_real_t value=0;
  Py_ssize_t i, j, k;
  int ok;

  if (PyBaseString_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing floats");
    return 1;
  }

  if (!PySequence_Check(list)) {
    /* try to use an iterator */
    PyObject *it = PyObject_GetIter(list);
    if (it) {
      PyObject *item;
      igraph_vector_init(v, 0);
      while ((item = PyIter_Next(it)) != 0) {
        ok = 1;
        if (!PyNumber_Check(item)) {
          PyErr_SetString(PyExc_TypeError, "iterable must return numbers");
          ok=0;
        } else {
          PyObject *item2 = PyNumber_Float(item);
          if (item2 == 0) {
            PyErr_SetString(PyExc_TypeError, "can't convert a list item to float");
            ok = 0;
          } else {
            value=(igraph_real_t)PyFloat_AsDouble(item);
            Py_DECREF(item2);
          }
        }
       
        if (ok == 0) {
          igraph_vector_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        if (igraph_vector_push_back(v, value)) {
          igraphmodule_handle_igraph_error();
          igraph_vector_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        Py_DECREF(item);
      }
      Py_DECREF(it);
      return 0;
    } else {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      return 1;
    }
    return 0;
  }

  j=PySequence_Size(list);
  igraph_vector_init(v, j);
  for (i=0, k=0; i<j; i++) {
    item=PySequence_GetItem(list, i);
    if (item) {
      ok=1;
      if (!PyNumber_Check(item)) {
        PyErr_SetString(PyExc_TypeError, "sequence elements must be integers");
        ok=0;
      } else {
        PyObject *item2 = PyNumber_Float(item);
        if (item2 == 0) {
          PyErr_SetString(PyExc_TypeError, "can't convert sequence element to float");
          ok=0;
        } else {
          value=(igraph_real_t)PyFloat_AsDouble(item2);
          Py_DECREF(item2);
        }
      }
      Py_XDECREF(item);
      if (!ok) {
        igraph_vector_destroy(v);
        return 1;
      }
      VECTOR(*v)[k]=value;
      k++;
    } else {
      /* this should not happen, but we return anyway.
       * an IndexError exception was set by PyList_GetItem
       * at this point */
      igraph_vector_destroy(v);
      return 1;
    }
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of ints to an igraph \c igraph_vector_int_t
 * The incoming \c igraph_vector_int_t should be uninitialized.
 * Raises suitable Python exceptions when needed.
 *
 * This function is almost identical to
 * \ref igraphmodule_PyObject_to_vector_t . Make sure you fix bugs
 * in both cases (if any).
 *
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_int_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vector_int_t(PyObject *list, igraph_vector_int_t *v) {
  PyObject *item;
  int value=0;
  Py_ssize_t i, j, k;
  int ok, retval;

  if (PyBaseString_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing integers");
    return 1;
  }

  if (!PySequence_Check(list)) {
    /* try to use an iterator */
    PyObject *it = PyObject_GetIter(list);
    if (it) {
      PyObject *item;
      igraph_vector_int_init(v, 0);
      while ((item = PyIter_Next(it)) != 0) {
        ok = 1;
        if (!PyNumber_Check(item)) {
          PyErr_SetString(PyExc_TypeError, "iterable must return numbers");
          ok=0;
        } else {
          PyObject *item2 = PyNumber_Int(item);
          if (item2 == 0) {
            PyErr_SetString(PyExc_TypeError, "can't convert a list item to integer");
            ok = 0;
          } else {
            ok = (PyInt_AsInt(item, &value) == 0);
            Py_DECREF(item2);
          }
        }
       
        if (ok == 0) {
          igraph_vector_int_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        if (igraph_vector_int_push_back(v, value)) {
          igraphmodule_handle_igraph_error();
          igraph_vector_int_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        Py_DECREF(item);
      }
      Py_DECREF(it);
      return 0;
    } else {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      return 1;
    }
    return 0;
  }

  j=PySequence_Size(list);
  igraph_vector_int_init(v, j);
  for (i=0, k=0; i<j; i++) {
    item=PySequence_GetItem(list, i);
    if (item) {
      ok=1;
      if (!PyNumber_Check(item)) {
        PyErr_SetString(PyExc_TypeError, "sequence elements must be integers");
        ok=0;
      } else {
        PyObject *item2 = PyNumber_Int(item);
        if (item2 == 0) {
          PyErr_SetString(PyExc_TypeError, "can't convert sequence element to int");
          ok=0;
        } else {
          retval = PyInt_AsInt(item2, &value);
          if (retval)
            ok = 0;
          Py_DECREF(item2);
        }
      }
      Py_XDECREF(item);
      if (!ok) {
        igraph_vector_int_destroy(v);
        return 1;
      }
      VECTOR(*v)[k]=value;
      k++;
    } else {
      /* this should not happen, but we return anyway.
       * an IndexError exception was set by PyList_GetItem
       * at this point */
      igraph_vector_int_destroy(v);
      return 1;
    }
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of ints to an igraph \c igraph_vector_long_t
 * The incoming \c igraph_vector_long_t should be uninitialized.
 * Raises suitable Python exceptions when needed.
 *
 * This function is almost identical to
 * \ref igraphmodule_PyObject_to_vector_t . Make sure you fix bugs
 * in both cases (if any).
 *
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_long_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vector_long_t(PyObject *list, igraph_vector_long_t *v) {
  PyObject *item;
  long value=0;
  Py_ssize_t i, j, k;
  int ok;

  if (PyBaseString_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing integers");
    return 1;
  }

  if (!PySequence_Check(list)) {
    /* try to use an iterator */
    PyObject *it = PyObject_GetIter(list);
    if (it) {
      PyObject *item;
      igraph_vector_long_init(v, 0);
      while ((item = PyIter_Next(it)) != 0) {
        ok = 1;
        if (!PyNumber_Check(item)) {
          PyErr_SetString(PyExc_TypeError, "iterable must return numbers");
          ok=0;
        } else {
          PyObject *item2 = PyNumber_Int(item);
          if (item2 == 0) {
            PyErr_SetString(PyExc_TypeError, "can't convert a list item to integer");
            ok = 0;
          } else {
            value=(long)PyInt_AsLong(item);
            Py_DECREF(item2);
          }
        }
       
        if (ok == 0) {
          igraph_vector_long_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        if (igraph_vector_long_push_back(v, value)) {
          igraphmodule_handle_igraph_error();
          igraph_vector_long_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        Py_DECREF(item);
      }
      Py_DECREF(it);
      return 0;
    } else {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      return 1;
    }
    return 0;
  }

  j=PySequence_Size(list);
  igraph_vector_long_init(v, j);
  for (i=0, k=0; i<j; i++) {
    item=PySequence_GetItem(list, i);
    if (item) {
      ok=1;
      if (!PyNumber_Check(item)) {
        PyErr_SetString(PyExc_TypeError, "sequence elements must be integers");
        ok=0;
      } else {
        PyObject *item2 = PyNumber_Int(item);
        if (item2 == 0) {
          PyErr_SetString(PyExc_TypeError, "can't convert sequence element to integer");
          ok=0;
        } else {
          value=(long)PyInt_AsLong(item2);
          Py_DECREF(item2);
        }
      }
      Py_XDECREF(item);
      if (!ok) {
        igraph_vector_long_destroy(v);
        return 1;
      }
      VECTOR(*v)[k]=value;
      k++;
    } else {
      /* this should not happen, but we return anyway.
       * an IndexError exception was set by PyList_GetItem
       * at this point */
      igraph_vector_long_destroy(v);
      return 1;
    }
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of objects to an igraph \c igraph_vector_bool_t
 * The incoming \c igraph_vector_bool_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * 
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_bool_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vector_bool_t(PyObject *list,
    igraph_vector_bool_t *v) {
  PyObject *item;
  Py_ssize_t i, j;

  if (PyBaseString_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable");
    return 1;
  }

  if (!PySequence_Check(list)) {
    /* try to use an iterator */
    PyObject *it = PyObject_GetIter(list);
    if (it) {
      PyObject *item;
      igraph_vector_bool_init(v, 0);
      while ((item = PyIter_Next(it)) != 0) {
        if (igraph_vector_bool_push_back(v, PyObject_IsTrue(item))) {
          igraphmodule_handle_igraph_error();
          igraph_vector_bool_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        Py_DECREF(item);
      }
      Py_DECREF(it);
      return 0;
    } else {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      return 1;
    }
    return 0;
  }

  j=PySequence_Size(list);
  igraph_vector_bool_init(v, j);
  for (i=0; i<j; i++) {
    item=PySequence_GetItem(list, i);
    if (item) {
      VECTOR(*v)[i]=PyObject_IsTrue(item);
	  Py_DECREF(item);
    } else {
      /* this should not happen, but we return anyway.
       * an IndexError exception was set by PySequence_GetItem
       * at this point */
      igraph_vector_bool_destroy(v);
      return 1;
    }
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_bool_t to a Python boolean list
 * 
 * \param v the \c igraph_vector_bool_t containing the vector to be converted
 * \return the Python boolean list as a \c PyObject*, or \c NULL if an
 * error occurred
 */
PyObject* igraphmodule_vector_bool_t_to_PyList(const igraph_vector_bool_t *v) {
  PyObject *list, *item;
  Py_ssize_t n, i;
   
  n=igraph_vector_bool_size(v);
  if (n<0)
    return igraphmodule_handle_igraph_error();

  list=PyList_New(n);
  for (i=0; i<n; i++) {
    item = VECTOR(*v)[i] ? Py_True : Py_False;
    Py_INCREF(item);
    PyList_SET_ITEM(list, i, item);
  }

  return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python integer list
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyList(const igraph_vector_t *v,
    igraphmodule_conv_t type) {
  PyObject *list, *item;
  Py_ssize_t n, i;
   
  n=igraph_vector_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();

  list=PyList_New(n);
  for (i=0; i<n; i++) {
    if (type == IGRAPHMODULE_TYPE_INT) {
	  if (!igraph_finite(VECTOR(*v)[i]))
		item=PyFloat_FromDouble((double)VECTOR(*v)[i]);
	  else
        item=PyInt_FromLong((long)VECTOR(*v)[i]);
    } else if (type == IGRAPHMODULE_TYPE_FLOAT) {
      item=PyFloat_FromDouble((double)VECTOR(*v)[i]);
    } else {
      item=Py_None;
      Py_INCREF(item);
    }
    if (!item) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SET_ITEM(list, i, item);
  }

  return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_long_t to a Python integer list
 * 
 * \param v the \c igraph_vector_long_t containing the vector to be converted
 * \return the Python integer list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_long_t_to_PyList(const igraph_vector_long_t *v) {
  PyObject *list, *item;
  Py_ssize_t n, i;
   
  n = igraph_vector_long_size(v);
  if (n<0)
    return igraphmodule_handle_igraph_error();

  list=PyList_New(n);
  for (i=0; i<n; i++) {
    item = PyInt_FromLong(VECTOR(*v)[i]);
    if (!item) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SET_ITEM(list, i, item);
  }

  return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python integer tuple
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer tuple as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyTuple(const igraph_vector_t *v) {
  PyObject* tuple;
  Py_ssize_t n, i;
  
  n=igraph_vector_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();
  
  tuple=PyTuple_New(n);
  for (i=0; i<n; i++) {
    PyObject *item=PyInt_FromLong((long)VECTOR(*v)[i]);
    if (!item) {
      Py_DECREF(tuple);
      return NULL;
    }
    PyTuple_SET_ITEM(tuple, i, item);
  }

  return tuple;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python list of integer pairs
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer pair list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyList_pairs(const igraph_vector_t *v) {
   PyObject *list, *pair;
   long n, i, j;
   
   n=igraph_vector_size(v);
   if (n<0) return igraphmodule_handle_igraph_error();
   if (n%2) return igraphmodule_handle_igraph_error();
   
   /* create a new Python list */
   n>>=1;
   list=PyList_New(n);
   
   /* populate the list with data */
   for (i=0, j=0; i<n; i++, j+=2) {
     pair=Py_BuildValue("(ll)", (long)VECTOR(*v)[j], (long)VECTOR(*v)[j+1]);
     if (pair==NULL || PyList_SetItem(list, i, pair)) {
       /* error occurred while populating the list, return immediately */
       Py_DECREF(pair);
       Py_DECREF(list);
       return NULL;
     }
   }

   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python iterable of non-negative integer pairs (i.e. an
 * edge list) to an igraph \c igraph_vector_t
 *
 * The incoming \c igraph_vector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * 
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_t containing the result
 * \param graph  the graph that will be used to interpret vertex names
 *               if a string is yielded by the Python iterable
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_edgelist(PyObject *list, igraph_vector_t *v,
    igraph_t *graph) {
  PyObject *item, *i1, *i2, *it;
  int ok;
  igraph_integer_t idx1=0, idx2=0;

  if (PyBaseString_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing integer or string pairs");
    return 1;
  }

  it = PyObject_GetIter(list);
  if (!it)
    return 1;

  igraph_vector_init(v, 0);
  while ((item = PyIter_Next(it)) != 0) {
    ok = 1;
    if (!PySequence_Check(item) || PySequence_Size(item) != 2) {
      PyErr_SetString(PyExc_TypeError, "iterable must return pairs of integers or strings");
      ok=0;
    } else {
      i1 = PySequence_ITEM(item, 0);
      if (i1 == 0) {
        i2 = 0;
      } else {
        i2 = PySequence_ITEM(item, 1);
      }
      ok = (i1 != 0 && i2 != 0);
      ok = ok && !igraphmodule_PyObject_to_vid(i1, &idx1, graph);
      ok = ok && !igraphmodule_PyObject_to_vid(i2, &idx2, graph);
      Py_XDECREF(i1); Py_XDECREF(i2); /* PySequence_ITEM returned new ref */
    }

    Py_DECREF(item);

    if (ok) {
      if (igraph_vector_push_back(v, idx1)) {
        igraphmodule_handle_igraph_error();
        ok = 0;
      }
      if (ok && igraph_vector_push_back(v, idx2)) {
        igraphmodule_handle_igraph_error();
        ok = 0;
      }
    }

    if (!ok) {
      igraph_vector_destroy(v);
      Py_DECREF(it);
      return 1;
    }
  }

  Py_DECREF(it);
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an attribute name or a sequence to a vector_t
 *
 * This function is useful for the interface of igraph C calls accepting
 * edge or vertex weights. The function checks the given Python object. If
 * it is None, returns a null pointer instead of an \c igraph_vector_t.
 * If it is a sequence, it converts the sequence to a newly allocated
 * \c igraph_vector_t and return a pointer to it. Otherwise it interprets the
 * object as an attribute name and returns the attribute values corresponding
 * to the name as an \c igraph_vector_t, or returns a null pointer if the attribute
 * does not exist.
 * 
 * Note that if the function returned a pointer to an \c igraph_vector_t,
 * it is the caller's responsibility to destroy the object and free its
 * pointer after having finished using it.
 *
 * \param o the Python object being converted.
 * \param self a Python Graph object being used when attributes are queried
 * \param vptr the pointer to the allocated vector is returned here.
 * \param attr_type the type of the attribute being handled
 * \return 0 if everything was OK, nonzero otherwise.
 */
int igraphmodule_attrib_to_vector_t(PyObject *o, igraphmodule_GraphObject *self,
    igraph_vector_t **vptr, int attr_type) {
  igraph_vector_t *result;

  *vptr = 0;
  if (attr_type != ATTRIBUTE_TYPE_EDGE && attr_type != ATTRIBUTE_TYPE_VERTEX)
    return 1;
  if (o == Py_None) return 0;
  if (PyString_Check(o)) {
    /* Check whether the attribute exists and is numeric */
    igraph_attribute_type_t at;
    igraph_attribute_elemtype_t et;
    long int n;
    char *name = PyString_CopyAsString(o);

    if (attr_type == ATTRIBUTE_TYPE_VERTEX) {
      et = IGRAPH_ATTRIBUTE_VERTEX;
      n = igraph_vcount(&self->g);
    } else {
      et = IGRAPH_ATTRIBUTE_EDGE;
      n = igraph_ecount(&self->g);
    }

    if (igraphmodule_i_attribute_get_type(&self->g, &at, et, name)) {
      /* exception was set by igraphmodule_i_attribute_get_type */
      free(name);
      return 1;
    }
    if (at != IGRAPH_ATTRIBUTE_NUMERIC) {
      PyErr_SetString(PyExc_ValueError, "attribute values must be numeric");
      free(name);
      return 1;
    }
    /* Now that the attribute type has been checked, allocate the target
     * vector */
    result = (igraph_vector_t*)calloc(1, sizeof(igraph_vector_t));
    if (result==0) {
      PyErr_NoMemory();
      free(name);
      return 1;
    }
    igraph_vector_init(result, n);
    if (attr_type == ATTRIBUTE_TYPE_VERTEX) {
      if (igraphmodule_i_get_numeric_vertex_attr(&self->g, name,
          igraph_vss_all(), result)) {
        /* exception has already been set, so return */
        igraph_vector_destroy(result);
        free(name);
        free(result);
        return 1;
      }
    } else {
      if (igraphmodule_i_get_numeric_edge_attr(&self->g, name,
          igraph_ess_all(IGRAPH_EDGEORDER_ID), result)) {
        /* exception has already been set, so return */
        igraph_vector_destroy(result);
        free(name);
        free(result);
        return 1;
      }
    }
    free(name);
    *vptr = result;
  } else if (PySequence_Check(o)) {
    result = (igraph_vector_t*)calloc(1, sizeof(igraph_vector_t));
    if (result==0) {
      PyErr_NoMemory();
      return 1;
    }
    if (igraphmodule_PyObject_float_to_vector_t(o, result)) {
      igraph_vector_destroy(result);
      free(result);
      return 1;
    }
    *vptr = result;
  } else {
    PyErr_SetString(PyExc_TypeError, "unhandled type");
    return 1;
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an attribute name or a sequence to a vector_int_t
 *
 * Similar to igraphmodule_attrib_to_vector_t and
 * igraphmodule_attrib_to_vector_long_t. Make sure you fix bugs
 * in all three places (if any).
 * 
 * Note that if the function returned a pointer to an \c igraph_vector_int_t,
 * it is the caller's responsibility to destroy the object and free its
 * pointer after having finished using it.
 *
 * \param o the Python object being converted.
 * \param self a Python Graph object being used when attributes are queried
 * \param vptr the pointer to the allocated vector is returned here.
 * \param attr_type the type of the attribute being handled
 * \return 0 if everything was OK, nonzero otherwise.
 */
int igraphmodule_attrib_to_vector_int_t(PyObject *o, igraphmodule_GraphObject *self,
    igraph_vector_int_t **vptr, int attr_type) {
  igraph_vector_int_t *result;

  *vptr = 0;

  if (attr_type != ATTRIBUTE_TYPE_EDGE && attr_type != ATTRIBUTE_TYPE_VERTEX)
    return 1;

  if (o == Py_None)
    return 0;

  if (PyString_Check(o)) {
    igraph_vector_t* dummy = 0;
    long int i, n;

    if (igraphmodule_attrib_to_vector_t(o, self, &dummy, attr_type))
      return 1;

    if (dummy == 0)
      return 0;

    n = igraph_vector_size(dummy);

    result = (igraph_vector_int_t*)calloc(1, sizeof(igraph_vector_int_t));
    igraph_vector_int_init(result, n);
    if (result==0) {
      igraph_vector_destroy(dummy); free(dummy);
      PyErr_NoMemory();
      return 1;
    }
    for (i=0; i<n; i++)
      VECTOR(*result)[i] = (int)VECTOR(*dummy)[i];
    igraph_vector_destroy(dummy); free(dummy);
    *vptr = result;
  } else if (PySequence_Check(o)) {
    result = (igraph_vector_int_t*)calloc(1, sizeof(igraph_vector_int_t));
    if (result==0) {
      PyErr_NoMemory();
      return 1;
    }
    if (igraphmodule_PyObject_to_vector_int_t(o, result)) {
      igraph_vector_int_destroy(result);
      free(result);
      return 1;
    }
    *vptr = result;
  } else {
    PyErr_SetString(PyExc_TypeError, "unhandled type");
    return 1;
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an attribute name or a sequence to a vector_long_t
 *
 * Similar to igraphmodule_attrib_to_vector_t and
 * igraphmodule_attrib_to_vector_int_t. Make sure you fix bugs
 * in all three places (if any).
 * 
 * Note that if the function returned a pointer to an \c igraph_vector_long_t,
 * it is the caller's responsibility to destroy the object and free its
 * pointer after having finished using it.
 *
 * \param o the Python object being converted.
 * \param self a Python Graph object being used when attributes are queried
 * \param vptr the pointer to the allocated vector is returned here.
 * \param attr_type the type of the attribute being handled
 * \return 0 if everything was OK, nonzero otherwise.
 */
int igraphmodule_attrib_to_vector_long_t(PyObject *o, igraphmodule_GraphObject *self,
    igraph_vector_long_t **vptr, int attr_type) {
  igraph_vector_long_t *result;

  *vptr = 0;
  if (attr_type != ATTRIBUTE_TYPE_EDGE && attr_type != ATTRIBUTE_TYPE_VERTEX)
    return 1;
  if (o == Py_None)
    return 0;
  if (PyString_Check(o)) {
    igraph_vector_t* dummy = 0;
    long int i, n;

    if (igraphmodule_attrib_to_vector_t(o, self, &dummy, attr_type))
      return 1;

    if (dummy == 0)
      return 0;

    n = igraph_vector_size(dummy);

    result = (igraph_vector_long_t*)calloc(1, sizeof(igraph_vector_long_t));
    igraph_vector_long_init(result, n);
    if (result==0) {
      igraph_vector_destroy(dummy); free(dummy);
      PyErr_NoMemory();
      return 1;
    }
    for (i=0; i<n; i++)
      VECTOR(*result)[i] = (long int)VECTOR(*dummy)[i];
    igraph_vector_destroy(dummy); free(dummy);
    *vptr = result;
  } else if (PySequence_Check(o)) {
    result = (igraph_vector_long_t*)calloc(1, sizeof(igraph_vector_long_t));
    if (result==0) {
      PyErr_NoMemory();
      return 1;
    }
    if (igraphmodule_PyObject_to_vector_long_t(o, result)) {
      igraph_vector_long_destroy(result);
      free(result);
      return 1;
    }
    *vptr = result;
  } else {
    PyErr_SetString(PyExc_TypeError, "unhandled type");
    return 1;
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an attribute name or a sequence to a vector_bool_t
 *
 * This function is useful for the interface of igraph C calls accepting
 * bipartite type vectors. The function checks the given Python object. If
 * it is None, returns a null pointer instead of an \c igraph_vector_bool_t.
 * If it is a sequence, it converts the sequence to a newly allocated
 * \c igraph_vector_bool_t and return a pointer to it. Otherwise it interprets the
 * object as an attribute name and returns the attribute values corresponding
 * to the name as an \c igraph_vector_bool_t, or returns a null pointer if the attribute
 * does not exist.
 *
 * Anything that evaluates to true using PyObject_IsTrue will be converted to
 * true in the resulting vector. Only numeric attributes are supported.
 * 
 * Note that if the function returned a pointer to an \c igraph_vector_bool_t,
 * it is the caller's responsibility to destroy the object and free its
 * pointer after having finished using it.
 *
 * \param o the Python object being converted.
 * \param self a Python Graph object being used when attributes are queried
 * \param vptr the pointer to the allocated vector is returned here.
 * \param attr_type the type of the attribute being handled
 * \return 0 if everything was OK, nonzero otherwise.
 */
int igraphmodule_attrib_to_vector_bool_t(PyObject *o, igraphmodule_GraphObject *self,
    igraph_vector_bool_t **vptr, int attr_type) {
  igraph_vector_bool_t *result;

  *vptr = 0;
  if (attr_type != ATTRIBUTE_TYPE_EDGE && attr_type != ATTRIBUTE_TYPE_VERTEX)
    return 1;
  if (o == Py_None) return 0;
  if (PyString_Check(o)) {
    igraph_vector_t *dummy = 0;
    long int i, n;

    if (igraphmodule_attrib_to_vector_t(o, self, &dummy, attr_type))
      return 1;

    if (dummy == 0)
      return 0;

    n = igraph_vector_size(dummy);

    result = (igraph_vector_bool_t*)calloc(1, sizeof(igraph_vector_bool_t));
    igraph_vector_bool_init(result, n);
    if (result==0) {
      igraph_vector_destroy(dummy); free(dummy);
      PyErr_NoMemory();
      return 1;
    }
    for (i=0; i<n; i++)
      VECTOR(*result)[i] = (VECTOR(*dummy)[i] != 0);
    igraph_vector_destroy(dummy); free(dummy);
    *vptr = result;
  } else if (PySequence_Check(o)) {
    result = (igraph_vector_bool_t*)calloc(1, sizeof(igraph_vector_bool_t));
    if (result==0) {
      PyErr_NoMemory();
      return 1;
    }
    if (igraphmodule_PyObject_to_vector_bool_t(o, result)) {
      free(result);
      return 1;
    }
    *vptr = result;
  } else {
    PyErr_SetString(PyExc_TypeError, "unhandled type");
    return 1;
  }

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts two igraph \c igraph_vector_t vectors to a Python list of integer pairs
 * 
 * \param v1 the \c igraph_vector_t containing the 1st vector to be converted
 * \param v2 the \c igraph_vector_t containing the 2nd vector to be converted
 * \return the Python integer pair list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_pair_to_PyList(const igraph_vector_t *v1,
    const igraph_vector_t *v2) {
   PyObject *list, *pair;
   long n, i;
   
   n=igraph_vector_size(v1);
   if (n<0) return igraphmodule_handle_igraph_error();
   if (igraph_vector_size(v2) != n) return igraphmodule_handle_igraph_error();

   /* create a new Python list */
   list=PyList_New(n);
   
   /* populate the list with data */
   for (i=0; i<n; i++) {
     pair=Py_BuildValue("(ll)", (long)VECTOR(*v1)[i], (long)VECTOR(*v2)[i]);
     if (pair==NULL || PyList_SetItem(list, i, pair)) {
       /* error occurred while populating the list, return immediately */
       Py_DECREF(pair);
       Py_DECREF(list);
       return NULL;
     }
   }

   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_matrix_t to a Python list of lists
 * 
 * \param m the \c igraph_matrix_t containing the matrix to be converted
 * \param type the type of conversion. If equals to IGRAPHMODULE_TYPE_INT,
 *        returns an integer matrix, else returns a float matrix.
 * \return the Python list of lists as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_matrix_t_to_PyList(const igraph_matrix_t *m,
    igraphmodule_conv_t type) {
   PyObject *list, *row, *item;
   Py_ssize_t nr, nc, i, j;
   
   nr = igraph_matrix_nrow(m);
   nc = igraph_matrix_ncol(m);
   if (nr<0 || nc<0)
     return igraphmodule_handle_igraph_error();

   // create a new Python list
   list=PyList_New(nr);
   // populate the list with data
   for (i=0; i<nr; i++) 
     {
    row=PyList_New(nc);
    for (j=0; j<nc; j++) 
      {
         if (type==IGRAPHMODULE_TYPE_INT) {
	       if (!igraph_finite(MATRIX(*m, i, j)))
		     item=PyFloat_FromDouble((double)MATRIX(*m, i, j));
	       else
             item=PyInt_FromLong((long)MATRIX(*m, i, j));
		 } else
           item=PyFloat_FromDouble(MATRIX(*m, i, j));
           
         if (PyList_SetItem(row, j, item))
           {
          // error occurred while populating the list, return immediately
          Py_DECREF(row);
          Py_DECREF(list);
          return NULL;
           }
      }
    if (PyList_SetItem(list, i, row)) 
      {
         Py_DECREF(row);
         Py_DECREF(list);
         return NULL;
      }
     }
   // return the list
   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_ptr_t to a Python list of lists
 * 
 * \param v the \c igraph_vector_ptr_t containing the vector to be converted
 * \return the Python list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_ptr_t_to_PyList(const igraph_vector_ptr_t *v,
    igraphmodule_conv_t type) {
  PyObject *list, *item;
  Py_ssize_t n, i;
   
  n=igraph_vector_ptr_size(v);
  if (n<0)
    return igraphmodule_handle_igraph_error();

  list=PyList_New(n);
  for (i=0; i<n; i++) {
    item=igraphmodule_vector_t_to_PyList((igraph_vector_t*)VECTOR(*v)[i], type);
    if (item == NULL) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SET_ITEM(list, i, item);
  }

  return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of lists to an \c igraph_matrix_t
 * 
 * \param o the Python object representing the list of lists
 * \param m the address of an uninitialized \c igraph_matrix_t
 * \return 0 if everything was OK, 1 otherwise. Sets appropriate exceptions.
 */
int igraphmodule_PyList_to_matrix_t(PyObject* o, igraph_matrix_t *m) {
  Py_ssize_t nr, nc, n, i, j;
  PyObject *row, *item;
  int was_warned=0;

  /* calculate the matrix dimensions */
  if (!PySequence_Check(o) || PyString_Check(o)) {
    PyErr_SetString(PyExc_TypeError, "matrix expected (list of sequences)");
    return 1;
  }

  nr = PySequence_Size(o);
  nc = 0;
  for (i=0; i<nr; i++) {
    row=PySequence_GetItem(o, i);
    if (!PySequence_Check(row)) {
      Py_DECREF(row);
      PyErr_SetString(PyExc_TypeError, "matrix expected (list of sequences)");
      return 1;
    }
    n=PySequence_Size(row);
    Py_DECREF(row);
    if (n>nc) nc=n;
  }
  
  igraph_matrix_init(m, nr, nc);
  for (i=0; i<nr; i++) {
    row=PySequence_GetItem(o, i);
    n=PySequence_Size(row);
    for (j=0; j<n; j++) {
      item=PySequence_GetItem(row, j);
      if (PyInt_Check(item)) {
        MATRIX(*m, i, j) = (igraph_real_t)PyInt_AsLong(item);
      } else if (PyLong_Check(item)) {
        MATRIX(*m, i, j) = (igraph_real_t)PyLong_AsLong(item);
      } else if (PyFloat_Check(item)) {
        MATRIX(*m, i, j) = (igraph_real_t)PyFloat_AsDouble(item);
      } else if (!was_warned) {
        PyErr_Warn(PyExc_Warning, "non-numeric value in matrix ignored");
        was_warned=1;
      }
      Py_DECREF(item);
    }
    Py_DECREF(row);
  }

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an \c igraph_strvector_t to a Python string list
 * 
 * \param v the \c igraph_strvector_t containing the vector to be converted
 * \return the Python string list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_strvector_t_to_PyList(igraph_strvector_t *v) {
  PyObject* list;
  Py_ssize_t n, i;
  char* ptr;
  
  n=igraph_strvector_size(v);
  if (n<0)
    return igraphmodule_handle_igraph_error();
  
  // create a new Python list
  list=PyList_New(n);
  /* populate the list with data */
  for (i=0; i<n; i++) {
    igraph_strvector_get(v, i, &ptr);
    if (PyList_SetItem(list, i, PyString_FromString(ptr))) {
      /* error occurred while populating the list, return immediately */
      Py_DECREF(list);
      return NULL;
    }
  }
  
   /* return the list */
   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python string list to an \c igraph_strvector_t
 * 
 * \param v the Python list as a \c PyObject*
 * \param result an \c igraph_strvector_t containing the result
 * The incoming \c igraph_strvector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyList_to_strvector_t(PyObject* v, igraph_strvector_t *result) {
  Py_ssize_t n, i;
  PyObject *o;
  
  if (!PyList_Check(v)) {
    PyErr_SetString(PyExc_TypeError, "expected list");
    return 1;
  }
  
  n=PyList_Size(v);

  /* initialize the string vector */
  if (igraph_strvector_init(result, n)) return 1;
  
  /* populate the vector with data */
  for (i=0; i<n; i++) {
    PyObject *item = PyList_GetItem(v, i);
    char* ptr;
    igraph_bool_t will_free = 0;

    if (PyUnicode_Check(item)) {
      ptr = PyString_CopyAsString(item);
      if (ptr == 0) {
        igraph_strvector_destroy(result);
        return 1;
      }
      will_free = 1;
#ifndef IGRAPH_PYTHON3
    } else if (PyString_Check(item)) {
      ptr = PyString_AS_STRING(item);
#endif
    } else {
      o = PyObject_Str(item);
      if (o == 0) {
        igraph_strvector_destroy(result);
        return 1;
      }
      ptr = PyString_CopyAsString(o);
      Py_DECREF(o);
      if (ptr == 0) {
        igraph_strvector_destroy(result);
        return 1;
      }
      will_free = 1;
    }

    if (igraph_strvector_set(result, i, ptr)) {
      if (will_free)
        free(ptr);
      igraph_strvector_destroy(result);
      return 1;
    }
    if (will_free)
      free(ptr);
  }

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Appends the contents of a Python iterator returning graphs to
 * an \c igraph_vectorptr_t
 *
 * The incoming \c igraph_vector_ptr_t should be INITIALIZED.
 * Raises suitable Python exceptions when needed.
 * 
 * \param it the Python iterator
 * \param v the \c igraph_vector_ptr_t which will contain the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_append_PyIter_to_vector_ptr_t(PyObject *it, igraph_vector_ptr_t *v) {
  PyObject *t;
  
  while ((t=PyIter_Next(it))) {
    if (!PyObject_TypeCheck(t, &igraphmodule_GraphType)) {
      PyErr_SetString(PyExc_TypeError, "iterable argument must contain graphs");
      Py_DECREF(t);
      return 1;
    }
    igraph_vector_ptr_push_back(v, &((igraphmodule_GraphObject*)t)->g);
    Py_DECREF(t);
  }  
  
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as a single vertex ID
 * 
 * \param o      the Python object
 * \param vid    the vertex ID will be stored here
 * \param graph  the graph that will be used to interpret vertex names
 *               if a string was given in o. It may also be a null pointer
 *               if we don't need name lookups.
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vid(PyObject *o, igraph_integer_t *vid, igraph_t *graph) {
  int retval, tmp;

  if (o == Py_None || o == 0) {
    *vid = 0;
  } else if (PyInt_Check(o)) {
    /* Single vertex ID */
    if (PyInt_AsInt(o, &tmp))
      return 1;
    *vid = tmp;
  } else if (PyLong_Check(o)) {
    /* Single vertex ID */
    if (PyLong_AsInt(o, &tmp))
      return 1;
    *vid = tmp;
  } else if (graph != 0 && PyBaseString_Check(o)) {
    /* Single vertex ID from vertex name */
    if (igraphmodule_get_vertex_id_by_name(graph, o, vid))
      return 1;
  } else if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_VertexType)) {
    /* Single vertex ID from Vertex object */
    igraphmodule_VertexObject *vo = (igraphmodule_VertexObject*)o;
    *vid = igraphmodule_Vertex_get_index_igraph_integer(vo);
  } else if (PyIndex_Check(o)) {
    /* Other numeric type that can be converted to an index */
    PyObject* num = PyNumber_Index(o);
    if (num) {
      if (PyInt_Check(num)) {
        retval = PyInt_AsInt(num, &tmp);
        if (retval) {
          Py_DECREF(num);
          return 1;
        }
        *vid = tmp;
      } else if (PyLong_Check(num)) {
        retval = PyLong_AsInt(num, &tmp);
        if (retval) {
          Py_DECREF(num);
          return 1;
        }
        *vid = tmp;
      } else {
        PyErr_SetString(PyExc_TypeError, "PyNumber_Index returned invalid type");
        Py_DECREF(num);
        return 1;
      }
      Py_DECREF(num);
    } else
      return 1;
  } else {
    PyErr_SetString(PyExc_TypeError, "only numbers, vertex names or igraph.Vertex objects can be converted to vertex IDs");
    return 1;
  }

  if (*vid < 0) {
    PyErr_Format(PyExc_ValueError, "vertex IDs must be positive, got: %ld", (long)(*vid));
    return 1;
  }

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as a vertex selector
 * 
 * \param o      the Python object
 * \param vs     the \c igraph_vs_t which will contain the result
 * \param graph  an \c igraph_t object which will be used to interpret vertex
 *               names (if the supplied Python object contains strings)
 * \param return_single will be 1 if the selector selected only a single vertex,
 *                      0 otherwise
 * \param single_vid    if the selector selected only a single vertex, the ID
 *                      of the selected vertex will also be returned here.
 *
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vs_t(PyObject *o, igraph_vs_t *vs,
    igraph_t *graph, igraph_bool_t *return_single, igraph_integer_t *single_vid) {
  igraph_integer_t vid;
  igraph_vector_t vector;

  if (o == 0 || o == Py_None) {
    /* Returns a vertex sequence for all vertices */
    if (return_single)
      *return_single = 0;
    igraph_vs_all(vs);
    return 0;
  }

  if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_VertexSeqType)) {
    /* Returns a vertex sequence from a VertexSeq object */
    igraphmodule_VertexSeqObject *vso = (igraphmodule_VertexSeqObject*)o;
    if (igraph_vs_copy(vs, &vso->vs)) {
      igraphmodule_handle_igraph_error();
      return 1;
    }
    if (return_single)
      *return_single = 0;
    return 0;
  }

  if (PySlice_Check(o) && graph != 0) {
    /* Returns a vertex sequence from a slice */
    Py_ssize_t no_of_vertices = igraph_vcount(graph);
    Py_ssize_t start, stop, step, slicelength, i;

    /* Casting to void* because Python 2.x expects PySliceObject*
     * but Python 3.x expects PyObject* */
    if (PySlice_GetIndicesEx((void*)o, no_of_vertices,
          &start, &stop, &step, &slicelength))
      return 1;

    if (start == 0 && slicelength == no_of_vertices) {
      igraph_vs_all(vs);
    } else {
      IGRAPH_CHECK(igraph_vector_init(&vector, slicelength));
      IGRAPH_FINALLY(igraph_vector_destroy, &vector);

      for (i = 0; i < slicelength; i++, start += step) {
        VECTOR(vector)[i] = start;
      }

      IGRAPH_CHECK(igraph_vs_vector_copy(vs, &vector));

      igraph_vector_destroy(&vector);
      IGRAPH_FINALLY_CLEAN(1);
    }

    if (return_single)
      *return_single = 0;

    return 0;
  }

  if (igraphmodule_PyObject_to_vid(o, &vid, graph)) {
    /* Object cannot be converted to a single vertex ID,
     * assume it is a sequence or iterable */

    PyObject *iterator;
    PyObject *item;

    if (PyBaseString_Check(o)) {
      /* Special case: strings and unicode objects are sequences, but they
       * will not yield valid vertex IDs */
      return 1;
    }

    /* Clear the exception set by igraphmodule_PyObject_to_vid */
    PyErr_Clear();

    iterator = PyObject_GetIter(o);

    if (iterator == NULL) {
      PyErr_SetString(PyExc_TypeError, "conversion to vertex sequence failed");
      return 1;
    }

    IGRAPH_CHECK(igraph_vector_init(&vector, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &vector);
    IGRAPH_CHECK(igraph_vector_reserve(&vector, 20));

    while ((item = PyIter_Next(iterator))) {
      vid = -1;

      if (igraphmodule_PyObject_to_vid(item, &vid, graph))
        break;

      Py_DECREF(item);
      igraph_vector_push_back(&vector, vid);
    }
    Py_DECREF(iterator);

    if (PyErr_Occurred()) {
      igraph_vector_destroy(&vector);
      IGRAPH_FINALLY_CLEAN(1);
      return 1;
    }

    IGRAPH_CHECK(igraph_vs_vector_copy(vs, &vector));
    igraph_vector_destroy(&vector);
    IGRAPH_FINALLY_CLEAN(1);

    if (return_single)
      *return_single = 0;
    
    return 0;
  }

  /* The object can be converted into a single vertex ID */
  if (return_single)
    *return_single = 1;
  if (single_vid)
    *single_vid = vid;

  igraph_vs_1(vs, vid);

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as a single edge ID
 * 
 * \param o      the Python object
 * \param eid    the edge ID will be stored here
 * \param graph  the graph that will be used to interpret vertex names and
 *               indices if o is a tuple. It may also be a null pointer
 *               if we don't want to handle tuples.
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_eid(PyObject *o, igraph_integer_t *eid, igraph_t *graph) {
  int retval, tmp;
  igraph_integer_t vid1, vid2;

  if (o == Py_None || o == 0) {
    *eid = 0;
  } else if (PyInt_Check(o)) {
    /* Single edge ID */
    if (PyInt_AsInt(o, &tmp))
      return 1;
    *eid = tmp;
  } else if (PyLong_Check(o)) {
    /* Single edge ID */
    if (PyLong_AsInt(o, &tmp))
      return 1;
    *eid = tmp;
  } else if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_EdgeType)) {
    /* Single edge ID from Edge object */
    igraphmodule_EdgeObject *eo = (igraphmodule_EdgeObject*)o;
    *eid = igraphmodule_Edge_get_index_igraph_integer(eo);
  } else if (PyIndex_Check(o)) {
    /* Other numeric type that can be converted to an index */
    PyObject* num = PyNumber_Index(o);
    if (num) {
      if (PyInt_Check(num)) {
        retval = PyInt_AsInt(num, &tmp);
        if (retval) {
          Py_DECREF(num);
          return 1;
        }
        *eid = tmp;
      } else if (PyLong_Check(num)) {
        retval = PyLong_AsInt(num, &tmp);
        if (retval) {
          Py_DECREF(num);
          return 1;
        }
        *eid = tmp;
      } else {
        PyErr_SetString(PyExc_TypeError, "PyNumber_Index returned invalid type");
        Py_DECREF(num);
        return 1;
      }
      Py_DECREF(num);
    } else
      return 1;
  } else if (graph != 0 && PyTuple_Check(o)) {
    PyObject *o1, *o2;
    
    o1 = PyTuple_GetItem(o, 0);
    if (!o1)
      return 1;

    o2 = PyTuple_GetItem(o, 1);
    if (!o2)
      return 1;

    if (igraphmodule_PyObject_to_vid(o1, &vid1, graph))
      return 1;
    if (igraphmodule_PyObject_to_vid(o2, &vid2, graph))
      return 1;

    igraph_get_eid(graph, eid, vid1, vid2, 1, 0);
    if (*eid < 0) {
      PyErr_Format(PyExc_ValueError, "no edge from vertex #%ld to #%ld",
          (long int)vid1, (long int)vid2);
      return 1;
    }
  } else {
    PyErr_SetString(PyExc_TypeError,
        "only numbers, igraph.Edge objects or tuples of vertex IDs can be "
        "converted to edge IDs");
    return 1;
  }

  if (*eid < 0) {
    PyErr_Format(PyExc_ValueError, "edge IDs must be positive, got: %ld", (long)(*eid));
    return 1;
  }

  return 0;
}


/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as an edge selector
 * 
 * \param o the Python object
 * \param vs the \c igraph_es_t which will contain the result
 * \param graph  an \c igraph_t object which will be used to interpret vertex
 *               names and tuples (if the supplied Python object contains them)
 * \param return_single will be 1 if the selector selected only a single edge,
 * 0 otherwise
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_es_t(PyObject *o, igraph_es_t *es, igraph_t *graph,
                  igraph_bool_t *return_single) {
  igraph_integer_t eid;
  igraph_vector_t vector;

  if (o == 0 || o == Py_None) {
    /* Returns an edge sequence for all edges */
    if (return_single)
      *return_single = 0;
    igraph_es_all(es, IGRAPH_EDGEORDER_ID);
    return 0;
  }

  if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_EdgeSeqType)) {
    /* Returns an edge sequence from an EdgeSeq object */
    igraphmodule_EdgeSeqObject *eso = (igraphmodule_EdgeSeqObject*)o;
    if (igraph_es_copy(es, &eso->es)) {
      igraphmodule_handle_igraph_error();
      return 1;
    }
    if (return_single)
      *return_single = 0;
    return 0;
  }

  if (igraphmodule_PyObject_to_eid(o, &eid, graph)) {
    /* Object cannot be converted to a single edge ID,
     * assume it is a sequence or iterable */

    PyObject *iterator;
    PyObject *item;

    /* Clear the exception set by igraphmodule_PyObject_to_eid */
    PyErr_Clear();

    iterator = PyObject_GetIter(o);

    if (iterator == NULL) {
      PyErr_SetString(PyExc_TypeError, "conversion to edge sequene failed");
      return 1;
    }

    IGRAPH_CHECK(igraph_vector_init(&vector, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &vector);
    IGRAPH_CHECK(igraph_vector_reserve(&vector, 20));

    while ((item = PyIter_Next(iterator))) {
      eid = -1;

      if (igraphmodule_PyObject_to_eid(item, &eid, graph))
        break;

      Py_DECREF(item);
      igraph_vector_push_back(&vector, eid);
    }
    Py_DECREF(iterator);

    if (PyErr_Occurred()) {
      igraph_vector_destroy(&vector);
      IGRAPH_FINALLY_CLEAN(1);
      return 1;
    }
    
    if (igraph_vector_size(&vector) > 0) {
      igraph_es_vector_copy(es, &vector);
    } else {
      igraph_es_none(es);
    }

    igraph_vector_destroy(&vector);
    IGRAPH_FINALLY_CLEAN(1);

    if (return_single)
      *return_single = 0;
    
    return 0;
  }

  /* The object can be converted into a single edge ID */
  if (return_single)
    *return_single = 1;
  /*
  if (single_eid)
    *single_eid = eid;
  */

  igraph_es_1(es, eid);

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as a numeric attribute value list
 * 
 * \param o the Python object
 * \param v the \c igraph_vector_t which will contain the result
 * \param g a \c igraphmodule_GraphObject object or \c NULL - used when the
 * provided Python object is not a list and we're trying to interpret it as
 * an attribute name.
 * \param type the attribute type (graph = 0, vertex = 1, edge = 2) to be used
 * \param def default value if the attribute name supplied is \c None
 * if \c o is not a list.
 * \return 0 if everything was OK, 1 otherwise
 *
 * If the Python object is not a list, tries to interpret it as an attribute
 * name.
 */
int igraphmodule_PyObject_to_attribute_values(PyObject *o,
                          igraph_vector_t *v,
                          igraphmodule_GraphObject* g,
                          int type, igraph_real_t def) {
  PyObject* list = o;
  long i, n;

  if (o==NULL) return 1;
  
  if (o == Py_None) {
    if (type == ATTRHASH_IDX_VERTEX) n=igraph_vcount(&g->g);
    else if (type == ATTRHASH_IDX_EDGE) n=igraph_ecount(&g->g);
    else n=1;

    if (igraph_vector_init(v, n)) return 1;
    for (i=0; i<n; i++) VECTOR(*v)[i] = def;
    return 0;
  }

  if (!PyList_Check(o)) {
    list = PyDict_GetItem(((PyObject**)g->g.attr)[type], o);
    if (!list) {
      if (!PyErr_Occurred())
    PyErr_SetString(PyExc_KeyError, "Attribute does not exist");
      return 1;
    }
  }

  n=PyList_Size(list);
  if (igraph_vector_init(v, n)) return 1;

  for (i=0; i<n; i++) {
    PyObject *item = PyList_GetItem(list, i);
    if (!item) {
      igraph_vector_destroy(v);
      return 1;
    }

    if (PyInt_Check(item))
      VECTOR(*v)[i] = PyInt_AsLong(item);
    else if (PyLong_Check(item))
      VECTOR(*v)[i] = PyLong_AsLong(item);
    else if (PyFloat_Check(item))
      VECTOR(*v)[i] = PyFloat_AsDouble(item);
    else
      VECTOR(*v)[i] = def;
  }

  return 0;
}


int igraphmodule_PyObject_to_drl_options_t(PyObject *obj,
    igraph_layout_drl_options_t *options) {
  if (obj == Py_None) {
    igraph_layout_drl_options_init(options, IGRAPH_LAYOUT_DRL_DEFAULT);
  } else if (PyString_Check(obj)) {
    /* We have a string, so we are using a preset */
    igraph_layout_drl_default_t def=IGRAPH_LAYOUT_DRL_DEFAULT;
    if (PyString_IsEqualToASCIIString(obj, "default"))
      def=IGRAPH_LAYOUT_DRL_DEFAULT;
    else if (PyString_IsEqualToASCIIString(obj, "coarsen"))
      def=IGRAPH_LAYOUT_DRL_COARSEN;
    else if (PyString_IsEqualToASCIIString(obj, "coarsest"))
      def=IGRAPH_LAYOUT_DRL_COARSEST;
    else if (PyString_IsEqualToASCIIString(obj, "refine"))
      def=IGRAPH_LAYOUT_DRL_REFINE;
    else if (PyString_IsEqualToASCIIString(obj, "final"))
      def=IGRAPH_LAYOUT_DRL_FINAL;
    else {
      PyErr_SetString(PyExc_ValueError, "unknown DrL template name. Must be one of: default, coarsen, coarsest, refine, final");
      return 1;
    }
    if (igraph_layout_drl_options_init(options, def)) {
      igraphmodule_handle_igraph_error();
      return 1;
    }
  } else {
    igraph_layout_drl_options_init(options, IGRAPH_LAYOUT_DRL_DEFAULT);
#define CONVERT_DRL_OPTION(OPTION, TYPE) do { \
      PyObject *o1; \
      if (PyMapping_Check(obj)) { \
        o1 = PyMapping_GetItemString(obj, #OPTION); \
        igraphmodule_PyObject_to_##TYPE##_t(o1, &options->OPTION); \
        Py_XDECREF(o1); \
      } \
      o1 = PyObject_GetAttrString(obj, #OPTION); \
      igraphmodule_PyObject_to_##TYPE##_t(o1, &options->OPTION); \
      Py_XDECREF(o1); \
    } while (0)
#define CONVERT_DRL_OPTION_BLOCK(NAME) do { \
	  CONVERT_DRL_OPTION(NAME##_iterations, integer); \
	  CONVERT_DRL_OPTION(NAME##_temperature, real); \
	  CONVERT_DRL_OPTION(NAME##_attraction, real); \
	  CONVERT_DRL_OPTION(NAME##_damping_mult, real); \
    } while (0)
	
    CONVERT_DRL_OPTION(edge_cut, real);
    CONVERT_DRL_OPTION_BLOCK(init);
    CONVERT_DRL_OPTION_BLOCK(liquid);
    CONVERT_DRL_OPTION_BLOCK(expansion);
    CONVERT_DRL_OPTION_BLOCK(cooldown);
    CONVERT_DRL_OPTION_BLOCK(crunch);
    CONVERT_DRL_OPTION_BLOCK(simmer);

#undef CONVERT_DRL_OPTION
#undef CONVERT_DRL_OPTION_BLOCK

    PyErr_Clear();
	return 0;
  }
  return 0;
}


int igraphmodule_i_PyObject_pair_to_attribute_combination_record_t(
    PyObject* name, PyObject* value,
    igraph_attribute_combination_record_t *result) {
  if (igraphmodule_PyObject_to_attribute_combination_type_t(value, &result->type))
    return 1;

  if (result->type == IGRAPH_ATTRIBUTE_COMBINE_FUNCTION) {
    result->func = value;
  } else {
    result->func = 0;
  }

  if (name == Py_None)
    result->name = 0;
  else if (!PyString_Check(name)) {
    PyErr_SetString(PyExc_TypeError, "keys must be strings or None in attribute combination specification dicts");
    return 1;
  } else {
#ifdef IGRAPH_PYTHON3
    result->name = PyString_CopyAsString(name);
#else
    result->name = PyString_AS_STRING(name);
#endif
  }

  return 0;
}

/**
 * \brief Converts a Python object to an \c igraph_attribute_combination_t
 *
 * Raises suitable Python exceptions when needed.
 *
 * An \c igraph_attribute_combination_t specifies how the attributes of multiple
 * vertices/edges should be combined when they are collapsed into a single vertex
 * or edge (e.g., when simplifying a graph). For each attribute, one can specify
 * a Python callable object to call or one of a list of recognised strings which
 * map to simple functions. The recognised strings are as follows:
 *
 *   - \c "ignore"  - the attribute will be ignored
 *   - \c "sum"     - the attribute values will be added
 *   - \c "prod"    - the product of the attribute values will be taken
 *   - \c "min"     - the minimum attribute value will be used
 *   - \c "max"     - the maximum attribute value will be used
 *   - \c "random"  - a random value will be selected
 *   - \c "first"   - the first value encountered will be selected
 *   - \c "last"    - the last value encountered will be selected
 *   - \c "mean"    - the mean of the attributes will be selected
 *   - \c "median"  - the median of the attributes will be selected
 *   - \c "concat"  - the attribute values will be concatenated
 *
 * The Python object being converted must either be a string, a callable or a dict.
 * If a string is given, it is considered as an \c igraph_attribute_combination_t
 * object that combines all attributes according to the function given by that
 * string. If a callable is given, it is considered as an
 * \c igraph_attribute_combination_t that combines all attributes by calling the
 * callable and taking its return value. If a dict is given, its key-value pairs
 * are iterated, the keys specify the attribute names (a key of None means all
 * explicitly not specified attributes), the values specify the functions to
 * call for those attributes.
 *
 * \param object the Python object to be converted
 * \param result the result is returned here. It must be an uninitialized
 *   \c igraph_attribute_combination_t object, it will be initialized accordingly.
 *   It is the responsibility of the caller to 
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_attribute_combination_t(PyObject* object,
    igraph_attribute_combination_t *result) {
  igraph_attribute_combination_record_t rec;

  if (igraph_attribute_combination_init(result)) {
    igraphmodule_handle_igraph_error();
    return 1;
  }

  if (object == Py_None) {
    return 0;
  }

  if (PyDict_Check(object)) {
    /* a full-fledged dict was passed */
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(object, &pos, &key, &value)) {
      if (igraphmodule_i_PyObject_pair_to_attribute_combination_record_t(key, value, &rec)) {
        igraph_attribute_combination_destroy(result);
        return 1;
      }
      igraph_attribute_combination_add(result, rec.name, rec.type, rec.func);
#ifdef IGRAPH_PYTHON3
      free((char*)rec.name);   /* was allocated in pair_to_attribute_combination_record_t above */
#endif
    }
  } else {
    /* assume it is a string or callable */
    if (igraphmodule_i_PyObject_pair_to_attribute_combination_record_t(Py_None, object, &rec)) {
      igraph_attribute_combination_destroy(result);
      return 1;
    }

    igraph_attribute_combination_add(result, 0, rec.type, rec.func);
#ifdef IGRAPH_PYTHON3
    free((char*)rec.name);   /* was allocated in pair_to_attribute_combination_record_t above */
#endif
  }

  return 0;
}

