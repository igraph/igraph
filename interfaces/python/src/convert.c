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

/************************ Miscellaneous functions *************************/

#include "graphobject.h"
#include "vertexseqobject.h"
#include "vertexobject.h"
#include "edgeseqobject.h"
#include "edgeobject.h"
#include "convert.h"
#include "error.h"
#include "memory.h"

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
    
    if (o == 0 || o==Py_None) return 0;
    if (PyInt_Check(o)) { *result = (int)PyInt_AsLong(o); return 0; }
    if (PyLong_Check(o)) { *result = (int)PyLong_AsLong(o); return 0; }
    if (!PyString_Check(o)) {
        PyErr_SetString(PyExc_TypeError, "int, long or string expected");
        return -1;
    }
    s=PyString_AsString(o);
    /* Convert string to lowercase */
    for (s2=s; *s2; s2++) *s2 = tolower(*s2);
    best = 0; best_unique = 0; best_result = -1;
    /* Search for matches */
    while (table->name != 0) {
        if (strcmp(s, table->name) == 0) { *result = table->value; return 0; }
        for (i=0; s[i] == table->name[i]; i++);
        if (i > best) {
            best = i; best_unique = 1; best_result = table->value;
        } else if (i == best) best_unique = 0;
        table++;
    }
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
        {"infinity", IGRAPH_VCONN_NEI_INFINITY},
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
 * \brief Converts a Python object to an igraph \c igraph_degseq_t
 */
int igraphmodule_PyObject_to_degseq_t(PyObject *o,
  igraph_degseq_t *result) {
  static igraphmodule_enum_translation_table_entry_t degseq_tt[] = {
        {"simple", IGRAPH_DEGSEQ_SIMPLE},
        {"vl", IGRAPH_DEGSEQ_VL},
        {"viger-latapy", IGRAPH_DEGSEQ_VL},
        {0,0}
    };

  return igraphmodule_PyObject_to_enum(o, degseq_tt, (int*)result);
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
  if (object == NULL) {
  } else if (PyInt_Check(object)) {
    long l = PyInt_AS_LONG((PyIntObject*)object);
    *v=l;
    return 0;
  } else if (PyLong_Check(object)) {
    double d = PyLong_AsDouble(object);
    *v=d;
    return 0;
  } else if (PyNumber_Check(object)) {
    PyObject *i = PyNumber_Int(object);
    long l;
    if (i == NULL) return 1;
    l = PyInt_AS_LONG((PyIntObject*)i);
    Py_DECREF(i);
    *v = l;
    return 0;
  }
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
  } else if (PyInt_Check(object)) {
    long l = PyInt_AS_LONG((PyIntObject*)object);
    *v=l;
    return 0;
  } else if (PyLong_Check(object)) {
    double d = PyLong_AsDouble(object);
    *v=d;
    return 0;
  } else if (PyFloat_Check(object)) {
    double d = PyFloat_AS_DOUBLE((PyFloatObject*)object);
    *v=d;
    return 0;
  } else if (PyNumber_Check(object)) {
    PyObject *i = PyNumber_Int(object);
    long l;
    if (i == NULL) return 1;
    l = PyInt_AS_LONG((PyIntObject*)i);
    Py_DECREF(i);
    *v = l;
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
  int i, j, k, ok;
  long idx=0, idx2=0;

  if (PyString_Check(list) || PyUnicode_Check(list)) {
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
  int i, j, k, ok;

  if (PyString_Check(list) || PyUnicode_Check(list)) {
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
  int i, j;

  if (PyString_Check(list) || PyUnicode_Check(list)) {
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
PyObject* igraphmodule_vector_bool_t_to_PyList(igraph_vector_bool_t *v) {
  PyObject *list, *item;
  int n, i;
   
  n=igraph_vector_bool_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();

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
PyObject* igraphmodule_vector_t_to_PyList(igraph_vector_t *v,
    igraphmodule_conv_t type) {
  PyObject *list, *item;
  int n, i;
   
  n=igraph_vector_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();

  list=PyList_New(n);
  for (i=0; i<n; i++) {
    if (type == IGRAPHMODULE_TYPE_INT) {
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
 * \brief Converts an igraph \c igraph_vector_t to a Python integer tuple
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer tuple as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyTuple(igraph_vector_t *v) {
  PyObject* tuple;
  int n, i;
  
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
PyObject* igraphmodule_vector_t_to_PyList_pairs(igraph_vector_t *v) {
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
    char *name = PyString_AsString(o);
    et = (attr_type == ATTRIBUTE_TYPE_VERTEX) ?
      IGRAPH_ATTRIBUTE_VERTEX : IGRAPH_ATTRIBUTE_EDGE;
    if (igraphmodule_i_attribute_get_type(&self->g, &at, et, name)) {
      /* exception was set by igraphmodule_i_attribute_get_type */
      return 1;
    }
    if (at != IGRAPH_ATTRIBUTE_NUMERIC) {
      PyErr_SetString(PyExc_ValueError, "attribute values must be numeric");
      return 1;
    }
    /* Now that the attribute type has been checked, allocate the target
     * vector */
    result = (igraph_vector_t*)calloc(1, sizeof(igraph_vector_t));
    if (result==0) {
      PyErr_NoMemory();
      return 1;
    }
    igraph_vector_init(result, 1);
    if (attr_type == ATTRIBUTE_TYPE_VERTEX) {
      if (igraphmodule_i_get_numeric_vertex_attr(&self->g, name,
          igraph_vss_all(), result)) {
        /* exception has already been set, so return */
		igraph_vector_destroy(result);
        free(result);
        return 1;
      }
    } else {
      if (igraphmodule_i_get_numeric_edge_attr(&self->g, name,
          igraph_ess_all(IGRAPH_EDGEORDER_ID), result)) {
        /* exception has already been set, so return */
		igraph_vector_destroy(result);
        free(result);
        return 1;
      }
    }
    *vptr = result;
  } else if (PySequence_Check(o)) {
    result = (igraph_vector_t*)calloc(1, sizeof(igraph_vector_t));
    if (result==0) {
      PyErr_NoMemory();
      return 1;
    }
    if (igraphmodule_PyObject_to_vector_t(o, result, 0, 0)) {
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
    /* Check whether the attribute exists and is numeric */
    igraph_attribute_type_t at;
    igraph_attribute_elemtype_t et;
    igraph_vector_t tmpvec;
	long i, n;
    char *name = PyString_AsString(o);
    et = (attr_type == ATTRIBUTE_TYPE_VERTEX) ?
      IGRAPH_ATTRIBUTE_VERTEX : IGRAPH_ATTRIBUTE_EDGE;
    if (igraphmodule_i_attribute_get_type(&self->g, &at, et, name)) {
      /* exception was set by igraphmodule_i_attribute_get_type */
      return 1;
    }
    if (at != IGRAPH_ATTRIBUTE_NUMERIC) {
      PyErr_SetString(PyExc_ValueError, "attribute values must be numeric");
      return 1;
    }
    /* Now that the attribute type has been checked, allocate the target
     * vector */
    result = (igraph_vector_bool_t*)calloc(1, sizeof(igraph_vector_bool_t));
    if (result==0) {
      PyErr_NoMemory();
      return 1;
    }
    igraph_vector_init(&tmpvec, 1);
    if (attr_type == ATTRIBUTE_TYPE_VERTEX) {
      if (igraphmodule_i_get_numeric_vertex_attr(&self->g, name,
          igraph_vss_all(), &tmpvec)) {
        /* exception has already been set, so return */
        free(result);
        return 1;
      }
    } else {
      if (igraphmodule_i_get_numeric_edge_attr(&self->g, name,
          igraph_ess_all(IGRAPH_EDGEORDER_ID), &tmpvec)) {
        /* exception has already been set, so return */
        free(result);
        return 1;
      }
    }
    /* Now tmpvec contains the original numeric attributes and we must convert
     * everything to boolean */
    n = igraph_vector_size(&tmpvec);
    if (igraph_vector_bool_init(result, n)) {
      PyErr_NoMemory();
      igraph_vector_destroy(&tmpvec);
      free(result);
    }
	for (i=0; i<n; i++) VECTOR(*result)[i] = (VECTOR(tmpvec)[i] != 0);
    *vptr = result;
    igraph_vector_destroy(&tmpvec);
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
PyObject* igraphmodule_vector_t_pair_to_PyList(igraph_vector_t *v1,
                           igraph_vector_t *v2) {
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
PyObject* igraphmodule_matrix_t_to_PyList(igraph_matrix_t *m,
                         igraphmodule_conv_t type) {
   PyObject *list, *row, *item;
   int nr, nc, i, j;
   
   
   nr=igraph_matrix_nrow(m); nc=igraph_matrix_ncol(m);
   if (nr<0 || nc<0) return igraphmodule_handle_igraph_error();

   // create a new Python list
   list=PyList_New(nr);
   // populate the list with data
   for (i=0; i<nr; i++) 
     {
    row=PyList_New(nc);
    for (j=0; j<nc; j++) 
      {
         if (type==IGRAPHMODULE_TYPE_INT)
           item=PyInt_FromLong(MATRIX(*m,i,j));
         else
           item=PyFloat_FromDouble(MATRIX(*m,i,j));
           
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
PyObject* igraphmodule_vector_ptr_t_to_PyList(igraph_vector_ptr_t *v,
    igraphmodule_conv_t type) {
  PyObject *list, *item;
  int n, i;
   
  n=igraph_vector_ptr_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();

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
  int nr, nc, n, i, j;
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
  int n, i;
  char* ptr;
  
  n=igraph_strvector_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();
  
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
  int n, i;
  static char* empty = "";
  
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
    if (PyString_Check(item))
      ptr=PyString_AS_STRING(item);
    else
      ptr=empty;
    if (igraph_strvector_set(result, i, ptr)) {
      igraph_strvector_destroy(result);
      return 1;
    }
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
 * \brief Tries to interpret a Python object as a vertex selector
 * 
 * \param o the Python object
 * \param vs the \c igraph_vs_t which will contain the result
 * \param return_single will be 1 if the selector selected only a single vertex,
 * 0 otherwise
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vs_t(PyObject *o, igraph_vs_t *vs,
                  igraph_bool_t *return_single) {
  if (return_single) *return_single=0;
  if (o==NULL || o == Py_None) {
    /* Returns a vertex sequence for all vertices */
    igraph_vs_all(vs);
  } else if (PyInt_Check(o)) {
    /* Returns a vertex sequence for a single vertex ID */
    igraph_vs_1(vs, PyInt_AsLong(o));
    if (return_single) *return_single=1;
  } else if (PyLong_Check(o)) {
    /* Returns a vertex sequence for a single vertex ID */
    igraph_vs_1(vs, PyLong_AsLong(o));
    if (return_single) *return_single=1;
  } else if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_VertexSeqType)) {
    igraphmodule_VertexSeqObject *vso = (igraphmodule_VertexSeqObject*)o;
    if (igraph_vs_copy(vs, &vso->vs)) {
      igraphmodule_handle_igraph_error();
      return 1;
    }
  } else if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_VertexType)) {
    igraphmodule_VertexObject *vo = (igraphmodule_VertexObject*)o;
    igraph_vs_1(vs, igraphmodule_Vertex_get_index_long(vo));
    if (return_single) *return_single=1;
  } else {
    /* Returns a vertex sequence with the IDs returned by the iterator */
    PyObject *iterator = PyObject_GetIter(o);
    PyObject *item;
    igraph_vector_t vector;

    if (iterator == NULL) {
      PyErr_SetString(PyExc_TypeError, "integer, long, iterable, Vertex, VertexSeq or None expected");
      return 1;
    }

    IGRAPH_CHECK(igraph_vector_init(&vector, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &vector);
    IGRAPH_CHECK(igraph_vector_reserve(&vector, 20));

    while ((item = PyIter_Next(iterator))) {
      long val=-1;
      if (PyInt_Check(item)) val=PyInt_AsLong(item);
      else if (PyLong_Check(item)) val=PyLong_AsLong(item);
      Py_DECREF(item);

      if (val >= 0) igraph_vector_push_back(&vector, val);
      else {
        PyErr_SetString(PyExc_TypeError, "integer or long expected");
      }

      if (PyErr_Occurred()) break;
    }
    Py_DECREF(iterator);

    if (PyErr_Occurred()) {
      igraph_vector_destroy(&vector);
      IGRAPH_FINALLY_CLEAN(1);
      return 1;
    } else {
      igraph_vs_vector_copy(vs, &vector);
      igraph_vector_destroy(&vector);
      IGRAPH_FINALLY_CLEAN(1);
    }
  }
  return 0;
}


/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as an edge selector
 * 
 * \param o the Python object
 * \param vs the \c igraph_es_t which will contain the result
 * \param return_single will be 1 if the selector selected only a single edge,
 * 0 otherwise
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_es_t(PyObject *o, igraph_es_t *es,
                  igraph_bool_t *return_single) {
  if (return_single) *return_single=0;
  if (o==NULL || o == Py_None) {
    /* Returns a vertex sequence for all edges */
    igraph_es_all(es, IGRAPH_EDGEORDER_ID);
  } else if (PyInt_Check(o)) {
    /* Returns an edge sequence for a single edge ID */
    igraph_es_1(es, PyInt_AsLong(o));
    if (return_single) *return_single=1;
  } else if (PyLong_Check(o)) {
    /* Returns an edge sequence for a single edge ID */
    igraph_es_1(es, PyLong_AsLong(o));
    if (return_single) *return_single=1;
  } else if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_EdgeSeqType)) {
    igraphmodule_EdgeSeqObject *eso = (igraphmodule_EdgeSeqObject*)o;
    if (igraph_es_copy(es, &eso->es)) {
      igraphmodule_handle_igraph_error();
      return 1;
    }
  } else if (PyObject_IsInstance(o, (PyObject*)&igraphmodule_EdgeType)) {
    igraphmodule_EdgeObject *eo = (igraphmodule_EdgeObject*)o;
    igraph_es_1(es, igraphmodule_Edge_get_index_long(eo));
    if (return_single) *return_single=1;
  } else {
    /* Returns an edge sequence with the IDs returned by the iterator */
    PyObject *iterator = PyObject_GetIter(o);
    PyObject *item;
    igraph_vector_t vector;
	igraph_vector_t tuples;

    if (iterator == NULL) {
      PyErr_SetString(PyExc_TypeError, "integer, long, iterable, Edge, EdgeSeq or None expected");
      return 1;
    }

    IGRAPH_CHECK(igraph_vector_init(&vector, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &vector);
    IGRAPH_CHECK(igraph_vector_init(&tuples, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &tuples);

    while ((item = PyIter_Next(iterator))) {
      long val=-1;
      if (PyInt_Check(item)) val=PyInt_AsLong(item);
      else if (PyLong_Check(item)) val=PyLong_AsLong(item);
	  else if (PyTuple_Check(item) && PyTuple_Size(item) >= 2) {
		PyObject *o1, *o2;
		o1=PyTuple_GetItem(item, 0);
		o2=PyTuple_GetItem(item, 1);
		if (PyInt_Check(o1) && PyInt_Check(o2)) {
		  if (igraph_vector_push_back(&tuples, PyInt_AsLong(o1)))
			PyErr_NoMemory();
		  else if (igraph_vector_push_back(&tuples, PyInt_AsLong(o2)))
			PyErr_NoMemory();
		}
		val=-2;
	  }

      Py_DECREF(item);

      if (val >= 0) {
		if (igraph_vector_push_back(&vector, val)) PyErr_NoMemory();
	  } else if (val == -1) {
        PyErr_SetString(PyExc_TypeError, "integer, long or integer tuple expected");
      }

      if (PyErr_Occurred()) break;
    }
    Py_DECREF(iterator);

    if (PyErr_Occurred()) {
      igraph_vector_destroy(&vector);
      igraph_vector_destroy(&tuples);
      IGRAPH_FINALLY_CLEAN(2);
      return 1;
    } else {
      if (igraph_vector_size(&vector) > 0 &&
		  igraph_vector_size(&tuples) == 0)
		igraph_es_vector_copy(es, &vector);
	  else if (igraph_vector_size(&tuples) > 0 &&
		  igraph_vector_size(&vector) == 0)
		igraph_es_pairs(es, &tuples, 1);
	  else if (igraph_vector_size(&tuples) == 0 &&
		  igraph_vector_size(&vector) == 0)
		igraph_es_none(es);
	  else
		PyErr_SetString(PyExc_TypeError, "edge IDs and from-to tuples can not be mixed");
      igraph_vector_destroy(&vector);
      igraph_vector_destroy(&tuples);
      IGRAPH_FINALLY_CLEAN(2);
    }
  }
  if (PyErr_Occurred()) return 1;
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
	char* s=PyString_AsString(obj);
	igraph_layout_drl_default_t def=IGRAPH_LAYOUT_DRL_DEFAULT;
	if (strcmp(s, "default") == 0) def=IGRAPH_LAYOUT_DRL_DEFAULT;
	else if (strcmp(s, "coarsen") == 0) def=IGRAPH_LAYOUT_DRL_COARSEN;
	else if (strcmp(s, "coarsest") == 0) def=IGRAPH_LAYOUT_DRL_COARSEST;
	else if (strcmp(s, "refine") == 0) def=IGRAPH_LAYOUT_DRL_REFINE;
	else if (strcmp(s, "final") == 0) def=IGRAPH_LAYOUT_DRL_FINAL;
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

    printf("cooldown_iterations = %.4f\n", options->cooldown_iterations);
    printf("cooldown_damping_mult = %.4f\n", options->cooldown_damping_mult);
    PyErr_Clear();
	return 0;
  }
  return 0;
}

