/* -*- mode: C -*-  */
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

#include "py2compat.h"
#include "random.h"
#include <limits.h>
#include <igraph_random.h>

/**
 * \ingroup python_interface_rng
 * \brief Internal data structure for storing references to the
 *        functions used from Python's random number generator.
 */
typedef struct {
  PyObject* randint_func;
  PyObject* random_func;
  PyObject* gauss_func;
} igraph_i_rng_Python_state_t;

static igraph_i_rng_Python_state_t igraph_rng_Python_state = {0, 0, 0};
static igraph_rng_t igraph_rng_Python = {0, 0, 0};

int igraph_rng_Python_init(void **state) {
  IGRAPH_ERROR("Python RNG error, unsupported function called",
      IGRAPH_EINTERNAL);
  return 0;
}

void igraph_rng_Python_destroy(void *state) {
  igraph_error("Python RNG error, unsupported function called",
      __FILE__, __LINE__, IGRAPH_EINTERNAL);
}

/**
 * \ingroup python_interface_rng
 * \brief Sets the random number generator used by igraph.
 */
PyObject* igraph_rng_Python_set_generator(PyObject* self, PyObject* object) {
  igraph_i_rng_Python_state_t new_state, old_state;
  PyObject* func;

  if (object == Py_None) {
    /* Reverting to the default igraph random number generator instead
     * of the Python-based one */
    igraph_rng_set_default(igraph_rng_default());
    Py_RETURN_NONE;
  }

#define GET_FUNC(name) {\
  func = PyObject_GetAttrString(object, name); \
  if (func == 0) \
    return NULL; \
  if (!PyCallable_Check(func)) {\
    PyErr_SetString(PyExc_TypeError, name "attribute must be callable"); \
    return NULL; \
  } \
}

  GET_FUNC("randint"); new_state.randint_func = func;
  GET_FUNC("random"); new_state.random_func = func;
  GET_FUNC("gauss"); new_state.gauss_func = func;

  old_state = igraph_rng_Python_state;
  igraph_rng_Python_state = new_state;
  Py_XDECREF(old_state.randint_func);
  Py_XDECREF(old_state.random_func);
  Py_XDECREF(old_state.gauss_func);

  igraph_rng_set_default(&igraph_rng_Python);

  Py_RETURN_NONE;
}

/**
 * \ingroup python_interface_rng
 * \brief Sets the seed of the random generator.
 */
int igraph_rng_Python_seed(void *state, unsigned long int seed) {
  IGRAPH_ERROR("Python RNG error, unsupported function called",
      IGRAPH_EINTERNAL);
  return 0;
}

/**
 * \ingroup python_interface_rng
 * \brief Generates an unsigned long integer using the Python random number generator.
 */
unsigned long int igraph_rng_Python_get(void *state) {
  PyObject* result = PyObject_CallFunction(igraph_rng_Python_state.randint_func, "kk", 0, LONG_MAX);
  unsigned long int retval;

  if (result == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    /* Fallback to the C random generator */
    return rand() * LONG_MAX;
  }
  retval = PyInt_AsLong(result);
  Py_DECREF(result);
  return retval;
}

/**
 * \ingroup python_interface_rng
 * \brief Generates a real number between 0 and 1 using the Python random number generator.
 */
igraph_real_t igraph_rng_Python_get_real(void *state) {
  PyObject* result = PyObject_CallFunction(igraph_rng_Python_state.random_func, NULL);
  double retval;

  if (result == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    /* Fallback to the C random generator */
    return rand();
  }

  retval = PyFloat_AsDouble(result);
  Py_DECREF(result);
  return retval;
}

/**
 * \ingroup python_interface_rng
 * \brief Generates a real number distributed according to the normal distribution
 *        around zero with unit variance.
 */
igraph_real_t igraph_rng_Python_get_norm(void *state) {
  PyObject* result = PyObject_CallFunction(igraph_rng_Python_state.gauss_func, "dd", 0.0, 1.0);
  double retval;

  if (result == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    /* Fallback to the C random generator */
    return 0;
  }

  retval = PyFloat_AsDouble(result);
  Py_DECREF(result);
  return retval;
}

/**
 * \ingroup python_interface_rng
 * \brief Specification table for Python's random number generator.
 *        This tells igraph which functions to call to obtain random numbers.
 */
igraph_rng_type_t igraph_rngtype_Python = {
  /* name= */      "Python random generator",
  /* min=  */      0,
  /* max=  */      LONG_MAX,
  /* init= */      igraph_rng_Python_init,
  /* destroy= */   igraph_rng_Python_destroy,
  /* seed= */      igraph_rng_Python_seed,
  /* get= */       igraph_rng_Python_get,
  /* get_real */   igraph_rng_Python_get_real,
  /* get_norm= */  igraph_rng_Python_get_norm,
  /* get_geom= */  0,
  /* get_binom= */ 0
};

void igraphmodule_init_rng(PyObject* igraph_module) {
  PyObject* random_module;

  if (igraph_rng_Python.state != 0)
    return;

  random_module = PyImport_ImportModule("random");
  if (random_module == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    return;
  }

  igraph_rng_Python.type = &igraph_rngtype_Python;
  igraph_rng_Python.state = &igraph_rng_Python_state;

  if (igraph_rng_Python_set_generator(igraph_module, random_module) == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    return;
  }
  Py_DECREF(random_module);
}
