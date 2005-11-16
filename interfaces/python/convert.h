/************************ Miscellaneous functions *************************/

/** \defgroup python_interface_conversion Converting between Python and igraph data types
 * \ingroup python_interface */

#ifndef PYTHON_CONVERT_H
#define PYTHON_CONVERT_H

#include <Python.h>
#include "types.h"

typedef enum { IGRAPHMODULE_TYPE_INT=0, IGRAPHMODULE_TYPE_FLOAT }
igraphmodule_conv_t;

int igraphmodule_PyList_to_vector_t(PyObject *list, vector_t *v, bool_t need_non_negative, bool_t pairs);
PyObject* igraphmodule_vector_t_to_PyList(vector_t *v);
PyObject* igraphmodule_vector_t_to_PyList_pairs(vector_t *v);
PyObject* igraphmodule_vector_t_to_float_PyList(vector_t *v);
PyObject* igraphmodule_matrix_t_to_PyList(matrix_t *m,
						 igraphmodule_conv_t type);
PyObject* igraphmodule_strvector_t_to_PyList(igraph_strvector_t *v);
#endif
