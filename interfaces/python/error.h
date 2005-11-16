#ifndef PYTHON_ERROR_H
#define PYTHON_ERROR_H

#include <Python.h>

/** \defgroup python_interface_errors Error handling
 * \ingroup python_interface */

/** \ingroup python_interface_errors
 * \brief Exception type to be returned when an internal \c igraph error occurs.
 */
PyObject* igraphmodule_InternalError;

PyObject* igraphmodule_handle_igraph_error();
void igraphmodule_igraph_error_hook(const char *reason, const char *file,
				    int line, int igraph_errno);
#endif
