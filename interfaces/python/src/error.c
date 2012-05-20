/* -*- mode: C -*-  */
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

#include "error.h"
#include <igraph.h>

/** \ingroup python_interface_errors
 * \brief Exception type to be returned when an internal \c igraph error occurs.
 */
PyObject* igraphmodule_InternalError;

/**
 * \ingroup python_interface_errors
 * \brief Generic error handler for internal \c igraph errors.
 * 
 * Since now \c igraph supports error handler functions, a special
 * function called \c igraphmodule_igraph_error_hook is responsible
 * for providing a meaningful error message. If it fails (or it isn't
 * even called), this function will provide a default error message.
 * 
 * \return Always returns \c NULL, and all callers are advised to pass this
 * \c NULL value to their callers until it is propagated to the Python
 * interpreter.
 */
PyObject* igraphmodule_handle_igraph_error() 
{
  if (!PyErr_Occurred()) {
    PyErr_SetString(igraphmodule_InternalError,
		    "Internal igraph error. Please contact the author!");
  }

  return NULL;
}

/**
 * \ingroup python_interface_errors
 * \brief Warning hook for \c igraph
 */
void igraphmodule_igraph_warning_hook(const char *reason, const char *file,
				    int line, int igraph_errno) {
  char buf[4096];
  sprintf(buf, "%s at %s:%i", reason, file, line);
  PyErr_Warn(PyExc_RuntimeWarning, buf);
}

/**
 * \ingroup python_interface_errors
 * \brief Error hook for \c igraph
 */
void igraphmodule_igraph_error_hook(const char *reason, const char *file,
				    int line, int igraph_errno) {
  char buf[4096];
  PyObject *exc = igraphmodule_InternalError;

  if (igraph_errno == IGRAPH_UNIMPLEMENTED)
      exc = PyExc_NotImplementedError;

  if (igraph_errno == IGRAPH_ENOMEM)
      exc = PyExc_MemoryError;

  sprintf(buf, "Error at %s:%i: %s, %s", file, line, reason,
	  igraph_strerror(igraph_errno));
  IGRAPH_FINALLY_FREE();

  /* make sure we are not masking already thrown exceptions */
  if (!PyErr_Occurred())
    PyErr_SetString(exc, buf);
}
