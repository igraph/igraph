#include "error.h"

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
 * \brief Error hook for \c igraph
 */
void igraphmodule_igraph_error_hook(const char *reason, const char *file,
				    int line, int igraph_errno) {
  char buf[4096];
  sprintf(buf, "Error at %s:%i: %s, %s", file, line, reason,
	  igraph_strerror(igraph_errno));
  PyErr_SetString(igraphmodule_InternalError, buf);
}
