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

#ifndef PYTHON_ERROR_H
#define PYTHON_ERROR_H

#include <Python.h>
#include <igraph_error.h>

/** \defgroup python_interface_errors Error handling
 * \ingroup python_interface */

PyObject* igraphmodule_handle_igraph_error(void);
void igraphmodule_igraph_warning_hook(const char *reason, const char *file,
				    int line, int igraph_errno);
void igraphmodule_igraph_error_hook(const char *reason, const char *file,
				    int line, int igraph_errno);

extern PyObject* igraphmodule_InternalError;

#define IGRAPH_PYCHECK(a) do { \
	int igraph_i_pyret=(a); \
	if (IGRAPH_UNLIKELY(igraph_i_pyret != 0)) {\
		igraphmodule_handle_igraph_error();  \
		IGRAPH_FINALLY_FREE();               \
		return 0;                            \
	} } while (0)

#endif
