/* -*- mode: C -*-  */
/* vim: set ts=2 sw=2 sts=2 et: */

/* 
   IGraph library.
   Copyright (C) 2006-2010  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include "py2compat.h"

#ifdef IGRAPH_PYTHON3

int PyString_IsEqualToUTF8String(PyObject* py_string,
		const char* c_string) {
	PyObject* c_string_conv;
	int result;

	if (!PyUnicode_Check(py_string))
		return 0;

	c_string_conv = PyUnicode_FromString(c_string);
	if (c_string_conv == 0)
		return 0;

	result = (PyUnicode_Compare(py_string, c_string_conv) == 0);
	Py_DECREF(c_string_conv);

	return result;
}

#endif
