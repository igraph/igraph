/* vim:set ts=4 sw=2 sts=2 et:  */
/* 
   IGraph library - Python interface.
   Copyright (C) 2006-2011  Tamas Nepusz <ntamas@gmail.com>
   5 Avenue Road, Staines, Middlesex, TW18 3AW, United Kingdom
   
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

#ifndef PYTHON_HELPERS_H
#define PYTHON_HELPERS_H

#include <Python.h>

PyObject* PyList_NewFill(Py_ssize_t len, PyObject* item);
PyObject* PyList_Zeroes(Py_ssize_t len);

char* PyObject_ConvertToCString(PyObject* string);

#define PY_IGRAPH_DEPRECATED(msg) \
  PyErr_WarnEx(PyExc_DeprecationWarning, (msg), 1)
#define PY_IGRAPH_WARN(msg) \
  PyErr_WarnEx(PyExc_RuntimeWarning, (msg), 1)

#endif
