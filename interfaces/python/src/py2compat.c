/* -*- mode: C -*-  */
/* vim: set ts=2 sw=2 sts=2 et: */

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

/* Common utility functions that are useful both in Python 2.x and 3.x */

int PyFile_Close(PyObject* fileObj) {
  PyObject *result;

  result = PyObject_CallMethod(fileObj, "close", 0);
  if (result) {
    Py_DECREF(result);
    return 0;
  } else {
    /* Exception raised already */
    return 1;
  }
}


#ifdef IGRAPH_PYTHON3

/* Python 3.x functions */

PyObject* PyFile_FromObject(PyObject* filename, const char* mode) {
  PyObject *ioModule, *fileObj;

  ioModule = PyImport_ImportModule("io");
  if (ioModule == 0)
    return 0;

  fileObj = PyObject_CallMethod(ioModule, "open", "Os", filename, mode);
  Py_DECREF(ioModule);

  return fileObj;
}

char* PyString_CopyAsString(PyObject* string) {
  PyObject* bytes;
  char* result;

  if (PyBytes_Check(string)) {
    bytes = string;
    Py_INCREF(bytes);
  } else {
    bytes = PyUnicode_AsUTF8String(string);
  }

  if (bytes == 0)
    return 0;
  
  result = strdup(PyBytes_AS_STRING(bytes));
  Py_DECREF(bytes);

  if (result == 0)
    PyErr_NoMemory();

  return result;
}

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

#else

/* Python 2.x functions */

char* PyString_CopyAsString(PyObject* string) {
  char* result;

  if (!PyBaseString_Check(string)) {
    PyErr_SetString(PyExc_TypeError, "string or unicode object expected");
    return 0;
  }

  result = PyString_AsString(string);
  if (result == 0)
    return 0;

  result = strdup(result);
  if (result == 0)
    PyErr_NoMemory();

  return result;
}

int PyString_IsEqualToASCIIString(PyObject* py_string,
		const char* c_string) {
	PyObject* c_string_conv;
	int result;

  if (PyString_Check(py_string)) {
    return strcmp(PyString_AS_STRING(py_string), c_string) == 0;
  }

	if (!PyUnicode_Check(py_string))
		return 0;

	c_string_conv = PyUnicode_DecodeASCII(c_string, strlen(c_string), "strict");
	if (c_string_conv == 0)
		return 0;

	result = (PyUnicode_Compare(py_string, c_string_conv) == 0);
	Py_DECREF(c_string_conv);

	return result;
}

#endif
