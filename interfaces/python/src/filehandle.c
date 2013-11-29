/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010-2012  Tamas Nepusz <ntamas@gmail.com>
   
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

#include "filehandle.h"
#include "py2compat.h"

/**
 * \ingroup python_interface_filehandle
 * \brief Constructs a new file handle object from a Python object.
 *
 * \return 0 if everything was OK, 1 otherwise. An appropriate Python
 *   exception is raised in this case.
 */
int igraphmodule_filehandle_init(igraphmodule_filehandle_t* handle,
        PyObject* object, char* mode) {
    handle->need_close = 0;

#ifdef IGRAPH_PYTHON3
    int fp;
    if (object == 0 || PyLong_Check(object)) {
        PyErr_SetString(PyExc_TypeError, "string or file-like object expected");
        return 1;
    }
#else
    if (object == 0 ||
        (!PyBaseString_Check(object) && !PyFile_Check(object))) {
        PyErr_SetString(PyExc_TypeError, "string or file handle expected");
        return 1;
    }
#endif

    if (PyBaseString_Check(object)) {
        /* We have received a string; we need to open the file denoted by this
         * string now and mark that we opened the file ourselves (so we need
         * to close it when igraphmodule_filehandle_destroy is invoked). */
#ifdef IGRAPH_PYTHON3
        handle->object = PyFile_FromObject(object, mode);
#else
        handle->object = PyFile_FromString(PyString_AsString(object), mode);
#endif
        if (handle->object == 0) {
            /* Could not open the file; just return an error code because an
             * exception was raised already */
            return 1;
        }
        /* Remember that we need to close the file ourselves */
        handle->need_close = 1;
    } else {
        /* This is probably a file-like object; store a reference for it and
         * we will handle it later */
        handle->object = object;
        Py_INCREF(handle->object);
    }

    /* At this stage, handle->object is something we can handle.
     * In Python 2, we get here only if object is a file object so we
     * can safely call PyFile_AsFile to get a FILE* object.
     * In Python 3, we have to call PyObject_AsFileDescriptor instead
     * and then fdopen() it to get the corresponding FILE* object.
     */
#ifdef IGRAPH_PYTHON3
    fp = PyObject_AsFileDescriptor(handle->object);
    if (fp == -1) {
        igraphmodule_filehandle_destroy(handle);
        /* This already called Py_DECREF(handle->object), no need to call it */
        return 1;
    }
    handle->fp = fdopen(fp, mode);
    if (handle->fp == 0) {
        igraphmodule_filehandle_destroy(handle);
        /* This already called Py_DECREF(handle->object), no need to call it */
        PyErr_SetString(PyExc_RuntimeError, "fdopen() failed unexpectedly");
        return 1;
    }
#else
    handle->fp = PyFile_AsFile(handle->object);
    if (handle->fp == 0) {
        igraphmodule_filehandle_destroy(handle);
        /* This already called Py_DECREF(handle->object), no need to call it */
        PyErr_SetString(PyExc_RuntimeError, "PyFile_AsFile() failed unexpectedly");
        return 1;
    }
#endif

    return 0;
}

/**
 * \ingroup python_interface_filehandle
 * \brief Destroys the file handle object.
 */
void igraphmodule_filehandle_destroy(igraphmodule_filehandle_t* handle) {
	if (handle->fp != 0) {
		fflush(handle->fp);
	}
    handle->fp = 0;
    
    if (handle->object != 0) {
        if (handle->need_close) {
            if (PyFile_Close(handle->object)) {
                PyErr_WriteUnraisable(0);
            }
        }
        Py_DECREF(handle->object);
    }

    handle->object = 0;
    handle->need_close = 0;
}

/**
 * \ingroup python_interface_filehandle
 * \brief Returns the file encapsulated by the given \c igraphmodule_filehandle_t.
 */
FILE* igraphmodule_filehandle_get(const igraphmodule_filehandle_t* handle) {
    return handle->fp;
}

