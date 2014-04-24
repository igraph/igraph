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

#ifndef PYTHON_FILEHANDLE_H
#define PYTHON_FILEHANDLE_H

#include <Python.h>
#include <stdio.h>

/**
 * \defgroup python_interface_filehandle File handle object
 */

/**
 * \ingroup python_interface_filehandle
 * \brief A structure encapsulating a Python object and a \c FILE* pointer
 * created out of it.
 */
typedef struct {
    PyObject* object;
    FILE* fp;
    unsigned short int need_close;
} igraphmodule_filehandle_t;


int igraphmodule_filehandle_init(igraphmodule_filehandle_t* handle,
        PyObject* object, char* mode);
FILE* igraphmodule_filehandle_get(const igraphmodule_filehandle_t* handle);
void igraphmodule_filehandle_destroy(igraphmodule_filehandle_t* handle);

#endif
