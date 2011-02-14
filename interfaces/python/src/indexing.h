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

#ifndef PYTHON_INDEXING_H
#define PYTHON_INDEXING_H

#include <Python.h>
#include <igraph_datatype.h>

PyObject* igraphmodule_Graph_adjmatrix_get_index(igraph_t* graph,
        PyObject* row_index, PyObject* column_index, PyObject* attr_name);
int igraphmodule_Graph_adjmatrix_set_index(igraph_t* graph,
        PyObject* row_index, PyObject* column_index, PyObject* attr_name,
        PyObject* value);

#endif
