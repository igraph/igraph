/* vim:set ts=2 sw=2 sts=2 et: */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef PY_IGRAPH_ATTRIBUTES_H
#define PY_IGRAPH_ATTRIBUTES_H

#include <Python.h>
#include <igraph/igraph_attributes.h>
#include <igraph/igraph_datatype.h>
#include <igraph/igraph_iterators.h>
#include <igraph/igraph_strvector.h>
#include <igraph/igraph_vector.h>

#define ATTRHASH_IDX_GRAPH  0
#define ATTRHASH_IDX_VERTEX 1
#define ATTRHASH_IDX_EDGE   2

typedef struct {
  PyObject* attrs[3];
  PyObject* vertex_name_index;
} igraphmodule_i_attribute_struct;

#define ATTR_STRUCT(graph) ((igraphmodule_i_attribute_struct*)((graph)->attr))
#define ATTR_STRUCT_DICT(graph) ((igraphmodule_i_attribute_struct*)((graph)->attr))->attrs

int igraphmodule_i_attribute_get_type(const igraph_t *graph,
				      igraph_attribute_type_t *type,
				      igraph_attribute_elemtype_t elemtype,
				      const char *name);
int igraphmodule_i_get_numeric_graph_attr(const igraph_t *graph,
					  const char *name, igraph_vector_t *value);
int igraphmodule_i_get_numeric_vertex_attr(const igraph_t *graph,
					   const char *name,
					   igraph_vs_t vs,
					   igraph_vector_t *value);
int igraphmodule_i_get_numeric_edge_attr(const igraph_t *graph,
					 const char *name,
					 igraph_es_t es,
					 igraph_vector_t *value);
int igraphmodule_i_get_string_graph_attr(const igraph_t *graph,
					 const char *name, igraph_strvector_t *value);
int igraphmodule_i_get_string_vertex_attr(const igraph_t *graph,
					  const char *name,
					  igraph_vs_t vs,
					  igraph_strvector_t *value);
int igraphmodule_i_get_string_edge_attr(const igraph_t *graph,
					const char *name,
					igraph_es_t es,
					igraph_strvector_t *value);

void igraphmodule_initialize_attribute_handler(void);
void igraphmodule_invalidate_vertex_name_index(igraph_t *graph);
int igraphmodule_get_vertex_id_by_name(igraph_t *graph, PyObject* o, long int* id);

PyObject* igraphmodule_get_edge_attribute_values(const igraph_t* graph,
    const char* name);

igraph_bool_t igraphmodule_has_graph_attribute(const igraph_t *graph, const char* name);
igraph_bool_t igraphmodule_has_vertex_attribute(const igraph_t *graph, const char* name);
igraph_bool_t igraphmodule_has_edge_attribute(const igraph_t *graph, const char* name);

#endif

