/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_GRAPH_ATTRIBUTES_H
#define IGRAPH_GRAPH_ATTRIBUTES_H

#include "igraph_attributes.h"
#include "igraph_strvector.h"
#include "igraph_types.h"
#include "igraph_vector_ptr.h"

#define IGRAPH_I_ATTRIBUTE_DESTROY(graph) \
    do {if ((graph)->attr) igraph_i_attribute_destroy(graph);} while(0)
#define IGRAPH_I_ATTRIBUTE_COPY(to,from,ga,va,ea) do { \
        int igraph_i_ret2=0; \
        if ((from)->attr) { \
            IGRAPH_CHECK(igraph_i_ret2=igraph_i_attribute_copy((to),(from),(ga),(va),(ea))); \
        } else { \
            (to)->attr = 0; \
        } \
        if (igraph_i_ret2 != 0) { \
            IGRAPH_ERROR("", igraph_i_ret2); \
        } \
    } while(0)

int igraph_i_attribute_init(igraph_t *graph, void *attr);
void igraph_i_attribute_destroy(igraph_t *graph);
int igraph_i_attribute_copy(igraph_t *to, const igraph_t *from,
                            igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea);
int igraph_i_attribute_add_vertices(igraph_t *graph, long int nv, void *attr);
int igraph_i_attribute_permute_vertices(const igraph_t *graph,
                                        igraph_t *newgraph,
                                        const igraph_vector_t *idx);
int igraph_i_attribute_combine_vertices(const igraph_t *graph,
                                        igraph_t *newgraph,
                                        const igraph_vector_ptr_t *merges,
                                        const igraph_attribute_combination_t *comb);
int igraph_i_attribute_add_edges(igraph_t *graph,
                                 const igraph_vector_t *edges, void *attr);
int igraph_i_attribute_permute_edges(const igraph_t *graph,
                                     igraph_t *newgraph,
                                     const igraph_vector_t *idx);
int igraph_i_attribute_combine_edges(const igraph_t *graph,
                                     igraph_t *newgraph,
                                     const igraph_vector_ptr_t *merges,
                                     const igraph_attribute_combination_t *comb);

int igraph_i_attribute_get_info(const igraph_t *graph,
                                igraph_strvector_t *gnames,
                                igraph_vector_t *gtypes,
                                igraph_strvector_t *vnames,
                                igraph_vector_t *vtypes,
                                igraph_strvector_t *enames,
                                igraph_vector_t *etypes);
igraph_bool_t igraph_i_attribute_has_attr(const igraph_t *graph,
                                          igraph_attribute_elemtype_t type,
                                          const char *name);
int igraph_i_attribute_gettype(const igraph_t *graph,
                               igraph_attribute_type_t *type,
                               igraph_attribute_elemtype_t elemtype,
                               const char *name);

int igraph_i_attribute_get_numeric_graph_attr(const igraph_t *graph,
                                              const char *name,
                                              igraph_vector_t *value);
int igraph_i_attribute_get_numeric_vertex_attr(const igraph_t *graph,
                                               const char *name,
                                               igraph_vs_t vs,
                                               igraph_vector_t *value);
int igraph_i_attribute_get_numeric_edge_attr(const igraph_t *graph,
                                             const char *name,
                                             igraph_es_t es,
                                             igraph_vector_t *value);
int igraph_i_attribute_get_string_graph_attr(const igraph_t *graph,
                                             const char *name,
                                             igraph_strvector_t *value);
int igraph_i_attribute_get_string_vertex_attr(const igraph_t *graph,
                                              const char *name,
                                              igraph_vs_t vs,
                                              igraph_strvector_t *value);
int igraph_i_attribute_get_string_edge_attr(const igraph_t *graph,
                                            const char *name,
                                            igraph_es_t es,
                                            igraph_strvector_t *value);
int igraph_i_attribute_get_bool_graph_attr(const igraph_t *graph,
                                           const char *name,
                                           igraph_vector_bool_t *value);
int igraph_i_attribute_get_bool_vertex_attr(const igraph_t *graph,
                                            const char *name,
                                            igraph_vs_t vs,
                                            igraph_vector_bool_t *value);
int igraph_i_attribute_get_bool_edge_attr(const igraph_t *graph,
                                          const char *name,
                                          igraph_es_t es,
                                          igraph_vector_bool_t *value);

#endif /* IGRAPH_GRAPH_ATTRIBUTES_H */
