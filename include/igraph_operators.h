/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#ifndef IGRAPH_OPERATORS_H
#define IGRAPH_OPERATORS_H

#include "igraph_decls.h"

#include "igraph_attributes.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Graph operators                                    */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_add_edge(igraph_t *graph, igraph_integer_t from, igraph_integer_t to);
IGRAPH_EXPORT igraph_error_t igraph_disjoint_union(igraph_t *res,
                                        const igraph_t *left, const igraph_t *right);
IGRAPH_EXPORT igraph_error_t igraph_disjoint_union_many(igraph_t *res,
                                             const igraph_vector_ptr_t *graphs);
IGRAPH_EXPORT igraph_error_t igraph_union(igraph_t *res, const igraph_t *left, const igraph_t *right,
                               igraph_vector_int_t *edge_map1, igraph_vector_int_t *edge_map2);
IGRAPH_EXPORT igraph_error_t igraph_union_many(igraph_t *res, const igraph_vector_ptr_t *graphs,
                                    igraph_vector_int_list_t *edgemaps);
IGRAPH_EXPORT igraph_error_t igraph_join(igraph_t *res, const igraph_t *left,
                                         const igraph_t *right);
IGRAPH_EXPORT igraph_error_t igraph_intersection(igraph_t *res,
                                      const igraph_t *left, const igraph_t *right,
                                      igraph_vector_int_t *edge_map1,
                                      igraph_vector_int_t *edge_map2);
IGRAPH_EXPORT igraph_error_t igraph_intersection_many(igraph_t *res,
                                           const igraph_vector_ptr_t *graphs,
                                           igraph_vector_int_list_t *edgemaps);
IGRAPH_EXPORT igraph_error_t igraph_difference(igraph_t *res,
                                    const igraph_t *orig, const igraph_t *sub);
IGRAPH_EXPORT igraph_error_t igraph_complementer(igraph_t *res, const igraph_t *graph,
                                      igraph_bool_t loops);
IGRAPH_EXPORT igraph_error_t igraph_compose(igraph_t *res, const igraph_t *g1, const igraph_t *g2,
                                 igraph_vector_int_t *edge_map1, igraph_vector_int_t *edge_map2);
IGRAPH_EXPORT igraph_error_t igraph_contract_vertices(igraph_t *graph,
                                           const igraph_vector_int_t *mapping,
                                           const igraph_attribute_combination_t *vertex_comb);
IGRAPH_EXPORT igraph_error_t igraph_permute_vertices(const igraph_t *graph, igraph_t *res,
                                          const igraph_vector_int_t *permutation);
IGRAPH_EXPORT igraph_error_t igraph_connect_neighborhood(igraph_t *graph, igraph_integer_t order,
                                              igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_graph_power(const igraph_t *graph, igraph_t *res,
                                                igraph_integer_t order, igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode);
IGRAPH_EXPORT igraph_error_t igraph_simplify(igraph_t *graph, igraph_bool_t multiple,
                                  igraph_bool_t loops,
                                  const igraph_attribute_combination_t *edge_comb);
IGRAPH_EXPORT igraph_error_t igraph_induced_subgraph_map(const igraph_t *graph, igraph_t *res,
                                              const igraph_vs_t vids,
                                              igraph_subgraph_implementation_t impl,
                                              igraph_vector_int_t *map,
                                              igraph_vector_int_t *invmap);
IGRAPH_EXPORT igraph_error_t igraph_induced_subgraph(const igraph_t *graph, igraph_t *res,
                                          const igraph_vs_t vids, igraph_subgraph_implementation_t impl);
IGRAPH_EXPORT igraph_error_t igraph_induced_subgraph_edges(
        const igraph_t *graph, igraph_vs_t vids, igraph_vector_int_t *edges);
IGRAPH_EXPORT igraph_error_t igraph_subgraph_from_edges(const igraph_t *graph, igraph_t *res,
                                        const igraph_es_t eids, igraph_bool_t delete_vertices);
IGRAPH_EXPORT igraph_error_t igraph_reverse_edges(igraph_t *graph, const igraph_es_t eids);

IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_subgraph_edges(
        const igraph_t *graph, igraph_t *res, const igraph_es_t eids,
        igraph_bool_t delete_vertices
);
__END_DECLS

#endif
