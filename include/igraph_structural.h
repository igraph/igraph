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

#ifndef IGRAPH_STRUCTURAL_H
#define IGRAPH_STRUCTURAL_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_attributes.h"
#include "igraph_sparsemat.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

DECLDIR int igraph_are_connected(const igraph_t *graph, igraph_integer_t v1, igraph_integer_t v2, igraph_bool_t *res);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

DECLDIR int igraph_minimum_spanning_tree(const igraph_t *graph, igraph_vector_t *res,
        const igraph_vector_t *weights);
DECLDIR int igraph_minimum_spanning_tree_unweighted(const igraph_t *graph,
        igraph_t *mst);
DECLDIR int igraph_minimum_spanning_tree_prim(const igraph_t *graph, igraph_t *mst,
        const igraph_vector_t *weights);
DECLDIR int igraph_random_spanning_tree(const igraph_t *graph, igraph_vector_t *res,
                                        igraph_integer_t vid);

DECLDIR int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, igraph_real_t vid,
                                igraph_neimode_t mode);
DECLDIR int igraph_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode);
DECLDIR int igraph_subgraph(const igraph_t *graph, igraph_t *res,
                            const igraph_vs_t vids);
DECLDIR int igraph_induced_subgraph_map(const igraph_t *graph, igraph_t *res,
                                        const igraph_vs_t vids,
                                        igraph_subgraph_implementation_t impl,
                                        igraph_vector_t *map,
                                        igraph_vector_t *invmap);
DECLDIR int igraph_induced_subgraph(const igraph_t *graph, igraph_t *res,
                                    const igraph_vs_t vids, igraph_subgraph_implementation_t impl);
DECLDIR int igraph_subgraph_edges(const igraph_t *graph, igraph_t *res,
                                  const igraph_es_t eids, igraph_bool_t delete_vertices);
DECLDIR int igraph_simplify(igraph_t *graph, igraph_bool_t multiple,
                            igraph_bool_t loops,
                            const igraph_attribute_combination_t *edge_comb);
DECLDIR int igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
                               igraph_bool_t ignore_loops,
                               igraph_reciprocity_t mode);

DECLDIR int igraph_maxdegree(const igraph_t *graph, igraph_integer_t *res,
                             igraph_vs_t vids, igraph_neimode_t mode,
                             igraph_bool_t loops);
DECLDIR int igraph_density(const igraph_t *graph, igraph_real_t *res,
                           igraph_bool_t loops);

DECLDIR int igraph_has_loop(const igraph_t *graph, igraph_bool_t *res);
DECLDIR int igraph_is_loop(const igraph_t *graph, igraph_vector_bool_t *res,
                           igraph_es_t es);
DECLDIR int igraph_is_simple(const igraph_t *graph, igraph_bool_t *res);
DECLDIR int igraph_has_multiple(const igraph_t *graph, igraph_bool_t *res);
DECLDIR int igraph_is_multiple(const igraph_t *graph, igraph_vector_bool_t *res,
                               igraph_es_t es);
DECLDIR int igraph_count_multiple(const igraph_t *graph, igraph_vector_t *res, igraph_es_t es);
DECLDIR int igraph_is_tree(const igraph_t *graph, igraph_bool_t *res, igraph_integer_t *root, igraph_neimode_t mode);
DECLDIR int igraph_girth(const igraph_t *graph, igraph_integer_t *girth,
                         igraph_vector_t *circle);
DECLDIR int igraph_add_edge(igraph_t *graph, igraph_integer_t from, igraph_integer_t to);

DECLDIR int igraph_unfold_tree(const igraph_t *graph, igraph_t *tree,
                               igraph_neimode_t mode, const igraph_vector_t *roots,
                               igraph_vector_t *vertex_index);

DECLDIR int igraph_is_mutual(igraph_t *graph, igraph_vector_bool_t *res, igraph_es_t es);

DECLDIR int igraph_maximum_cardinality_search(const igraph_t *graph,
        igraph_vector_t *alpha,
        igraph_vector_t *alpham1);
DECLDIR int igraph_is_chordal(const igraph_t *graph,
                              const igraph_vector_t *alpha,
                              const igraph_vector_t *alpham1,
                              igraph_bool_t *chordal,
                              igraph_vector_t *fill_in,
                              igraph_t *newgraph);
DECLDIR int igraph_avg_nearest_neighbor_degree(const igraph_t *graph,
        igraph_vs_t vids,
        igraph_neimode_t mode,
        igraph_neimode_t neighbor_degree_mode,
        igraph_vector_t *knn,
        igraph_vector_t *knnk,
        const igraph_vector_t *weights);
DECLDIR int igraph_contract_vertices(igraph_t *graph,
                                     const igraph_vector_t *mapping,
                                     const igraph_attribute_combination_t
                                     *vertex_comb);

DECLDIR int igraph_feedback_arc_set(const igraph_t *graph, igraph_vector_t *result,
                                    const igraph_vector_t *weights, igraph_fas_algorithm_t algo);

DECLDIR int igraph_diversity(igraph_t *graph, const igraph_vector_t *weights,
                             igraph_vector_t *res, const igraph_vs_t vs);

/* -------------------------------------------------- */
/* Spectral Properties                                */
/* -------------------------------------------------- */

DECLDIR int igraph_laplacian(const igraph_t *graph, igraph_matrix_t *res,
                             igraph_sparsemat_t *sparseres,
                             igraph_bool_t normalized,
                             const igraph_vector_t *weights);

/* -------------------------------------------------- */
/* Internal functions, may change any time            */
/* -------------------------------------------------- */

int igraph_i_feedback_arc_set_undirected(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights, igraph_vector_t *layering);
int igraph_i_feedback_arc_set_eades(const igraph_t *graph, igraph_vector_t *result,
                                    const igraph_vector_t *weights, igraph_vector_t *layering);

__END_DECLS

#endif
