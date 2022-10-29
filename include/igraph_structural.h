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
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_constants.h"
#include "igraph_iterators.h"
#include "igraph_matrix.h"
#include "igraph_sparsemat.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_are_connected(const igraph_t *graph, igraph_integer_t v1, igraph_integer_t v2, igraph_bool_t *res);
IGRAPH_EXPORT igraph_error_t igraph_count_multiple(const igraph_t *graph, igraph_vector_int_t *res, igraph_es_t es);
IGRAPH_EXPORT igraph_error_t igraph_count_multiple_1(const igraph_t *graph, igraph_integer_t *res, igraph_integer_t eid);
IGRAPH_EXPORT igraph_error_t igraph_density(const igraph_t *graph, igraph_real_t *res,
                                 igraph_bool_t loops);
IGRAPH_EXPORT igraph_error_t igraph_diversity(const igraph_t *graph, const igraph_vector_t *weights,
                                   igraph_vector_t *res, const igraph_vs_t vs);
IGRAPH_EXPORT igraph_error_t igraph_girth(const igraph_t *graph, igraph_real_t *girth,
                               igraph_vector_int_t *circle);
IGRAPH_EXPORT igraph_error_t igraph_has_loop(const igraph_t *graph, igraph_bool_t *res);
IGRAPH_EXPORT igraph_error_t igraph_has_multiple(const igraph_t *graph, igraph_bool_t *res);
IGRAPH_EXPORT igraph_error_t igraph_is_loop(const igraph_t *graph, igraph_vector_bool_t *res,
                                 igraph_es_t es);
IGRAPH_EXPORT igraph_error_t igraph_is_multiple(const igraph_t *graph, igraph_vector_bool_t *res,
                                     igraph_es_t es);
IGRAPH_EXPORT igraph_error_t igraph_is_mutual(const igraph_t *graph, igraph_vector_bool_t *res, igraph_es_t es, igraph_bool_t loops);
IGRAPH_EXPORT igraph_error_t igraph_has_mutual(const igraph_t *graph, igraph_bool_t *res, igraph_bool_t loops);
IGRAPH_EXPORT igraph_error_t igraph_is_simple(const igraph_t *graph, igraph_bool_t *res);
IGRAPH_EXPORT igraph_error_t igraph_is_tree(const igraph_t *graph, igraph_bool_t *res, igraph_integer_t *root, igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_is_acyclic(const igraph_t *graph, igraph_bool_t *res);
IGRAPH_EXPORT igraph_error_t igraph_is_forest(const igraph_t *graph, igraph_bool_t *res,
                                    igraph_vector_int_t *roots, igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_maxdegree(const igraph_t *graph, igraph_integer_t *res,
                                   igraph_vs_t vids, igraph_neimode_t mode,
                                   igraph_bool_t loops);
IGRAPH_EXPORT igraph_error_t igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
                                     igraph_bool_t ignore_loops,
                                     igraph_reciprocity_t mode);
IGRAPH_EXPORT igraph_error_t igraph_strength(const igraph_t *graph, igraph_vector_t *res,
                                  const igraph_vs_t vids, igraph_neimode_t mode,
                                  igraph_bool_t loops, const igraph_vector_t *weights);
IGRAPH_EXPORT igraph_error_t igraph_sort_vertex_ids_by_degree(const igraph_t *graph,
                                                   igraph_vector_int_t *outvids,
                                                   igraph_vs_t vids,
                                                   igraph_neimode_t mode,
                                                   igraph_bool_t loops,
                                                   igraph_order_t order,
                                                   igraph_bool_t only_indices);
IGRAPH_EXPORT igraph_error_t igraph_is_perfect(const igraph_t *graph, igraph_bool_t *perfect);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_minimum_spanning_tree(const igraph_t *graph, igraph_vector_int_t *res,
                                               const igraph_vector_t *weights);
IGRAPH_EXPORT igraph_error_t igraph_minimum_spanning_tree_unweighted(const igraph_t *graph,
                                                          igraph_t *mst);
IGRAPH_EXPORT igraph_error_t igraph_minimum_spanning_tree_prim(const igraph_t *graph, igraph_t *mst,
                                                    const igraph_vector_t *weights);
IGRAPH_EXPORT igraph_error_t igraph_random_spanning_tree(const igraph_t *graph, igraph_vector_int_t *res,
                                              igraph_integer_t vid);

IGRAPH_EXPORT igraph_error_t igraph_subcomponent(const igraph_t *graph, igraph_vector_int_t *res, igraph_integer_t vid,
                                      igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_unfold_tree(const igraph_t *graph, igraph_t *tree,
                                     igraph_neimode_t mode, const igraph_vector_int_t *roots,
                                     igraph_vector_int_t *vertex_index);

IGRAPH_EXPORT igraph_error_t igraph_maximum_cardinality_search(const igraph_t *graph,
                                                    igraph_vector_int_t *alpha,
                                                    igraph_vector_int_t *alpham1);
IGRAPH_EXPORT igraph_error_t igraph_is_chordal(const igraph_t *graph,
                                    const igraph_vector_int_t *alpha,
                                    const igraph_vector_int_t *alpham1,
                                    igraph_bool_t *chordal,
                                    igraph_vector_int_t *fill_in,
                                    igraph_t *newgraph);
IGRAPH_EXPORT igraph_error_t igraph_avg_nearest_neighbor_degree(const igraph_t *graph,
                                                     igraph_vs_t vids,
                                                     igraph_neimode_t mode,
                                                     igraph_neimode_t neighbor_degree_mode,
                                                     igraph_vector_t *knn,
                                                     igraph_vector_t *knnk,
                                                     const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_feedback_arc_set(const igraph_t *graph, igraph_vector_int_t *result,
                                          const igraph_vector_t *weights, igraph_fas_algorithm_t algo);

/* -------------------------------------------------- */
/* Spectral Properties                                */
/* -------------------------------------------------- */

/**
 * \typedef igraph_laplacian_normalization_t
 * \brief Normalization methods for a Laplacian matrix.
 *
 * Normalization methods for \ref igraph_get_laplacian() and
 * \ref igraph_get_laplacian_sparse(). In the following, \c A refers to the
 * (possibly weighted) adjacency matrix and \c D is a diagonal matrix containing
 * degrees (unweighted case) or strengths (weighted case). Out-, in- or total degrees
 * are used according to the \p mode parameter.
 *
 * \enumval IGRAPH_LAPLACIAN_UNNORMALIZED Unnormalized Laplacian, <code>L = D - A</code>.
 * \enumval IGRAPH_LAPLACIAN_SYMMETRIC Symmetric normalized Laplacian, <code>L = I - D^(-1/2) A D^(-1/2)</code>.
 * \enumval IGRAPH_LAPLACIAN_LEFT Left-stochastic normalized Laplacian, <code>L = I - D^-1 A</code>.
 * \enumval IGRAPH_LAPLACIAN_RIGHT Right-stochastic normalized Laplacian, <code>L = I - A D^-1</code>.
 */
typedef enum {
    IGRAPH_LAPLACIAN_UNNORMALIZED = 0,
    IGRAPH_LAPLACIAN_SYMMETRIC = 1,
    IGRAPH_LAPLACIAN_LEFT = 2,
    IGRAPH_LAPLACIAN_RIGHT = 3
} igraph_laplacian_normalization_t;

IGRAPH_EXPORT igraph_error_t igraph_get_laplacian(
   const igraph_t *graph, igraph_matrix_t *res, igraph_neimode_t mode,
   igraph_laplacian_normalization_t normalization,
   const igraph_vector_t *weights
);
IGRAPH_EXPORT igraph_error_t igraph_get_laplacian_sparse(
   const igraph_t *graph, igraph_sparsemat_t *sparseres, igraph_neimode_t mode,
   igraph_laplacian_normalization_t normalization,
   const igraph_vector_t *weights
);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_laplacian(
   const igraph_t *graph, igraph_matrix_t *res, igraph_sparsemat_t *sparseres,
   igraph_bool_t normalized, const igraph_vector_t *weights
);

__END_DECLS

#endif
