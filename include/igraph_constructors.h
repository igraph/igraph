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

#ifndef IGRAPH_CONSTRUCTORS_H
#define IGRAPH_CONSTRUCTORS_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_matrix.h"
#include "igraph_datatype.h"
#include "igraph_graphicality.h"
#include "igraph_sparsemat.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Constructors, deterministic                        */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_create(igraph_t *graph, const igraph_vector_int_t *edges, igraph_integer_t n,
                                igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_small(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                                          int first, ...);
IGRAPH_EXPORT igraph_error_t igraph_adjacency(
        igraph_t *graph, const igraph_matrix_t *adjmatrix, igraph_adjacency_t mode,
        igraph_loops_t loops);
IGRAPH_EXPORT igraph_error_t igraph_weighted_adjacency(
        igraph_t *graph, const igraph_matrix_t *adjmatrix, igraph_adjacency_t mode,
        igraph_vector_t *weights, igraph_loops_t loops);
IGRAPH_EXPORT igraph_error_t igraph_sparse_adjacency(igraph_t *graph, igraph_sparsemat_t *adjmatrix, igraph_adjacency_t mode, igraph_loops_t loops);
IGRAPH_EXPORT igraph_error_t igraph_sparse_weighted_adjacency(igraph_t *graph, igraph_sparsemat_t *adjmatrix, igraph_adjacency_t mode, igraph_vector_t *weights, igraph_loops_t loops);
IGRAPH_EXPORT igraph_error_t igraph_star(igraph_t *graph, igraph_integer_t n, igraph_star_mode_t mode,
                              igraph_integer_t center);
IGRAPH_EXPORT igraph_error_t igraph_wheel(igraph_t *graph, igraph_integer_t n, igraph_wheel_mode_t mode,
                              igraph_integer_t center);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_lattice(igraph_t *graph, const igraph_vector_int_t *dimvector, igraph_integer_t nei,
                                 igraph_bool_t directed, igraph_bool_t mutual, igraph_bool_t circular);
IGRAPH_EXPORT igraph_error_t igraph_square_lattice(igraph_t *graph, const igraph_vector_int_t *dimvector, igraph_integer_t nei,
                                 igraph_bool_t directed, igraph_bool_t mutual, const igraph_vector_bool_t *circular);
IGRAPH_EXPORT igraph_error_t igraph_ring(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                              igraph_bool_t mutual, igraph_bool_t circular);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children,
                              igraph_tree_mode_t type);
IGRAPH_EXPORT igraph_error_t igraph_kary_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children,
                                              igraph_tree_mode_t type);
IGRAPH_EXPORT igraph_error_t igraph_symmetric_tree(igraph_t *graph, const igraph_vector_int_t *branches,
                                                   igraph_tree_mode_t type);
IGRAPH_EXPORT igraph_error_t igraph_regular_tree(igraph_t *graph, igraph_integer_t h, igraph_integer_t k,
                                                 igraph_tree_mode_t type);
IGRAPH_EXPORT igraph_error_t igraph_tree_from_parent_vector(igraph_t *graph, const igraph_vector_int_t *parents,
                                                            igraph_tree_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_from_prufer(igraph_t *graph, const igraph_vector_int_t *prufer);
IGRAPH_EXPORT igraph_error_t igraph_full(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, igraph_bool_t loops);
IGRAPH_EXPORT igraph_error_t igraph_full_multipartite(igraph_t *graph, igraph_vector_int_t *types, const igraph_vector_int_t *n,
                                        igraph_bool_t directed, igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_turan(igraph_t *graph, igraph_vector_int_t *types, igraph_integer_t n, igraph_integer_t r);
IGRAPH_EXPORT igraph_error_t igraph_full_citation(igraph_t *graph, igraph_integer_t n,
                                       igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_atlas(igraph_t *graph, igraph_integer_t number);
IGRAPH_EXPORT igraph_error_t igraph_extended_chordal_ring(igraph_t *graph, igraph_integer_t nodes,
                                               const igraph_matrix_int_t *W, igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_linegraph(const igraph_t *graph, igraph_t *linegraph);

IGRAPH_EXPORT igraph_error_t igraph_de_bruijn(igraph_t *graph, igraph_integer_t m, igraph_integer_t n);
IGRAPH_EXPORT igraph_error_t igraph_circulant(igraph_t *graph, igraph_integer_t n, const igraph_vector_int_t *l, igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_generalized_petersen(igraph_t *graph, igraph_integer_t n, igraph_integer_t k);
IGRAPH_EXPORT igraph_error_t igraph_kautz(igraph_t *graph, igraph_integer_t m, igraph_integer_t n);
IGRAPH_EXPORT igraph_error_t igraph_famous(igraph_t *graph, const char *name);
IGRAPH_EXPORT igraph_error_t igraph_lcf_vector(igraph_t *graph, igraph_integer_t n,
                                    const igraph_vector_int_t *shifts,
                                    igraph_integer_t repeats);
IGRAPH_EXPORT igraph_error_t igraph_lcf(igraph_t *graph, igraph_integer_t n, ...);
IGRAPH_EXPORT igraph_error_t igraph_realize_degree_sequence(igraph_t *graph,
                                                 const igraph_vector_int_t *outdeg, const igraph_vector_int_t *indeg,
                                                 igraph_edge_type_sw_t allowed_edge_types,
                                                 igraph_realize_degseq_t method);
IGRAPH_EXPORT igraph_error_t igraph_triangular_lattice(igraph_t *graph, const igraph_vector_int_t *dims, igraph_bool_t directed, igraph_bool_t mutual);
IGRAPH_EXPORT igraph_error_t igraph_hexagonal_lattice(igraph_t *graph, const igraph_vector_int_t *dims, igraph_bool_t directed, igraph_bool_t mutual);
IGRAPH_EXPORT igraph_error_t igraph_realize_bipartite_degree_sequence(igraph_t *graph, const igraph_vector_int_t *deg1, const igraph_vector_int_t *deg2, const igraph_edge_type_sw_t allowed_edge_types, const igraph_realize_degseq_t method);

__END_DECLS

#endif
