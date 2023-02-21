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

#ifndef IGRAPH_CONVERSION_H
#define IGRAPH_CONVERSION_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_matrix.h"
#include "igraph_sparsemat.h"
#include "igraph_attributes.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Conversion                                         */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_get_adjacency(
   const igraph_t *graph, igraph_matrix_t *res, igraph_get_adjacency_t type,
   const igraph_vector_t *weights, igraph_loops_t loops
);
IGRAPH_EXPORT igraph_error_t igraph_get_adjacency_sparse(
   const igraph_t *graph, igraph_sparsemat_t *res, igraph_get_adjacency_t type,
   const igraph_vector_t *weights, igraph_loops_t loops
);

IGRAPH_EXPORT igraph_error_t igraph_get_stochastic(
   const igraph_t *graph, igraph_matrix_t *matrix, igraph_bool_t column_wise,
   const igraph_vector_t *weights
);

IGRAPH_EXPORT igraph_error_t igraph_get_stochastic_sparse(
   const igraph_t *graph, igraph_sparsemat_t *res, igraph_bool_t column_wise,
   const igraph_vector_t *weights
);

/* Deprecated, will be removed in 0.11. Use igraph_get_adjacency_sparse() instead, paying attention to differences. */
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_get_sparsemat(const igraph_t *graph, igraph_sparsemat_t *res);

/* Deprecated, will be removed in 0.11. Use igraph_get_stochastic_sparse() instead, paying attention to differences. */
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_get_stochastic_sparsemat(const igraph_t *graph,
                                                  igraph_sparsemat_t *res,
                                                  igraph_bool_t column_wise);

IGRAPH_EXPORT igraph_error_t igraph_get_edgelist(const igraph_t *graph, igraph_vector_int_t *res, igraph_bool_t bycol);

IGRAPH_EXPORT igraph_error_t igraph_to_directed(igraph_t *graph,
                                     igraph_to_directed_t flags);
IGRAPH_EXPORT igraph_error_t igraph_to_undirected(igraph_t *graph,
                                       igraph_to_undirected_t mode,
                                       const igraph_attribute_combination_t *edge_comb);
IGRAPH_EXPORT igraph_error_t igraph_to_prufer(const igraph_t *graph, igraph_vector_int_t *prufer);

__END_DECLS

#endif
