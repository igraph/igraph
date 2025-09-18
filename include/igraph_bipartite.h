/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_BIPARTITE_H
#define IGRAPH_BIPARTITE_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_graphicality.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Bipartite networks                                 */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_full_bipartite(igraph_t *graph,
                                        igraph_vector_bool_t *types,
                                        igraph_int_t n1, igraph_int_t n2,
                                        igraph_bool_t directed,
                                        igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_create_bipartite(igraph_t *g, const igraph_vector_bool_t *types,
                                          const igraph_vector_int_t *edges,
                                          igraph_bool_t directed);

IGRAPH_EXPORT igraph_error_t igraph_bipartite_projection_size(const igraph_t *graph,
                                                   const igraph_vector_bool_t *types,
                                                   igraph_int_t *vcount1,
                                                   igraph_int_t *ecount1,
                                                   igraph_int_t *vcount2,
                                                   igraph_int_t *ecount2);

IGRAPH_EXPORT igraph_error_t igraph_bipartite_projection(const igraph_t *graph,
                                              const igraph_vector_bool_t *types,
                                              igraph_t *proj1,
                                              igraph_t *proj2,
                                              igraph_vector_int_t *multiplicity1,
                                              igraph_vector_int_t *multiplicity2,
                                              igraph_int_t probe1);

IGRAPH_EXPORT igraph_error_t igraph_biadjacency(
    igraph_t *graph,
    igraph_vector_bool_t *types,
    const igraph_matrix_t *biadjmatrix,
    igraph_bool_t directed,
    igraph_neimode_t mode,
    igraph_bool_t multiple);

IGRAPH_EXPORT igraph_error_t igraph_weighted_biadjacency(
    igraph_t *graph,
    igraph_vector_bool_t *types,
    igraph_vector_t *weights,
    const igraph_matrix_t *biadjmatrix,
    igraph_bool_t directed,
    igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_get_biadjacency(const igraph_t *graph,
                                                    const igraph_vector_bool_t *types,
                                                    const igraph_vector_t *weights,
                                                    igraph_matrix_t *res,
                                                    igraph_vector_int_t *row_ids,
                                                    igraph_vector_int_t *col_ids);

IGRAPH_EXPORT igraph_error_t igraph_is_bipartite(const igraph_t *graph,
                                      igraph_bool_t *res,
                                      igraph_vector_bool_t *types);

IGRAPH_EXPORT igraph_error_t igraph_bipartite_game_gnp(
        igraph_t *graph,
        igraph_vector_bool_t *types,
        igraph_int_t n1, igraph_int_t n2, igraph_real_t p,
        igraph_bool_t directed, igraph_neimode_t mode,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_bool_t edge_labeled);

IGRAPH_EXPORT igraph_error_t igraph_bipartite_game_gnm(
        igraph_t *graph,
        igraph_vector_bool_t *types,
        igraph_int_t n1, igraph_int_t n2, igraph_int_t m,
        igraph_bool_t directed, igraph_neimode_t mode,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_bool_t edge_labeled);

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_bipartite_iea_game(
    igraph_t *graph, igraph_vector_bool_t *types,
    igraph_int_t n1, igraph_int_t n2, igraph_int_t m,
    igraph_bool_t directed, igraph_neimode_t mode);

IGRAPH_END_C_DECLS

#endif
