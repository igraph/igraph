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

#ifndef IGRAPH_COMPONENTS_H
#define IGRAPH_COMPONENTS_H

#include "igraph_decls.h"

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_graph_list.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_list.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Components                                         */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_connected_components(const igraph_t *graph, igraph_vector_int_t *membership,
                                  igraph_vector_int_t *csize, igraph_int_t *no,
                                  igraph_connectedness_t mode);
IGRAPH_EXPORT igraph_error_t igraph_is_connected(const igraph_t *graph, igraph_bool_t *res,
                                      igraph_connectedness_t mode);
IGRAPH_EXPORT igraph_error_t igraph_decompose(const igraph_t *graph, igraph_graph_list_t *components,
                                   igraph_connectedness_t mode,
                                   igraph_int_t maxcompno, igraph_int_t minelements);
IGRAPH_EXPORT igraph_error_t igraph_articulation_points(const igraph_t *graph,
                                             igraph_vector_int_t *res);
IGRAPH_EXPORT igraph_error_t igraph_biconnected_components(const igraph_t *graph,
                                                igraph_int_t *no,
                                                igraph_vector_int_list_t *tree_edges,
                                                igraph_vector_int_list_t *component_edges,
                                                igraph_vector_int_list_t *components,
                                                igraph_vector_int_t *articulation_points);
IGRAPH_EXPORT igraph_error_t igraph_is_biconnected(const igraph_t *graph, igraph_bool_t *result);
IGRAPH_EXPORT igraph_error_t igraph_bridges(const igraph_t *graph, igraph_vector_int_t *bridges);

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_bond_percolation(
        const igraph_t *graph,
        igraph_vector_int_t *giant_size,
        igraph_vector_int_t *vertex_count,
        const igraph_vector_int_t *edge_order);
IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_site_percolation(
        const igraph_t *graph,
        igraph_vector_int_t *giant_size,
        igraph_vector_int_t *edge_count,
        const igraph_vector_int_t *vertex_order);
IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_edgelist_percolation(
        const igraph_vector_int_t *edges,
        igraph_vector_int_t *giant_size,
        igraph_vector_int_t *vertex_count);

IGRAPH_END_C_DECLS

#endif
