/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_PATHS_INTERNAL_H
#define IGRAPH_PATHS_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_paths.h"

IGRAPH_BEGIN_C_DECLS

/* Helper functions for validating input to shortest path functions. */

igraph_error_t igraph_i_validate_distance_weights(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_bool_t *negative_weights);


/* The following shortest path functions skip most input validation,
 * and are meant for internal use in context where the validation
 * has already been done. */

/* Distances */

igraph_error_t igraph_i_distances_unweighted_cutoff(
        const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_vs_t from, igraph_vs_t to,
        igraph_neimode_t mode,
        igraph_real_t cutoff);

igraph_error_t igraph_i_distances_dijkstra_cutoff(
        const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_vs_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode,
        igraph_real_t cutoff);

igraph_error_t igraph_i_distances_bellman_ford(
        const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_vs_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode);

igraph_error_t igraph_i_distances_johnson(
        const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_vs_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode);

igraph_error_t igraph_i_distances_floyd_warshall(
        const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_vs_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode,
        igraph_floyd_warshall_algorithm_t method);

/* Get shortest paths */

igraph_error_t igraph_i_get_shortest_paths_unweighted(
        const igraph_t *graph,
        igraph_vector_int_list_t *vertices,
        igraph_vector_int_list_t *edges,
        igraph_int_t from, igraph_vs_t to,
        igraph_neimode_t mode,
        igraph_vector_int_t *parents,
        igraph_vector_int_t *inbound_edges);

igraph_error_t igraph_i_get_shortest_paths_dijkstra(
        const igraph_t *graph,
        igraph_vector_int_list_t *vertices,
        igraph_vector_int_list_t *edges,
        igraph_int_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode,
        igraph_vector_int_t *parents,
        igraph_vector_int_t *inbound_edges);

igraph_error_t igraph_i_get_shortest_paths_bellman_ford(
        const igraph_t *graph,
        igraph_vector_int_list_t *vertices,
        igraph_vector_int_list_t *edges,
        igraph_int_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode,
        igraph_vector_int_t *parents,
        igraph_vector_int_t *inbound_edges);

igraph_error_t igraph_i_get_all_shortest_paths_unweighted(
        const igraph_t *graph,
        igraph_vector_int_list_t *vertices,
        igraph_vector_int_list_t *edges,
        igraph_vector_int_t *nrgeo,
        igraph_int_t from, igraph_vs_t to,
        igraph_neimode_t mode);

IGRAPH_END_C_DECLS

#endif /* IGRAPH_PATHS_INTERNAL_H */
