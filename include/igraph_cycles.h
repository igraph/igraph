/*
   igraph library.
   Copyright (C) 2022-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Directed acyclic graphs                            */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_topological_sorting(
        const igraph_t *graph, igraph_vector_int_t *res, igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_is_dag(const igraph_t *graph, igraph_bool_t *res);

/* -------------------------------------------------- */
/* Cycle bases                                        */
/* -------------------------------------------------- */

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_fundamental_cycles(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_int_list_t *result,
        igraph_int_t start_vid, igraph_real_t bfs_cutoff);

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_minimum_cycle_basis(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_int_list_t *result,
        igraph_real_t bfs_cutoff,
        igraph_bool_t complete, igraph_bool_t use_cycle_order);

IGRAPH_EXPORT igraph_error_t igraph_find_cycle(
        const igraph_t *graph,
        igraph_vector_int_t *vertices,
        igraph_vector_int_t *edges,
        igraph_neimode_t mode);

/**
 * \typedef igraph_cycle_handler_t
 * \brief Type of cycle handler functions.
 *
 * Callback type, called by \ref igraph_simple_cycles_callback() when
 * a cycle is found.
 *
 * \param vertices The vertices of the current cycle. Must not be modified.
 * \param edges The edges of the current cycle. Must not be modified.
 * \param arg The extra parameter passed to \ref igraph_simple_cycles_callback()
 * \return Error code; \c IGRAPH_SUCCESS to continue the search or
 *   \c IGRAPH_STOP to stop the search without signaling an error.
 */
typedef igraph_error_t igraph_cycle_handler_t(
        const igraph_vector_int_t *vertices,
        const igraph_vector_int_t *edges,
        void *arg);

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_simple_cycles_callback(
        const igraph_t *graph,
        igraph_neimode_t mode,
        igraph_int_t min_cycle_length,
        igraph_int_t max_cycle_length,
        igraph_cycle_handler_t *callback, void *arg);

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_simple_cycles(
        const igraph_t *graph,
        igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
        igraph_neimode_t mode,
        igraph_int_t min_cycle_length, igraph_int_t max_cycle_length,
        igraph_int_t max_results);

IGRAPH_EXPORT igraph_error_t igraph_feedback_arc_set(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights, igraph_fas_algorithm_t algo);

IGRAPH_EXPORT igraph_error_t igraph_feedback_vertex_set(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *vertex_weights, igraph_fvs_algorithm_t algo);

IGRAPH_END_C_DECLS

#endif
