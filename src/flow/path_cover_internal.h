/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_PATH_COVER_INTERNAL_H
#define IGRAPH_PATH_COVER_INTERNAL_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_list.h"

IGRAPH_BEGIN_C_DECLS

/* Ported from https://github.com/algbio/PerformanceMPC (src/mpc/naive.cpp,
 * antichain.cpp, cc.cpp). The vertex-splitting reduction below is the
 * "standard reduction from the minimum path cover problem to a minimum flow
 * one" of Mäkinen, Tomescu, Kuosmanen, Paavilainen, Gagie & Chikhi, "Sparse
 * Dynamic Programming on DAGs with Small Width", ACM Transactions on
 * Algorithms 15(2):29, 2019, Section 2; see path_cover.c for the reduction
 * and solver implementations, each annotated with its exact source.
 *
 * Vertex-splitting layout of the minimum path cover flow network, built on
 * top of igraph_i_split_vertices() (src/flow/flow_conversion.c). For an
 * original vertex v in [0, n): the "output half" is v_out(v) = v and the
 * "input half" is v_in(v) = n + v, exactly as igraph_i_split_vertices()
 * lays them out. Two extra vertices are appended: source = 2n, sink = 2n+1.
 *
 * With m = original edge count, the network has m + 3n edges, laid out as:
 *   [0, m)          original DAG edges i->u,   as (i, n+u)
 *   [m, m+n)         vertex demand edges,       as (n+i, i),   id m+i
 *   [m+n, m+2n)      source edges,               as (2n, n+i), id m+n+i
 *   [m+2n, m+3n)     sink edges,                 as (i, 2n+1), id m+2n+i
 *
 * This layout is identical for every reduction strategy; only the flow
 * values differ. The demand vector only depends on vertex_weights.
 */

static inline igraph_int_t igraph_i_mpc_v_out(igraph_int_t v) {
    return v;
}

static inline igraph_int_t igraph_i_mpc_v_in(igraph_int_t v, igraph_int_t no_of_nodes) {
    return no_of_nodes + v;
}

static inline igraph_int_t igraph_i_mpc_v_r(igraph_int_t x, igraph_int_t no_of_nodes) {
    return x < no_of_nodes ? x : x - no_of_nodes;
}

igraph_error_t igraph_i_mpc_build_network(
        const igraph_t *graph, igraph_t *network,
        igraph_vector_int_t *demand,
        const igraph_vector_int_t *vertex_weights);

igraph_error_t igraph_i_mpc_is_valid_minflow(
        const igraph_t *network, const igraph_vector_int_t *flow,
        const igraph_vector_int_t *demand, igraph_bool_t *valid);

igraph_error_t igraph_i_mpc_naive_reduction(
        const igraph_t *graph, const igraph_vector_int_t *vertex_weights,
        igraph_vector_int_t *flow);

igraph_error_t igraph_i_mpc_greedy_reduction(
        const igraph_t *graph, igraph_vector_int_t *flow);

igraph_error_t igraph_i_mpc_solve_naive_dfs(
        const igraph_t *network, const igraph_vector_int_t *demand,
        igraph_int_t source, igraph_int_t sink, igraph_vector_int_t *flow);

igraph_error_t igraph_i_mpc_solve_maxflow_reduction(
        const igraph_t *network, const igraph_vector_int_t *demand,
        igraph_int_t source, igraph_int_t sink, igraph_vector_int_t *flow);

igraph_error_t igraph_i_mpc_recover_paths(
        const igraph_t *network, igraph_vector_int_t *flow,
        igraph_int_t source, igraph_int_t sink, igraph_int_t no_of_nodes,
        igraph_vector_int_list_t *cover);

IGRAPH_END_C_DECLS

#endif
