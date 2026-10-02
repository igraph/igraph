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

#ifndef IGRAPH_PATH_COVER_H
#define IGRAPH_PATH_COVER_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_vector.h"
#include "igraph_vector_list.h"

IGRAPH_BEGIN_C_DECLS

/* This module is a port of the Minimum Path Cover (MPC) algorithms from
 * https://github.com/algbio/PerformanceMPC (src/mpc/naive.cpp, antichain.cpp,
 * cc.cpp), reimplemented on top of igraph's own flow and graph primitives.
 * The underlying minimum-flow-with-lower-bounds reduction and the greedy
 * initial solution follow Mäkinen, Tomescu, Kuosmanen, Paavilainen, Gagie &
 * Chikhi, "Sparse Dynamic Programming on DAGs with Small Width", ACM
 * Transactions on Algorithms 15(2):29, 2019 (Section 2); see the .c files
 * for exact per-function references. */

/**
 * \typedef igraph_mpc_reduction_t
 * \brief Initial feasible solution strategy for minimum path cover.
 *
 * \ref igraph_minimum_path_cover() works by building an initial feasible
 * solution to an equivalent minimum-flow-with-lower-bounds problem, and
 * then shrinking it down to the true minimum. This type selects the
 * strategy used to build the initial solution.
 *
 * \enumval IGRAPH_MPC_REDUCTION_NAIVE Every vertex starts out as its own
 *    trivial one-vertex path. Cheapest to build, in O(|V|) time, but
 *    leaves the most work for the solver, as it starts from the largest
 *    possible width (equal to the number of vertices).
 * \enumval IGRAPH_MPC_REDUCTION_GREEDY A greedy longest-uncovered-chain
 *    heuristic is used to build a much smaller initial feasible solution,
 *    which usually reduces the amount of work left for the solver.
 *    Only supported together with uniform (or absent) \c vertex_weights.
 */
typedef enum {
    IGRAPH_MPC_REDUCTION_NAIVE = 0,
    IGRAPH_MPC_REDUCTION_GREEDY = 1
} igraph_mpc_reduction_t;

/**
 * \typedef igraph_mpc_solver_t
 * \brief Solver backend for minimum path cover.
 *
 * Selects the algorithm used to shrink a feasible solution of the
 * minimum-flow-with-lower-bounds problem underlying
 * \ref igraph_minimum_path_cover() down to the true minimum.
 *
 * \enumval IGRAPH_MPC_SOLVER_NAIVE_DFS A direct Ford-Fulkerson-style DFS
 *    search over the demand-residual graph, cancelling slack flow one
 *    augmenting path at a time.
 * \enumval IGRAPH_MPC_SOLVER_MAXFLOW_REDUCTION Reduces the problem to a
 *    single ordinary maximum flow instance, solved with
 *    \ref igraph_maxflow().
 */
typedef enum {
    IGRAPH_MPC_SOLVER_NAIVE_DFS = 0,
    IGRAPH_MPC_SOLVER_MAXFLOW_REDUCTION = 1
} igraph_mpc_solver_t;

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_minimum_path_cover(
        const igraph_t *graph,
        igraph_vector_int_list_t *cover,
        igraph_int_t *width,
        const igraph_vector_int_t *vertex_weights,
        igraph_mpc_reduction_t reduction,
        igraph_mpc_solver_t solver,
        igraph_t *flow_network,
        igraph_vector_int_t *flow_network_flow,
        igraph_vector_int_t *flow_network_demand);

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_maximum_antichain(
        const igraph_t *flow_network,
        const igraph_vector_int_t *flow,
        const igraph_vector_int_t *demand,
        igraph_vector_int_t *antichain);

IGRAPH_EXPERIMENTAL IGRAPH_EXPORT igraph_error_t igraph_minimum_chain_cover(
        const igraph_t *graph,
        const igraph_vector_int_list_t *cover,
        igraph_vector_int_list_t *chain_cover);

IGRAPH_END_C_DECLS

#endif
