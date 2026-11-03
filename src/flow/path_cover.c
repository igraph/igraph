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

#include "igraph_path_cover.h"

#include "igraph_adjlist.h"
#include "igraph_constants.h"
#include "igraph_constructors.h"
#include "igraph_cycles.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_flow.h"
#include "igraph_interface.h"

#include "flow/flow_internal.h"
#include "flow/path_cover_internal.h"

/* -------------------------------------------------------------------- */
/* Network construction                                                 */
/* -------------------------------------------------------------------- */

/* Standard MPC-to-minimum-flow reduction: split each vertex v into
 * (v_in, v_out) with demand 1 on the (v_in, v_out) edge, keep original
 * edges as (u_out, v_in), and attach a global source/sink. See Mäkinen,
 * Tomescu, Kuosmanen, Paavilainen, Gagie & Chikhi, "Sparse Dynamic
 * Programming on DAGs with Small Width", ACM TALG 15(2):29, 2019, Section
 * 2, paragraph "The standard reduction from the minimum path cover problem
 * to a minimum flow one...". Ported from naive_minflow_reduction() /
 * greedy_minflow_reduction() in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/naive.cpp,
 * which construct this same network inline for each reduction; here it is
 * factored out since the topology is reduction-independent. */
igraph_error_t igraph_i_mpc_build_network(
        const igraph_t *graph, igraph_t *network,
        igraph_vector_int_t *demand,
        const igraph_vector_int_t *vertex_weights) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_int_t total_edges = no_of_edges + 3 * no_of_nodes;
    igraph_int_t source = 2 * no_of_nodes;
    igraph_int_t sink = 2 * no_of_nodes + 1;
    igraph_vector_int_t new_edges;
    igraph_int_t i;

    IGRAPH_CHECK(igraph_i_split_vertices(graph, network));
    IGRAPH_FINALLY(igraph_destroy, network);

    IGRAPH_CHECK(igraph_add_vertices(network, 2, NULL));

    IGRAPH_VECTOR_INT_INIT_FINALLY(&new_edges, 4 * no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        /* source -> v_in(i) */
        VECTOR(new_edges)[2 * i] = source;
        VECTOR(new_edges)[2 * i + 1] = igraph_i_mpc_v_in(i, no_of_nodes);
        /* v_out(i) -> sink */
        VECTOR(new_edges)[2 * no_of_nodes + 2 * i] = igraph_i_mpc_v_out(i);
        VECTOR(new_edges)[2 * no_of_nodes + 2 * i + 1] = sink;
    }
    IGRAPH_CHECK(igraph_add_edges(network, &new_edges, NULL));
    igraph_vector_int_destroy(&new_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INT_INIT_FINALLY(demand, total_edges);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*demand)[no_of_edges + i] = vertex_weights ? VECTOR(*vertex_weights)[i] : 1;
    }
    IGRAPH_FINALLY_CLEAN(1); /* demand */

    IGRAPH_FINALLY_CLEAN(1); /* network */

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_mpc_is_valid_minflow(
        const igraph_t *network, const igraph_vector_int_t *flow,
        const igraph_vector_int_t *demand, igraph_bool_t *valid) {

    igraph_int_t no_of_nodes = igraph_vcount(network);
    igraph_int_t no_of_edges = igraph_ecount(network);
    igraph_vector_int_t balance;
    igraph_int_t i;

    *valid = true;

    for (i = 0; i < no_of_edges; i++) {
        if (VECTOR(*flow)[i] < VECTOR(*demand)[i]) {
            *valid = false;
            return IGRAPH_SUCCESS;
        }
    }

    /* The network always has source = no_of_nodes - 2 and sink = no_of_nodes - 1,
     * i.e. they are the last two vertices; flow conservation is only required
     * everywhere else. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&balance, no_of_nodes);
    for (i = 0; i < no_of_edges; i++) {
        VECTOR(balance)[IGRAPH_FROM(network, i)] -= VECTOR(*flow)[i];
        VECTOR(balance)[IGRAPH_TO(network, i)] += VECTOR(*flow)[i];
    }
    for (i = 0; i < no_of_nodes - 2; i++) {
        if (VECTOR(balance)[i] != 0) {
            *valid = false;
            break;
        }
    }
    igraph_vector_int_destroy(&balance);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Reductions: build an initial feasible flow                          */
/* -------------------------------------------------------------------- */

/* Trivial feasible flow: every vertex is its own one-vertex path, i.e. the
 * initial path cover has width |V|. Ported from naive_minflow_reduction()
 * in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/naive.cpp. */
igraph_error_t igraph_i_mpc_naive_reduction(
        const igraph_t *graph, const igraph_vector_int_t *vertex_weights,
        igraph_vector_int_t *flow) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_int_t i;

    igraph_vector_int_null(flow);
    for (i = 0; i < no_of_nodes; i++) {
        igraph_int_t w = vertex_weights ? VECTOR(*vertex_weights)[i] : 1;
        VECTOR(*flow)[no_of_edges + i] = w;                   /* v_in(i) -> v_out(i) */
        VECTOR(*flow)[no_of_edges + no_of_nodes + i] = w;      /* source -> v_in(i) */
        VECTOR(*flow)[no_of_edges + 2 * no_of_nodes + i] = w;  /* v_out(i) -> sink */
    }

    return IGRAPH_SUCCESS;
}

/* Greedy longest-uncovered-chain heuristic: repeatedly pick the path
 * covering the most still-uncovered vertices via a topological-order DP
 * (max_len[v] = 1{v uncovered} + max over out-neighbors u of max_len[u]),
 * mark it covered, and repeat. This is the O(k log|V|)-width greedy
 * set-cover algorithm and DP of Mäkinen, Tomescu, Kuosmanen, Paavilainen,
 * Gagie & Chikhi, "Sparse Dynamic Programming on DAGs with Small Width",
 * ACM TALG 15(2):29, 2019, Section 2, Lemma 2.1 and its proof. Ported from
 * greedy_minflow_reduction() in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/naive.cpp.
 * Only supports uniform (all-1) vertex weights; the caller
 * (igraph_minimum_path_cover) is responsible for rejecting non-uniform
 * weights before calling this function. */
igraph_error_t igraph_i_mpc_greedy_reduction(
        const igraph_t *graph, igraph_vector_int_t *flow) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t topo;
    igraph_inclist_t out_inc;
    igraph_vector_int_t max_len, from_vertex, from_edge;
    igraph_vector_bool_t not_covered;
    igraph_vector_int_t chain_edge_flow;
    igraph_vector_int_t node_flow, node_source, node_sink;
    igraph_int_t i, remaining = no_of_nodes;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&topo, 0);
    IGRAPH_CHECK(igraph_topological_sorting(graph, &topo, IGRAPH_OUT));

    IGRAPH_CHECK(igraph_inclist_init(graph, &out_inc, IGRAPH_OUT, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &out_inc);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&max_len, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&from_vertex, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&from_edge, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&not_covered, no_of_nodes);
    igraph_vector_bool_fill(&not_covered, true);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&chain_edge_flow, no_of_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&node_flow, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&node_source, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&node_sink, no_of_nodes);

    while (remaining > 0) {
        igraph_int_t best_vertex = -1, best_len = 0;

        igraph_vector_int_null(&max_len);
        igraph_vector_int_fill(&from_vertex, -1);
        igraph_vector_int_fill(&from_edge, -1);

        for (i = 0; i < no_of_nodes; i++) {
            igraph_int_t s = VECTOR(topo)[i];
            igraph_vector_int_t *neis;
            igraph_int_t j, nlen;

            /* Already-covered vertices break chains: they must not relay an
             * inherited max_len to their successors, or successors would
             * wrongly extend a chain "through" a vertex that contributes
             * nothing, causing that vertex to be walked (and its coverage
             * counter double-decremented) when the chain is later applied. */
            if (!VECTOR(not_covered)[s]) {
                continue;
            }

            VECTOR(max_len)[s]++;
            if (VECTOR(max_len)[s] > best_len) {
                best_len = VECTOR(max_len)[s];
                best_vertex = s;
            }

            neis = igraph_inclist_get(&out_inc, s);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_int_t eid = VECTOR(*neis)[j];
                igraph_int_t u = IGRAPH_TO(graph, eid);
                if (VECTOR(max_len)[s] > VECTOR(max_len)[u]) {
                    VECTOR(max_len)[u] = VECTOR(max_len)[s];
                    VECTOR(from_vertex)[u] = s;
                    VECTOR(from_edge)[u] = eid;
                }
            }
        }

        if (best_vertex == -1) {
            break;
        }

        {
            igraph_int_t cur = best_vertex;
            VECTOR(node_sink)[cur]++;
            while (VECTOR(from_vertex)[cur] != -1) {
                igraph_int_t eid = VECTOR(from_edge)[cur];
                igraph_int_t prev = VECTOR(from_vertex)[cur];
                VECTOR(not_covered)[cur] = false;
                remaining--;
                VECTOR(chain_edge_flow)[eid]++;
                VECTOR(node_flow)[cur]++;
                cur = prev;
            }
            VECTOR(not_covered)[cur] = false;
            remaining--;
            VECTOR(node_flow)[cur]++;
            VECTOR(node_source)[cur]++;
        }
    }

    igraph_vector_int_null(flow);
    for (i = 0; i < no_of_edges; i++) {
        VECTOR(*flow)[i] = VECTOR(chain_edge_flow)[i];
    }
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*flow)[no_of_edges + i] = VECTOR(node_flow)[i];
        VECTOR(*flow)[no_of_edges + no_of_nodes + i] = VECTOR(node_source)[i];
        VECTOR(*flow)[no_of_edges + 2 * no_of_nodes + i] = VECTOR(node_sink)[i];
    }

    igraph_vector_int_destroy(&node_sink);
    igraph_vector_int_destroy(&node_source);
    igraph_vector_int_destroy(&node_flow);
    igraph_vector_int_destroy(&chain_edge_flow);
    igraph_vector_bool_destroy(&not_covered);
    igraph_vector_int_destroy(&from_edge);
    igraph_vector_int_destroy(&from_vertex);
    igraph_vector_int_destroy(&max_len);
    igraph_inclist_destroy(&out_inc);
    igraph_vector_int_destroy(&topo);
    IGRAPH_FINALLY_CLEAN(10);

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Solvers: shrink a feasible flow down to the minimum                  */
/* -------------------------------------------------------------------- */

/* Direct Ford-Fulkerson-style DFS over the demand-residual graph. Forward
 * (out) edges are traversable while they carry slack above their demand
 * (flow > demand), and taking one decreases the flow; backward (in) edges
 * are always traversable, and taking one increases the flow. Both the
 * visited marks and the per-vertex edge cursors are reset at the start of
 * every augmenting-path search, since flow can both increase and decrease
 * across searches here (unlike igraph_i_mpc_recover_paths, where flow is
 * monotonically drained).
 *
 * Correctness (termination with a minimum flow once no more augmenting path
 * exists) rests on the max-flow/min-cut argument of Ford & Fulkerson,
 * "Maximal Flow Through a Network", Canadian Journal of Mathematics 8,
 * 1956, Theorem 1 (the Minimal Cut Theorem), applied to the demand-residual
 * graph -- matching PerformanceMPC's own citation for this routine (see its
 * README, reference [2]). Ported from naive_minflow_solve() in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/naive.cpp. */
igraph_error_t igraph_i_mpc_solve_naive_dfs(
        const igraph_t *network, const igraph_vector_int_t *demand,
        igraph_int_t source, igraph_int_t sink, igraph_vector_int_t *flow) {

    igraph_int_t no_of_nodes = igraph_vcount(network);
    igraph_inclist_t out_inc, in_inc;
    igraph_vector_int_t out_cursor, in_cursor;
    igraph_vector_bool_t visited;
    igraph_vector_int_t stack_vertex, stack_edge;
    igraph_vector_bool_t stack_reverse;

    IGRAPH_CHECK(igraph_inclist_init(network, &out_inc, IGRAPH_OUT, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &out_inc);
    IGRAPH_CHECK(igraph_inclist_init(network, &in_inc, IGRAPH_IN, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &in_inc);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&out_cursor, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&in_cursor, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&visited, no_of_nodes);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&stack_vertex, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&stack_edge, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&stack_reverse, 0);

    for ( ; ; ) {
        igraph_bool_t found = false;

        igraph_vector_int_null(&out_cursor);
        igraph_vector_int_null(&in_cursor);
        igraph_vector_bool_null(&visited);
        igraph_vector_int_clear(&stack_vertex);
        igraph_vector_int_clear(&stack_edge);
        igraph_vector_bool_clear(&stack_reverse);

        IGRAPH_CHECK(igraph_vector_int_push_back(&stack_vertex, source));
        IGRAPH_CHECK(igraph_vector_int_push_back(&stack_edge, -1));
        IGRAPH_CHECK(igraph_vector_bool_push_back(&stack_reverse, false));
        VECTOR(visited)[source] = true;

        while (!igraph_vector_int_empty(&stack_vertex)) {
            igraph_int_t s = igraph_vector_int_tail(&stack_vertex);
            igraph_bool_t advanced = false;
            igraph_vector_int_t *neis;
            igraph_int_t nlen;

            if (s == sink) {
                found = true;
                break;
            }

            neis = igraph_inclist_get(&out_inc, s);
            nlen = igraph_vector_int_size(neis);
            while (VECTOR(out_cursor)[s] < nlen) {
                igraph_int_t eid = VECTOR(*neis)[VECTOR(out_cursor)[s]];
                igraph_int_t u = IGRAPH_TO(network, eid);
                VECTOR(out_cursor)[s]++;
                if (!VECTOR(visited)[u] && VECTOR(*flow)[eid] > VECTOR(*demand)[eid]) {
                    VECTOR(visited)[u] = true;
                    IGRAPH_CHECK(igraph_vector_int_push_back(&stack_vertex, u));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&stack_edge, eid));
                    IGRAPH_CHECK(igraph_vector_bool_push_back(&stack_reverse, false));
                    advanced = true;
                    break;
                }
            }

            if (!advanced) {
                neis = igraph_inclist_get(&in_inc, s);
                nlen = igraph_vector_int_size(neis);
                while (VECTOR(in_cursor)[s] < nlen) {
                    igraph_int_t eid = VECTOR(*neis)[VECTOR(in_cursor)[s]];
                    igraph_int_t u = IGRAPH_FROM(network, eid);
                    VECTOR(in_cursor)[s]++;
                    if (!VECTOR(visited)[u]) {
                        VECTOR(visited)[u] = true;
                        IGRAPH_CHECK(igraph_vector_int_push_back(&stack_vertex, u));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&stack_edge, eid));
                        IGRAPH_CHECK(igraph_vector_bool_push_back(&stack_reverse, true));
                        advanced = true;
                        break;
                    }
                }
            }

            if (!advanced) {
                igraph_vector_int_pop_back(&stack_vertex);
                igraph_vector_int_pop_back(&stack_edge);
                igraph_vector_bool_pop_back(&stack_reverse);
            }
        }

        if (!found) {
            break;
        }

        {
            igraph_int_t len = igraph_vector_int_size(&stack_edge);
            igraph_int_t i;
            for (i = 1; i < len; i++) {
                igraph_int_t eid = VECTOR(stack_edge)[i];
                if (VECTOR(stack_reverse)[i]) {
                    VECTOR(*flow)[eid]++;
                } else {
                    VECTOR(*flow)[eid]--;
                }
            }
        }
    }

    igraph_vector_bool_destroy(&stack_reverse);
    igraph_vector_int_destroy(&stack_edge);
    igraph_vector_int_destroy(&stack_vertex);
    igraph_vector_bool_destroy(&visited);
    igraph_vector_int_destroy(&in_cursor);
    igraph_vector_int_destroy(&out_cursor);
    igraph_inclist_destroy(&in_inc);
    igraph_inclist_destroy(&out_inc);
    IGRAPH_FINALLY_CLEAN(8);

    return IGRAPH_SUCCESS;
}

/* Reduces the lower-bound flow problem to an ordinary maximum flow problem
 * on an auxiliary graph, solved with igraph_maxflow(). For every network
 * edge with slack above its demand, a forward edge is added with capacity
 * equal to that slack; for every network edge, a reverse edge is added
 * with capacity equal to the total flow leaving the source (a safe upper
 * bound for this instance). The resulting max-flow is translated back onto
 * the original flow: forward-tagged edges cancel flow, reverse-tagged
 * edges add it back.
 *
 * This is the minimum-flow-to-maximum-flow reduction of Mäkinen, Tomescu,
 * Kuosmanen, Paavilainen, Gagie & Chikhi, "Sparse Dynamic Programming on
 * DAGs with Small Width", ACM TALG 15(2):29, 2019, Section 2: "(i) find a
 * feasible flow f; (ii) transform this into a minimum feasible flow, by
 * finding a maximum flow f' in G in which every e now has capacity
 * f(e)-d(e). The final minimum flow solution is obtained as f(e)-f'(e)."
 * Ported from minflow_maxflow_reduction() in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/naive.cpp,
 * but calling igraph_maxflow() directly instead of a hand-rolled
 * Edmonds-Karp/Dinic backend. */
igraph_error_t igraph_i_mpc_solve_maxflow_reduction(
        const igraph_t *network, const igraph_vector_int_t *demand,
        igraph_int_t source, igraph_int_t sink, igraph_vector_int_t *flow) {

    igraph_int_t no_of_nodes = igraph_vcount(network);
    igraph_int_t no_of_edges = igraph_ecount(network);
    igraph_int_t total_source_flow = 0;
    igraph_t reduction_graph;
    igraph_vector_int_t red_edges;
    igraph_vector_t red_capacity;
    igraph_vector_int_t red_orig_edge;
    igraph_vector_bool_t red_is_reverse;
    igraph_vector_t maxflow_result;
    igraph_vector_int_t out_eids;
    igraph_int_t i;

    {
        igraph_bool_t valid;
        IGRAPH_CHECK(igraph_i_mpc_is_valid_minflow(network, flow, demand, &valid));
        IGRAPH_ASSERT(valid);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&out_eids, 0);
    IGRAPH_CHECK(igraph_incident(network, &out_eids, source, IGRAPH_OUT, IGRAPH_LOOPS));
    for (i = 0; i < igraph_vector_int_size(&out_eids); i++) {
        total_source_flow += VECTOR(*flow)[VECTOR(out_eids)[i]];
    }
    igraph_vector_int_destroy(&out_eids);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&red_edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&red_capacity, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&red_orig_edge, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&red_is_reverse, 0);

    for (i = 0; i < no_of_edges; i++) {
        igraph_int_t from = IGRAPH_FROM(network, i);
        igraph_int_t to = IGRAPH_TO(network, i);
        igraph_int_t slack = VECTOR(*flow)[i] - VECTOR(*demand)[i];

        if (slack > 0) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&red_edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(&red_edges, to));
            IGRAPH_CHECK(igraph_vector_push_back(&red_capacity, (igraph_real_t) slack));
            IGRAPH_CHECK(igraph_vector_int_push_back(&red_orig_edge, i));
            IGRAPH_CHECK(igraph_vector_bool_push_back(&red_is_reverse, false));
        }

        IGRAPH_CHECK(igraph_vector_int_push_back(&red_edges, to));
        IGRAPH_CHECK(igraph_vector_int_push_back(&red_edges, from));
        IGRAPH_CHECK(igraph_vector_push_back(&red_capacity, (igraph_real_t) total_source_flow));
        IGRAPH_CHECK(igraph_vector_int_push_back(&red_orig_edge, i));
        IGRAPH_CHECK(igraph_vector_bool_push_back(&red_is_reverse, true));
    }

    IGRAPH_CHECK(igraph_create(&reduction_graph, &red_edges, no_of_nodes, IGRAPH_DIRECTED));
    IGRAPH_FINALLY(igraph_destroy, &reduction_graph);

    IGRAPH_VECTOR_INIT_FINALLY(&maxflow_result, 0);
    IGRAPH_CHECK(igraph_maxflow(&reduction_graph, NULL, &maxflow_result, NULL, NULL, NULL,
                                 source, sink, &red_capacity, NULL));

    {
        igraph_int_t red_len = igraph_vector_int_size(&red_orig_edge);
        for (i = 0; i < red_len; i++) {
            igraph_int_t orig_eid = VECTOR(red_orig_edge)[i];
            igraph_int_t amount = (igraph_int_t) (VECTOR(maxflow_result)[i] + 0.5);
            if (VECTOR(red_is_reverse)[i]) {
                VECTOR(*flow)[orig_eid] += amount;
            } else {
                VECTOR(*flow)[orig_eid] -= amount;
            }
        }
    }

    igraph_vector_destroy(&maxflow_result);
    igraph_destroy(&reduction_graph);
    igraph_vector_bool_destroy(&red_is_reverse);
    igraph_vector_int_destroy(&red_orig_edge);
    igraph_vector_destroy(&red_capacity);
    igraph_vector_int_destroy(&red_edges);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Path recovery: decompose a solved minimum flow into paths             */
/* -------------------------------------------------------------------- */

/* Forward-only DFS decomposing a solved (minimum) flow into paths (which
 * may share vertices at merge points -- see igraph_minimum_path_cover()).
 * Unlike igraph_i_mpc_solve_naive_dfs(), the per-vertex edge cursors
 * are NOT reset between passes -- only the visited marks are. This is safe
 * because the network is a DAG and flow here only ever decreases (never
 * increases) across passes: once a vertex is exhaustively found (within a
 * pass) unable to reach the sink via any then-positive-flow edge, it can
 * never become able to later, since future passes only see equal-or-lower
 * flow values. So permanently skipping edges that are spent (flow == 0) or
 * that lead to an already-visited-this-pass vertex is correct, and gives
 * amortized O(total path length + |E|) time instead of the naive solver's
 * O(paths * (|V|+|E|)).
 *
 * The decomposition of a feasible flow into source-to-sink paths is
 * folklore flow theory (a direct consequence of Ford & Fulkerson's flow
 * conservation laws, "Maximal Flow Through a Network", Canadian Journal of
 * Mathematics 8, 1956); the persistent-cursor optimization used here is
 * ported from minflow_reduction_path_recover_faster() in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/naive.cpp
 * (chosen over that file's plain and "_fast" variants, which it strictly
 * dominates). */
igraph_error_t igraph_i_mpc_recover_paths(
        const igraph_t *network, igraph_vector_int_t *flow,
        igraph_int_t source, igraph_int_t sink, igraph_int_t no_of_nodes,
        igraph_vector_int_list_t *cover) {

    igraph_int_t network_nodes = igraph_vcount(network);
    igraph_inclist_t out_inc;
    igraph_vector_int_t cursor;
    igraph_vector_bool_t visited;
    igraph_vector_int_t stack_vertex, stack_edge;
    igraph_int_t i;

    IGRAPH_CHECK(igraph_inclist_init(network, &out_inc, IGRAPH_OUT, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &out_inc);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&cursor, network_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&visited, network_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&stack_vertex, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&stack_edge, 0);

    for ( ; ; ) {
        igraph_bool_t found = false;

        igraph_vector_bool_null(&visited);
        igraph_vector_int_clear(&stack_vertex);
        igraph_vector_int_clear(&stack_edge);

        IGRAPH_CHECK(igraph_vector_int_push_back(&stack_vertex, source));
        IGRAPH_CHECK(igraph_vector_int_push_back(&stack_edge, -1));
        VECTOR(visited)[source] = true;

        while (!igraph_vector_int_empty(&stack_vertex)) {
            igraph_int_t s = igraph_vector_int_tail(&stack_vertex);
            igraph_bool_t advanced = false;
            igraph_vector_int_t *neis;
            igraph_int_t nlen;

            if (s == sink) {
                found = true;
                break;
            }

            neis = igraph_inclist_get(&out_inc, s);
            nlen = igraph_vector_int_size(neis);
            while (VECTOR(cursor)[s] < nlen) {
                igraph_int_t eid = VECTOR(*neis)[VECTOR(cursor)[s]];

                if (VECTOR(*flow)[eid] <= 0) {
                    VECTOR(cursor)[s]++;
                    continue;
                }
                {
                    igraph_int_t u = IGRAPH_TO(network, eid);
                    if (VECTOR(visited)[u]) {
                        VECTOR(cursor)[s]++;
                        continue;
                    }
                    VECTOR(visited)[u] = true;
                    IGRAPH_CHECK(igraph_vector_int_push_back(&stack_vertex, u));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&stack_edge, eid));
                    advanced = true;
                    break;
                }
            }

            if (!advanced) {
                igraph_vector_int_pop_back(&stack_vertex);
                igraph_vector_int_pop_back(&stack_edge);
            }
        }

        if (!found) {
            break;
        }

        {
            igraph_int_t len = igraph_vector_int_size(&stack_vertex);
            igraph_vector_int_t path;

            IGRAPH_VECTOR_INT_INIT_FINALLY(&path, 0);
            /* stack_vertex alternates source, v_in(x1), v_out(x1)=x1, v_in(x2),
             * v_out(x2)=x2, ..., sink; the original vertices are the v_out
             * halves, at even indices starting from 2. */
            for (i = 2; i < len - 1; i += 2) {
                igraph_int_t node = VECTOR(stack_vertex)[i];
                IGRAPH_CHECK(igraph_vector_int_push_back(&path, igraph_i_mpc_v_r(node, no_of_nodes)));
            }
            for (i = 1; i < len; i++) {
                igraph_int_t eid = VECTOR(stack_edge)[i];
                VECTOR(*flow)[eid]--;
            }
            IGRAPH_CHECK(igraph_vector_int_list_push_back(cover, &path));
            IGRAPH_FINALLY_CLEAN(1); /* ownership of path transferred to cover */
        }
    }

    igraph_vector_int_destroy(&stack_edge);
    igraph_vector_int_destroy(&stack_vertex);
    igraph_vector_bool_destroy(&visited);
    igraph_vector_int_destroy(&cursor);
    igraph_inclist_destroy(&out_inc);
    IGRAPH_FINALLY_CLEAN(5);

    /* Strong correctness postcondition: every unit of flow has been
     * decomposed into some path. */
    {
        igraph_int_t total_edges = igraph_ecount(network);
        igraph_bool_t all_zero = true;
        for (i = 0; i < total_edges; i++) {
            if (VECTOR(*flow)[i] != 0) {
                all_zero = false;
                break;
            }
        }
        IGRAPH_ASSERT(all_zero);
    }

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Public API                                                            */
/* -------------------------------------------------------------------- */

/**
 * \function igraph_minimum_path_cover
 * \brief Minimum path cover of a directed acyclic graph.
 *
 * A path cover of a DAG is a set of directed paths, each following real
 * edges of the graph, that together cover every vertex at least once.
 * Unlike a path \em partition, paths in a cover are allowed to share
 * vertices at merge points; this is the notion of "minimum path cover"
 * used in flow-decomposition applications (e.g. transcript assembly),
 * where it coincides with the minimum-flow-decomposition width. This
 * function computes a path cover of minimum cardinality (minimum
 * "width"), which by Dilworth's theorem equals the size of a maximum
 * antichain of the graph's reachability order.
 *
 * </para><para>
 * The computation proceeds in three stages: (1) an initial feasible
 * solution to an equivalent minimum-flow-with-lower-bounds problem is
 * built via vertex splitting, using the strategy selected by \p reduction;
 * (2) the flow is reduced to the true minimum, using the algorithm
 * selected by \p solver; (3) the minimum flow is decomposed back into
 * paths.
 *
 * \param graph The input graph. Must be directed and acyclic.
 * \param cover Pointer to an initialized list of integer vectors. The
 *        paths of the minimum path cover are stored here on return, one
 *        vector per path, listing the vertex IDs of that path in
 *        traversal order.
 * \param width Pointer to an integer, the number of paths (the width of
 *        the cover) is stored here, unless it is a null pointer.
 * \param vertex_weights Optional vector of non-negative integer vertex
 *        weights. A weight of \c k means the vertex must be covered by
 *        \c k separate paths (a weighted path cover). A null pointer
 *        means every vertex has weight 1, which is the usual case.
 *        Non-uniform weights are only supported together with
 *        \c IGRAPH_MPC_REDUCTION_NAIVE.
 * \param reduction The strategy used to build the initial feasible flow.
 * \param solver The algorithm used to reduce the flow to the minimum.
 * \param flow_network Optional pointer to an uninitialized graph. If not
 *        a null pointer, the solved (but not yet decomposed) auxiliary
 *        flow network is stored here, for later use with
 *        \ref igraph_maximum_antichain().
 * \param flow_network_flow Optional pointer to an initialized vector. If
 *        not a null pointer, filled with the flow on each edge of
 *        \p flow_network, indexed by edge ID.
 * \param flow_network_demand Optional pointer to an initialized vector.
 *        If not a null pointer, filled with the demand (lower bound) on
 *        each edge of \p flow_network, indexed by edge ID.
 * \return Error code.
 *
 * Time complexity: depends on \p reduction and \p solver. Building the
 * initial flow is O(|V|+|E|) for the naive reduction and O((|V|+|E|) w)
 * for the greedy reduction, where w is the resulting width. The naive DFS
 * solver is also O((|V|+|E|) w); the max-flow reduction solver is bounded
 * by the complexity of \ref igraph_maxflow(), O(|V|^3).
 *
 * </para><para>
 * This is a port of the Minimum Path Cover algorithms from
 * https://github.com/algbio/PerformanceMPC (Cáceres, Grigorjew et al.); see
 * the internal functions in path_cover.c for the exact literature reference
 * backing each reduction/solver combination.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Mäkinen V, Tomescu AI, Kuosmanen A, Paavilainen T, Gagie T, Chikhi R:
 * Sparse Dynamic Programming on DAGs with Small Width.
 * ACM Transactions on Algorithms 15(2):29, 2019.
 * https://doi.org/10.1145/3301312
 *
 * \sa \ref igraph_maximum_antichain(), \ref igraph_minimum_chain_cover(),
 * \ref igraph_maxflow(), \ref igraph_is_dag().
 */
igraph_error_t igraph_minimum_path_cover(
        const igraph_t *graph,
        igraph_vector_int_list_t *cover,
        igraph_int_t *width,
        const igraph_vector_int_t *vertex_weights,
        igraph_mpc_reduction_t reduction,
        igraph_mpc_solver_t solver,
        igraph_t *flow_network,
        igraph_vector_int_t *flow_network_flow,
        igraph_vector_int_t *flow_network_demand) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t total_edges = igraph_ecount(graph) + 3 * no_of_nodes;
    igraph_int_t source = 2 * no_of_nodes;
    igraph_int_t sink = 2 * no_of_nodes + 1;
    igraph_t network;
    igraph_vector_int_t demand, flow;
    igraph_bool_t is_dag;

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("Minimum path cover requires a directed graph.", IGRAPH_EINVAL);
    }
    IGRAPH_CHECK(igraph_is_dag(graph, &is_dag));
    if (!is_dag) {
        IGRAPH_ERROR("Minimum path cover requires an acyclic graph.", IGRAPH_EINVAL);
    }
    if (vertex_weights) {
        igraph_int_t i;
        if (igraph_vector_int_size(vertex_weights) != no_of_nodes) {
            IGRAPH_ERROR("Vertex weight vector length must match the number of vertices.",
                         IGRAPH_EINVAL);
        }
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(*vertex_weights)[i] < 0) {
                IGRAPH_ERROR("Vertex weights must be non-negative.", IGRAPH_EINVAL);
            }
        }
        if (reduction == IGRAPH_MPC_REDUCTION_GREEDY) {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*vertex_weights)[i] != 1) {
                    IGRAPH_ERROR("IGRAPH_MPC_REDUCTION_GREEDY only supports uniform "
                                 "vertex weights of 1; use IGRAPH_MPC_REDUCTION_NAIVE "
                                 "for a weighted path cover.", IGRAPH_UNIMPLEMENTED);
                }
            }
        }
    }

    IGRAPH_CHECK(igraph_i_mpc_build_network(graph, &network, &demand, vertex_weights));
    IGRAPH_FINALLY(igraph_destroy, &network);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &demand);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&flow, total_edges);
    switch (reduction) {
    case IGRAPH_MPC_REDUCTION_NAIVE:
        IGRAPH_CHECK(igraph_i_mpc_naive_reduction(graph, vertex_weights, &flow));
        break;
    case IGRAPH_MPC_REDUCTION_GREEDY:
        IGRAPH_CHECK(igraph_i_mpc_greedy_reduction(graph, &flow));
        break;
    default:
        IGRAPH_ERROR("Invalid reduction strategy.", IGRAPH_EINVAL);
    }

    {
        igraph_bool_t valid;
        IGRAPH_CHECK(igraph_i_mpc_is_valid_minflow(&network, &flow, &demand, &valid));
        IGRAPH_ASSERT(valid);
    }

    switch (solver) {
    case IGRAPH_MPC_SOLVER_NAIVE_DFS:
        IGRAPH_CHECK(igraph_i_mpc_solve_naive_dfs(&network, &demand, source, sink, &flow));
        break;
    case IGRAPH_MPC_SOLVER_MAXFLOW_REDUCTION:
        IGRAPH_CHECK(igraph_i_mpc_solve_maxflow_reduction(&network, &demand, source, sink, &flow));
        break;
    default:
        IGRAPH_ERROR("Invalid solver.", IGRAPH_EINVAL);
    }

    {
        igraph_bool_t valid;
        IGRAPH_CHECK(igraph_i_mpc_is_valid_minflow(&network, &flow, &demand, &valid));
        IGRAPH_ASSERT(valid);
    }

    if (width) {
        igraph_vector_int_t out_eids;
        igraph_int_t i, w = 0;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&out_eids, 0);
        IGRAPH_CHECK(igraph_incident(&network, &out_eids, source, IGRAPH_OUT, IGRAPH_LOOPS));
        for (i = 0; i < igraph_vector_int_size(&out_eids); i++) {
            w += VECTOR(flow)[VECTOR(out_eids)[i]];
        }
        igraph_vector_int_destroy(&out_eids);
        IGRAPH_FINALLY_CLEAN(1);
        *width = w;
    }

    if (flow_network_flow) {
        IGRAPH_CHECK(igraph_vector_int_update(flow_network_flow, &flow));
    }
    if (flow_network_demand) {
        IGRAPH_CHECK(igraph_vector_int_update(flow_network_demand, &demand));
    }
    if (flow_network) {
        IGRAPH_CHECK(igraph_copy(flow_network, &network));
    }

    igraph_vector_int_list_clear(cover);
    IGRAPH_CHECK(igraph_i_mpc_recover_paths(&network, &flow, source, sink, no_of_nodes, cover));

    igraph_vector_int_destroy(&flow);
    igraph_vector_int_destroy(&demand);
    igraph_destroy(&network);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
