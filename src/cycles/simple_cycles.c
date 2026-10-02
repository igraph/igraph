/*
   igraph library.
   Copyright (C) 2024-2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_cycles.h"

#include "igraph_adjlist.h"
#include "igraph_bitset.h"
#include "igraph_components.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

#include "core/interruption.h"

/* Johnson's cycle detection algorithm
 *
 * Based on the original implementation in:
 * Johnson DB: Finding all the elementary circuits of a directed graph.
 * SIAM J Comput 4(1):77-84.
 * https://doi.org/10.1137/0204007
 */

/**
 * State of the cycle search algorithm. Storing all the state variables in a
 * single struct allows us to resume the algorithm from any point and yield the
 * cycles one by one in an iterator-like manner.
 */
typedef struct {
    /* Number of vertices in the graph */
    igraph_int_t N;

    /* The incidence list of the graph */
    igraph_inclist_t IK;

    /* The adjacency list of the graph */
    igraph_adjlist_t AK;

    /* The set B from Johnson's paper */
    igraph_adjlist_t B;

    /* Stack in which the vertices of the current cycle are pushed */
    igraph_vector_int_t vertex_stack;

    /* Stack in which the edges of the current cycle are pushed */
    igraph_vector_int_t edge_stack;

    /* Boolean vector indicating which vertices are blocked */
    igraph_bitset_t v_blocked;

    /* Whether the graph is directed */
    igraph_bool_t directed;

    /* Boolean indicating whether the algorithm should stop searching for cycles */
    igraph_bool_t stop_search;
} simple_cycle_search_state_t;

/**
 * A struct to store one cycle found by the algorithm.
 */
typedef struct {
    /* the vertices in the cycle */
    igraph_vector_int_list_t *vertices;
    /* the edges in the cycle */
    igraph_vector_int_list_t *edges;
    /* number of cycles found so far */
    igraph_int_t cycle_count;
    /* how many cycles to record at most? */
    igraph_int_t max_results;
} simple_cycle_results_t;

/**
 * The implementation of procedure UNBLOCK from Johnson's paper
 */
static igraph_error_t simple_cycles_unblock(
    simple_cycle_search_state_t *state,
    igraph_int_t u) {

    // TODO: introduce stack for w & neis in order to reduce the number of
    // iterations.
    igraph_vector_int_t *neis;
    igraph_stack_int_t u_stack;

    IGRAPH_STACK_INT_INIT_FINALLY(&u_stack, 0);
    IGRAPH_CHECK(igraph_stack_int_push(&u_stack, u));

    while (igraph_stack_int_size(&u_stack) > 0) {
        igraph_bool_t recurse_deeper = false;
        const igraph_int_t current_u = igraph_stack_int_top(&u_stack);

        IGRAPH_BIT_CLEAR(state->v_blocked, current_u);

        neis = igraph_adjlist_get(&state->B, current_u);
        while (!igraph_vector_int_empty(neis) && !recurse_deeper) {
            const igraph_int_t w = igraph_vector_int_pop_back(neis);
            if (IGRAPH_BIT_TEST(state->v_blocked, w)) {
                IGRAPH_CHECK(igraph_stack_int_push(&u_stack, w));
                recurse_deeper = true;
            }
        }

        if (!recurse_deeper) {
            igraph_stack_int_pop(&u_stack);
        }
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_stack_int_destroy(&u_stack);

    return IGRAPH_SUCCESS;
}

/**
 * The implementation of procedure CIRCUIT from Johnson's paper
 *
 * \param state Local state object of the search.
 * \param V Vertex to start the search from.
 * \param callback Callback function to handle the found cycles.
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 *   Pass a negative value for no limit
 * \param arg Argument to pass to the callback function.
 */
static igraph_error_t simple_cycles_circuit(
    simple_cycle_search_state_t *state,
    igraph_int_t V,
    const igraph_bitset_t *in_component,
    igraph_int_t max_cycle_length,
    igraph_int_t min_cycle_length,
    igraph_cycle_handler_t *callback,
    void *arg) {

    const igraph_vector_int_t *neighbors;
    const igraph_vector_int_t *incident_edges;
    igraph_int_t num_neighbors;
    igraph_int_t S = V;  // start
    igraph_int_t E = -1; // edge start

    igraph_bool_t local_found = false;
    igraph_bool_t loop_length_stop = false;

    // keep track of what we were doing (rather than recursing)
    igraph_stack_int_t neigh_iteration_progress;
    IGRAPH_STACK_INT_INIT_FINALLY(&neigh_iteration_progress, 0);

    igraph_stack_int_t v_stack;
    IGRAPH_STACK_INT_INIT_FINALLY(&v_stack, 0);

    igraph_stack_int_t e_stack;
    IGRAPH_STACK_INT_INIT_FINALLY(&e_stack, 0);

    igraph_bool_t recurse_deeper = true;
    while ((recurse_deeper ||
            igraph_stack_int_size(&neigh_iteration_progress) > 0) &&
            !state->stop_search) {

        IGRAPH_ASSERT(igraph_stack_int_size(&neigh_iteration_progress) ==
                      igraph_stack_int_size(&e_stack));
        IGRAPH_ASSERT(igraph_stack_int_size(&v_stack) ==
                      igraph_stack_int_size(&e_stack));

        igraph_int_t i0 = 0;
        if (recurse_deeper) {
            // stack v & e
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->vertex_stack, V));
            if (E >= 0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, E));
            }
            // printf("Pushing %" IGRAPH_PRId " to stack, stack size is %" IGRAPH_PRId
            //        ", result size is %" IGRAPH_PRId "\n",
            //        V, igraph_vector_int_size(&state->vertex_stack),
            //        igraph_stack_int_size(&v_stack));
            IGRAPH_BIT_SET(state->v_blocked, V);
        } else {
            // back to what we were doing before
            i0 = igraph_stack_int_pop(&neigh_iteration_progress);
            V = igraph_stack_int_pop(&v_stack);
            E = igraph_stack_int_pop(&e_stack);
        }
        recurse_deeper = false;

        // L1
        neighbors = igraph_adjlist_get(&state->AK, V);
        incident_edges = igraph_inclist_get(&state->IK, V);
        num_neighbors = igraph_vector_int_size(neighbors);
        IGRAPH_ASSERT(igraph_vector_int_size(incident_edges) == num_neighbors);
        for (igraph_int_t i = i0; i < num_neighbors; ++i) {
            igraph_int_t W = VECTOR(*neighbors)[i];
            igraph_int_t WE = VECTOR(*incident_edges)[i];

            if (!IGRAPH_BIT_TEST(*in_component, W)) {
                continue;
            }

            if (W == S) {
                igraph_error_t ret;

                // need to unblock no matter whether we store the result or not (in
                // undirected case, we may not necessarily)
                local_found = true;

                if ((!state->directed &&
                        igraph_vector_int_size(&state->edge_stack) == 1 &&
                        VECTOR(state->edge_stack)[0] == WE)) {
                    // printf("Skipping cycle to %" IGRAPH_PRId " via %" IGRAPH_PRId " to prevent self-loop.\n", W, WE);
                    continue;
                }

                // prevent duplicates in undirected graphs by forcing a direction for
                // the closing edge
                if ((!state->directed &&
                        igraph_vector_int_size(&state->edge_stack) > 0 &&
                        VECTOR(state->edge_stack)[0] > WE)) {
                    // printf("Skipping cycle to %" IGRAPH_PRId " via %" IGRAPH_PRId " to
                    // prevent duplicates.\n", W, WE);
                    continue;
                }

                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, WE));

                // output circuit composed of stack
                // printf("Found cycle with size %" IGRAPH_PRId "\n",
                //        igraph_vector_int_size(&state->vertex_stack));

                if (igraph_vector_int_size(&state->edge_stack) >= min_cycle_length) {
                    IGRAPH_CHECK_CALLBACK(
                        callback(&state->vertex_stack, &state->edge_stack, arg), &ret);
                    if (ret == IGRAPH_STOP) {
                        state->stop_search = true;
                        break;
                    }
                }
                igraph_vector_int_pop_back(&state->edge_stack);
            } else if (! IGRAPH_BIT_TEST(state->v_blocked, W)) {
                // printf("Recursing deeper from %" IGRAPH_PRId " to  %" IGRAPH_PRId
                // "\n", V, W);
                recurse_deeper = ((max_cycle_length < 0) ||
                                  (igraph_vector_int_size(&state->vertex_stack) <=
                                   max_cycle_length - 1));
                if (recurse_deeper) {
                    IGRAPH_CHECK(igraph_stack_int_push(&neigh_iteration_progress, i + 1));
                    IGRAPH_CHECK(igraph_stack_int_push(&v_stack, V));
                    IGRAPH_CHECK(igraph_stack_int_push(&e_stack, E));
                    V = W;
                    E = WE;
                    break;
                } else {
                    loop_length_stop = true;
                }
            } else {
                // printf("Vertex %" IGRAPH_PRId ", neighbour of %" IGRAPH_PRId
                //        ", is blocked\n",
                //        W, V);
            }
        }

        if (!recurse_deeper) {
            // L2
            if (local_found || loop_length_stop) {
                IGRAPH_CHECK(simple_cycles_unblock(state, V));
            } else {
                for (igraph_int_t i = 0; i < num_neighbors; ++i) {
                    const igraph_int_t W = VECTOR(*neighbors)[i];
                    if (!IGRAPH_BIT_TEST(*in_component, W)) {
                        continue;
                    }
                    if (!igraph_vector_int_contains(igraph_adjlist_get(&state->B, W), V)) {
                        IGRAPH_CHECK(
                            igraph_vector_int_push_back(igraph_adjlist_get(&state->B, W), V));
                    }
                }
            }

            IGRAPH_ASSERT(!igraph_vector_int_empty(&state->vertex_stack));

            // unstack v
            V = igraph_vector_int_pop_back(&state->vertex_stack);
            if (!igraph_vector_int_empty(&state->edge_stack)) {
                // can be empty for the starting point.
                // alternatively, V == S
                E = igraph_vector_int_pop_back(&state->edge_stack);
            }
            // printf("Unstacked v %" IGRAPH_PRId ", e %" IGRAPH_PRId "\n", V, E);
        }

        IGRAPH_ALLOW_INTERRUPTION();
    }

    if (! state->stop_search) {
        IGRAPH_ASSERT(igraph_stack_int_size(&v_stack) == 0);
        IGRAPH_ASSERT(igraph_stack_int_size(&e_stack) == 0);
        IGRAPH_ASSERT(igraph_stack_int_size(&neigh_iteration_progress) == 0);
    }

    igraph_stack_int_destroy(&neigh_iteration_progress);
    igraph_stack_int_destroy(&v_stack);
    igraph_stack_int_destroy(&e_stack);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function simple_cycle_search_state_init
 * \brief Initializes the cycle search state.
 *
 * \param state The state structure to initialize.
 * \param graph The graph object.
 * \param mode A constant specifying how edge directions are
 *   considered in directed graphs. Valid modes are:
 *   \c IGRAPH_OUT, follows edge directions;
 *   \c IGRAPH_IN, follows the opposite directions; and
 *   \c IGRAPH_ALL, ignores edge directions. This argument is
 *   ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V|*|E|*log(|V|*|E|))
 *
 * \ref simple_cycle_search_state_destroy()
 */
static igraph_error_t simple_cycle_search_state_init(
    simple_cycle_search_state_t *state,
    const igraph_t *graph,
    igraph_neimode_t mode) {

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for finding cycles.", IGRAPH_EINVAL);
    }

    state->N = igraph_vcount(graph);
    state->directed = igraph_is_directed(graph) && mode != IGRAPH_ALL;
    state->stop_search = false;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&state->vertex_stack, 8);
    igraph_vector_int_clear(&state->vertex_stack);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&state->edge_stack, 8);
    igraph_vector_int_clear(&state->edge_stack);

    // by default false
    IGRAPH_BITSET_INIT_FINALLY(&state->v_blocked, state->N);

    IGRAPH_CHECK(igraph_inclist_init(
                     graph, &state->IK, mode,
                     IGRAPH_LOOPS_ONCE // each self-loop counts as a single cycle
                 ));
    IGRAPH_FINALLY(igraph_inclist_destroy, &state->IK);

    IGRAPH_CHECK(igraph_adjlist_init_from_inclist(graph, &state->AK, &state->IK));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &state->AK);

    IGRAPH_CHECK(igraph_adjlist_init_empty(&state->B, state->N));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &state->B);

    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/**
 * \function simple_cycle_search_state_destroy
 * \brief Destroys the cycle search state.
 *
 * \param state The state structure to destroy
 * \return Error code.
 *
 * Time complexity: O(1).
 *
 * \ref simple_cycle_search_state_init()
 */
static void simple_cycle_search_state_destroy(simple_cycle_search_state_t *state) {

    igraph_adjlist_destroy(&state->B);
    igraph_adjlist_destroy(&state->AK);
    igraph_inclist_destroy(&state->IK);
    igraph_bitset_destroy(&state->v_blocked);
    igraph_vector_int_destroy(&state->edge_stack);
    igraph_vector_int_destroy(&state->vertex_stack);
}

/**
 * A cycle handler that simply appends cycles to a vector list.
 * Use by \ref igraph_simple_cycles()
 */
static igraph_error_t append_simple_cycle_result(
    const igraph_vector_int_t *vertices,
    const igraph_vector_int_t *edges,
    void *arg) {

    simple_cycle_results_t *res_list = (simple_cycle_results_t *) arg;

    if (res_list->vertices != NULL) {
        // copy output: from stack to vector. No need to reverse because
        // we were putting vertices in the stack in reverse order anyway.
        igraph_vector_int_t v_res;
        IGRAPH_CHECK(igraph_vector_int_init_copy(&v_res, vertices));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &v_res);
        /* v_res ownership transferred to 'vertices' */
        IGRAPH_CHECK(igraph_vector_int_list_push_back(res_list->vertices, &v_res));
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (res_list->edges != NULL) {
        // same for edges
        igraph_vector_int_t e_res;
        IGRAPH_CHECK(igraph_vector_int_init_copy(&e_res, edges));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &e_res);
        /* e_res ownership transferred to 'edges' */
        IGRAPH_CHECK(igraph_vector_int_list_push_back(res_list->edges, &e_res));
        IGRAPH_FINALLY_CLEAN(1);
    }

    res_list->cycle_count++;

    if (res_list->max_results >= 0 && res_list->cycle_count == res_list->max_results) {
        return IGRAPH_STOP;
    } else {
        return IGRAPH_SUCCESS;
    }
}

/**
 * \function simple_cycles_search_callback_from_one_vertex
 * \brief Search simple cycles starting from one vertex.
 *
 * \param state The state structure to search on.
 * \param s The vertex index to start search with.
 * \param in_component A bitset indicating which vertices are in the same component as s.
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 *   Pass a negative value for no limit.
 * \param callback The callback function to call when a cycle is found.
 *   See \ref igraph_cycle_handler_t() for details.
 * \param arg The additional argument(s) for the callback function.
 *
 * \return Error code.
 *
 * https://en.wikipedia.org/wiki/Johnson%27s_algorithm
 * https://stackoverflow.com/a/35922906/3909202
 * https://epubs.siam.org/doi/epdf/10.1137/0204007
 */
static igraph_error_t simple_cycles_search_callback_from_one_vertex_in_component(
    simple_cycle_search_state_t *state,
    igraph_int_t s,
    const igraph_bitset_t *in_component,
    igraph_int_t min_cycle_length,
    igraph_int_t max_cycle_length,
    igraph_cycle_handler_t *callback,
    void *arg) {

    for (igraph_int_t i = s; i < state->N; ++i) {
        IGRAPH_BIT_CLEAR(state->v_blocked, i);
        igraph_vector_int_clear(igraph_adjlist_get(&state->B, i));
    }

    IGRAPH_CHECK(simple_cycles_circuit(state, s, in_component, max_cycle_length,
                                       min_cycle_length, callback, arg));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycles_callback
 * \brief Finds all simple cycles (callback version).
 *
 * \experimental
 *
 * This function searches for all simple cycles using Johnson's cycle
 * detection algorithm, and calls a function for each.
 * A simple cycle is a cycle (i.e. closed path) without repeated vertices.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Johnson DB: Finding all the elementary circuits of a directed graph.
 * SIAM J. Comput. 4(1):77-84.
 * https://doi.org/10.1137/0204007
 *
 * \param graph The graph to search for
 * \param mode A constant specifying how edge directions are
 *   considered in directed graphs. Valid modes are:
 *   \c IGRAPH_OUT, follows edge directions;
 *   \c IGRAPH_IN, follows the opposite directions; and
 *   \c IGRAPH_ALL, ignores edge directions. This argument is
 *   ignored for undirected graphs.
 * \param min_cycle_length Limit the minimum length of cycles to search for.
 *   Pass a negative value to search for all cycles.
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 *   Pass a negative value to search for all cycles.
 * \param callback A function to call for each cycle that is found.
 *   See also \ref igraph_cycle_handler_t
 * \param arg This parameter will be passed to \p callback.
 * \return Error code.
 *
 * \sa \ref igraph_simple_cycles() to store the found cycles;
 * \ref igraph_find_cycle() to find a single cycle;
 * \ref igraph_fundamental_cycles() and igraph_minimum_cycle_basis()
 * to find a cycle basis, a compact representation of the cycle structure
 * of the graph.
 */
igraph_error_t igraph_simple_cycles_callback(
    const igraph_t *graph,
    igraph_neimode_t mode,
    igraph_int_t min_cycle_length,
    igraph_int_t max_cycle_length,
    igraph_cycle_handler_t *callback,
    void *arg) {

    if (max_cycle_length == 0) {
        return IGRAPH_SUCCESS;
    }

    simple_cycle_search_state_t state;
    IGRAPH_CHECK(simple_cycle_search_state_init(&state, graph, mode));
    IGRAPH_FINALLY(simple_cycle_search_state_destroy, &state);

    igraph_t working_graph;
    IGRAPH_CHECK(igraph_copy(&working_graph, graph));
    IGRAPH_FINALLY(igraph_destroy, &working_graph);

    igraph_vector_int_t membership;
    igraph_vector_int_t csize;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&csize, 0);

    igraph_bitset_t removed;
    IGRAPH_BITSET_INIT_FINALLY(&removed, state.N);

    igraph_bitset_t in_component;
    IGRAPH_BITSET_INIT_FINALLY(&in_component, state.N);

    igraph_connectedness_t scc_mode = state.directed ? IGRAPH_STRONG : IGRAPH_WEAK;

    while (!state.stop_search) {
        igraph_int_t no_comps;
        IGRAPH_CHECK(igraph_connected_components(&working_graph, &membership, &csize, &no_comps, scc_mode));

        if (!state.directed) {
            igraph_vector_int_t comp_edges;
            igraph_vector_int_t edges_to_delete;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&comp_edges, no_comps);
            IGRAPH_VECTOR_INT_INIT_FINALLY(&edges_to_delete, 0);

            igraph_int_t no_edges = igraph_ecount(&working_graph);
            for (igraph_int_t e = 0; e < no_edges; e++) {
                igraph_int_t u, v;
                igraph_edge(&working_graph, e, &u, &v);
                igraph_int_t c_id = VECTOR(membership)[u];
                VECTOR(comp_edges)[c_id]++;
            }

            for (igraph_int_t e = 0; e < no_edges; e++) {
                igraph_int_t u, v;
                igraph_edge(&working_graph, e, &u, &v);
                igraph_int_t c_id = VECTOR(membership)[u];
                if (VECTOR(comp_edges)[c_id] < VECTOR(csize)[c_id]) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges_to_delete, e));
                }
            }

            if (igraph_vector_int_size(&edges_to_delete) > 0) {
                IGRAPH_CHECK(igraph_delete_edges(&working_graph, igraph_ess_vector(&edges_to_delete)));
            }

            for (igraph_int_t i = 0; i < state.N; i++) {
                if (IGRAPH_BIT_TEST(removed, i)) {
                    continue;
                }
                igraph_int_t c_id = VECTOR(membership)[i];
                if (VECTOR(comp_edges)[c_id] < VECTOR(csize)[c_id]) {
                    IGRAPH_BIT_SET(removed, i);
                }
            }

            igraph_vector_int_destroy(&edges_to_delete);
            igraph_vector_int_destroy(&comp_edges);
            IGRAPH_FINALLY_CLEAN(2);

            IGRAPH_CHECK(igraph_connected_components(&working_graph, &membership, &csize, &no_comps, scc_mode));
        }

        // Find the first vertex in a non-trivial component
        igraph_int_t s = -1;
        igraph_int_t target_comp = -1;
        for (igraph_int_t i = 0; i < state.N; i++) {
            if (IGRAPH_BIT_TEST(removed, i)) {
                continue;
            }
            igraph_int_t c_id = VECTOR(membership)[i];
            igraph_int_t c_size = VECTOR(csize)[c_id];
            if (c_size >= 2) {
                s = i;
                target_comp = c_id;
                break;
            } else {
                // Check if it has a self-loop
                igraph_int_t pos;
                if (igraph_vector_int_binsearch(igraph_adjlist_get(&state.AK, i), i, &pos)) {
                    s = i;
                    target_comp = c_id;
                    break;
                }
            }
        }

        if (s == -1) {
            break;
        }

        // Set in_component bitset for the vertices in this component
        igraph_bitset_null(&in_component);
        for (igraph_int_t i = 0; i < state.N; i++) {
            if (!IGRAPH_BIT_TEST(removed, i) && VECTOR(membership)[i] == target_comp) {
                IGRAPH_BIT_SET(in_component, i);
            }
        }

        IGRAPH_CHECK(simple_cycles_search_callback_from_one_vertex_in_component(
                         &state, s, &in_component, min_cycle_length, max_cycle_length, callback, arg));

        // "Remove" s from the working graph by deleting its incident edges
        IGRAPH_BIT_SET(removed, s);
        igraph_vector_int_t incident_edges;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&incident_edges, 0);
        IGRAPH_CHECK(igraph_incident(&working_graph, &incident_edges, s, IGRAPH_ALL, IGRAPH_LOOPS));
        IGRAPH_CHECK(igraph_delete_edges(&working_graph, igraph_ess_vector(&incident_edges)));
        igraph_vector_int_destroy(&incident_edges);
        IGRAPH_FINALLY_CLEAN(1);

        IGRAPH_ALLOW_INTERRUPTION();
    }

    igraph_bitset_destroy(&in_component);
    igraph_bitset_destroy(&removed);
    igraph_vector_int_destroy(&csize);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&working_graph);
    simple_cycle_search_state_destroy(&state);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycles
 * \brief Finds all simple cycles.
 *
 * \experimental
 *
 * This function searches for all simple cycles using Johnson's cycle
 * detection algorithm, and stores them in the provided vector lists.
 * A simple cycle is a cycle (i.e. closed path) without repeated vertices.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Johnson DB: Finding all the elementary circuits of a directed graph.
 * SIAM J. Comput. 4(1):77-84.
 * https://doi.org/10.1137/0204007
 *
 * \param graph The graph to search for cycles in.
 * \param vertices The vertex IDs of each cycle will be stored here.
 * \param edges The edge IDs of each cycle will be stored here.
 * \param mode A constant specifying how edge directions are
 *   considered in directed graphs. Valid modes are:
 *   \c IGRAPH_OUT, follows edge directions;
 *   \c IGRAPH_IN, follows the opposite directions; and
 *   \c IGRAPH_ALL, ignores edge directions. This argument is
 *   ignored for undirected graphs.
 * \param min_cycle_length Limit the minimum length of cycles to search for.
 *   Pass a negative value to search for all cycles.
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 *   Pass a negative value to search for all cycles.
 * \param max_results At most this many cycles will be recorded. If
 *   negative, or \ref IGRAPH_UNLIMITED, no limit is applied.
 * \return Error code.
 *
 * \sa \ref igraph_simple_cycles_callback() to call a function for each found
 * cycle;
 * \ref igraph_find_cycle() to find a single cycle;
 * \ref igraph_fundamental_cycles() and \ref igraph_minimum_cycle_basis()
 * to find a cycle basis, a compact representation of the cycle structure
 * of the graph.
 */
igraph_error_t igraph_simple_cycles(
    const igraph_t *graph,
    igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
    igraph_neimode_t mode,
    igraph_int_t min_cycle_length, igraph_int_t max_cycle_length,
    igraph_int_t max_results) {

    simple_cycle_results_t result_list;
    result_list.vertices = vertices;
    result_list.edges = edges;
    result_list.cycle_count = 0;
    result_list.max_results = max_results;

    if (vertices) {
        igraph_vector_int_list_clear(vertices);
    }
    if (edges) {
        igraph_vector_int_list_clear(edges);
    }

    if (max_results == 0) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_simple_cycles_callback(graph, mode, min_cycle_length, max_cycle_length,
                 &append_simple_cycle_result,
                 &result_list));

    return IGRAPH_SUCCESS;
}
