/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

#include "core/interruption.h"

/* Johnson's cycle detection algorithm
 *
 * Based on the original implementation in:
 * Johnson DB: Finding all the elementary circuits of a directed graph.
 * SIAM J Comput 4(1):77-84.
 * https://epubs.siam.org/doi/10.1137/0204007
 */

/*
 * State of the cycle search algorithm. Storing all the state variables in a
 * single struct allows us to resume the algorithm from any point and yield the
 * cycles one by one in an iterator-like manner.
 */
struct igraph_simple_cycle_search_state_t {
    /* Number of vertices in the graph */
    igraph_integer_t N;

    /* Number of edges in the graph */
    igraph_integer_t NE;

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
    igraph_vector_bool_t v_blocked;

    /* Whether the graph is directed */
    igraph_bool_t directed;

    /* Boolean indicating whether the algorithm should stop searching for cycles
     */
    igraph_bool_t stop_search;
};

/**
 * A struct to store one cycle found by the algorithm.
 */
typedef struct {
    /* the vertices in the cycle */
    igraph_vector_int_list_t *vertices;
    /* the edges in the cycle */
    igraph_vector_int_list_t *edges;
} igraph_i_simple_cycle_results_t;

/**
 * \experimental
 *
 * The implementation of procedure UNBLOCK from Johnson's paper
 *
 */
static igraph_error_t igraph_i_simple_cycles_unblock(
        igraph_simple_cycle_search_state_t *state,
        igraph_integer_t u) {

    // TODO: introduce stack for w & neis in order to reduce the number of iterations.
    igraph_vector_int_t *neis;
    igraph_integer_t w;
    igraph_stack_int_t u_stack;
    IGRAPH_STACK_INT_INIT_FINALLY(&u_stack, 0);
    IGRAPH_CHECK(igraph_stack_int_push(&u_stack, u));

    while (igraph_stack_int_size(&u_stack) > 0) {
        igraph_integer_t current_u = igraph_stack_int_top(&u_stack);
        VECTOR(state->v_blocked)[current_u] = false;

        neis = igraph_adjlist_get(&state->B, current_u);
        igraph_bool_t recurse_deeper = false;
        while (!igraph_vector_int_empty(neis) && !recurse_deeper) {
            w = igraph_vector_int_pop_back(neis);
            if (VECTOR(state->v_blocked)[w]) {
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
 * \experimental
 *
 * The implementation of procedure CIRCUIT from Johnson's paper
 *
 * Arguments:
 *
 * \param state: local state object of the search
 * \param V: vertex to start the search from
 * \param result_callback: callback function to handle the found cycles
 * \param max_cycle_length Limit the maximum length of cycles to search
 * for. Negative for no limit
 * \param arg: argument to pass to the callback function
 */

static igraph_error_t igraph_i_simple_cycles_circuit(
        igraph_simple_cycle_search_state_t *state,
        igraph_integer_t V,
        igraph_simple_cycle_handler_t *result_callback,
        igraph_integer_t max_cycle_length, void *arg) {

    const igraph_vector_int_t *neighbors;
    const igraph_vector_int_t *incident_edges;
    igraph_integer_t num_neighbors;
    igraph_integer_t S = V;  // start
    igraph_integer_t E = -1; // edge start

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

        igraph_integer_t i0 = 0;
        if (recurse_deeper) {
            // stack v & e
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->vertex_stack, V));
            if (E >= 0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, E));
            }
            // printf("Pushing %" IGRAPH_PRId " to stack, stack size is %" IGRAPH_PRId
            //        ", result size is %" IGRAPH_PRId "\n",
            //        V, igraph_vector_int_size(&state->vertex_stack),
            //        igraph_vector_int_list_size(vertices));
            VECTOR(state->v_blocked)[V] = true;
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
        for (igraph_integer_t i = i0; i < num_neighbors; ++i) {
            igraph_integer_t W = VECTOR(*neighbors)[i];
            igraph_integer_t WE = VECTOR(*incident_edges)[i];

            if (W == S) {
                // need to unblock no matter whether we store the result or not (in
                // undirected case, we may not necessarily)
                local_found = true;

                if ((!state->directed &&
                igraph_vector_int_size(&state->edge_stack) == 1 &&
                igraph_vector_int_get(&state->edge_stack, 0) == WE)) {
                    // printf("Skipping cycle to prevent self-loop: \n");
                    continue;
                }

                // prevent duplicates in undirected graphs by forcing a direction for
                // the closing edge
                if ((!state->directed &&
                     igraph_vector_int_get(&state->edge_stack, 0) > WE)) {
                    // printf("Skipping cycle to prevent duplicates: \n");
                    continue;
                }

                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, WE));

                // output circuit composed of stack
                // printf("Found cycle with size %" IGRAPH_PRId "\n",
                //        igraph_vector_int_size(&state->vertex_stack));

                if ((*result_callback)(&state->vertex_stack, &state->edge_stack, arg) !=
                        IGRAPH_SUCCESS) {
                    state->stop_search = true;
                    break;
                }
                igraph_vector_int_pop_back(&state->edge_stack);
            } else if (!(VECTOR(state->v_blocked)[W])) {
                // printf("Recursing deeper from %" IGRAPH_PRId " to  %" IGRAPH_PRId
                // "\n",
                //        V, W);
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
                // IGRAPH_CHECK(igraph_i_simple_cycles_unblock_recursive(state, V));
                IGRAPH_CHECK(igraph_i_simple_cycles_unblock(state, V));
            } else {
                for (igraph_integer_t i = 0; i < num_neighbors; ++i) {
                    igraph_integer_t W = VECTOR(*neighbors)[i];
                    if (!igraph_vector_int_contains(igraph_adjlist_get(&state->B, W),
                                                    V)) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(
                                         igraph_adjlist_get(&state->B, W), V));
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

    IGRAPH_ASSERT(igraph_stack_int_size(&v_stack) == 0);
    IGRAPH_ASSERT(igraph_stack_int_size(&e_stack) == 0);
    IGRAPH_ASSERT(igraph_stack_int_size(&neigh_iteration_progress) == 0);
    igraph_stack_int_destroy(&neigh_iteration_progress);
    igraph_stack_int_destroy(&v_stack);
    igraph_stack_int_destroy(&e_stack);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycle_search_state_init
 * \brief Initialize the cycle search state
 *
 * \experimental
 *
 * \param state The state structure to initialize
 * \param graph The graph object
 * \return Error code.
 *
 * Time complexity: O(|V|*|E|*log(|V|*|E|))
 *
 * \ref igraph_simple_cycle_search_state_destroy
 */
igraph_error_t igraph_simple_cycle_search_state_init(
        igraph_simple_cycle_search_state_t *state,
        const igraph_t *graph) {

    state->N = igraph_vcount(graph);
    state->NE = igraph_ecount(graph);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&state->vertex_stack, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->vertex_stack, 8));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&state->edge_stack, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->edge_stack, 8));
    IGRAPH_CHECK(igraph_vector_bool_init(&state->v_blocked, state->N));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &state->v_blocked);
    IGRAPH_CHECK(igraph_inclist_init(
                     graph, &state->IK, IGRAPH_OUT,
                     IGRAPH_LOOPS_ONCE)); // each self-loop counts as a single cycle in
    // directed graphs
    IGRAPH_FINALLY(igraph_inclist_destroy, &state->IK);
    IGRAPH_CHECK(igraph_adjlist_init_from_inclist(graph, &state->AK, &state->IK));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &state->AK);
    state->directed = igraph_is_directed(graph);
    state->stop_search = false;
    IGRAPH_CHECK(igraph_adjlist_init_empty(&state->B, state->N));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &state->B);

    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycle_search_state_destroy
 * \brief Destroy the cycle search state
 *
 * \experimental
 *
 * \param state The state structure to destroy
 * \return Error code.
 *
 * Time complexity: O(1).
 *
 * \ref igraph_simple_cycle_search_state_init
 */
void igraph_simple_cycle_search_state_destroy(igraph_simple_cycle_search_state_t *state) {
    igraph_vector_int_destroy(&state->vertex_stack);
    igraph_vector_int_destroy(&state->edge_stack);
    igraph_vector_bool_destroy(&state->v_blocked);
    igraph_adjlist_destroy(&state->AK);
    igraph_inclist_destroy(&state->IK);
    igraph_adjlist_destroy(&state->B);
}

/**
 * \function igraph_i_append_simple_cycle_result
 * \brief Implements `igraph_simple_cycle_handler_t` to simply append the found cycle to a list
 *
 * \experimental
 *
 * \param vertices The vertices of the cycle
 * \param edges The edges of the cycle
 * \param arg The result list (pointer to `igraph_i_simple_cycle_results_t`) to append the cycle to
 *
 * \return IGRAPH_SUCCESS if the operation was successful, IGRAPH_STOP otherwise.
 * */
static igraph_error_t igraph_i_append_simple_cycle_result(
        const igraph_vector_int_t *vertices,
        const igraph_vector_int_t *edges,
        void *arg) {

    igraph_i_simple_cycle_results_t *res_list = (igraph_i_simple_cycle_results_t *) arg;
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

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycles_search_callback_from_one_vertex
 * \brief Search simple cycles starting from one vertex
 *
 * \experimental
 *
 * \param state The state structure to search on
 * \param s The vertex index to start search with
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 * Negative for no limit.
 * \param cycle_handler The callback function to call when a cycle is found.
 * See \ref igraph_simple_cycle_handler_t() for details.
 * \param arg The additional argument(s) for the callback function
 *
 * \return Error code.
 *
 * @see https://en.wikipedia.org/wiki/Johnson%27s_algorithm
 * @see https://stackoverflow.com/a/35922906/3909202
 * @see https://epubs.siam.org/doi/epdf/10.1137/0204007
 */
igraph_error_t igraph_simple_cycles_search_callback_from_one_vertex(
        igraph_simple_cycle_search_state_t *state,
        igraph_integer_t s,
        igraph_integer_t max_cycle_length,
        igraph_simple_cycle_handler_t *cycle_handler,
        void *arg) {

    // L3:
    for (igraph_integer_t i = s; i < state->N; ++i) {
        VECTOR(state->v_blocked)[i] = false;
        igraph_vector_int_clear(igraph_adjlist_get(&state->B, i));
    }

    IGRAPH_CHECK(igraph_i_simple_cycles_circuit(state, s, cycle_handler,
                 max_cycle_length, arg));

    for (igraph_integer_t i = 0; i < state->N; ++i) {
        // We want to remove the vertex with value s, not at position s.
        // It's fine to use binary search since we never add to, only remove from
        // an already sorted adjacency list.
        igraph_integer_t pos;
        if (igraph_vector_int_binsearch(igraph_adjlist_get(&state->AK, i), s,
                                        &pos)) {
            igraph_vector_int_remove(igraph_adjlist_get(&state->AK, i), pos);
            igraph_vector_int_remove(igraph_inclist_get(&state->IK, i), pos);
        }
    }
    igraph_vector_int_clear(igraph_adjlist_get(&state->AK, s));
    igraph_vector_int_clear(igraph_inclist_get(&state->IK, s));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycles_search_from_one_vertex
 * \brief Search simple cycles starting from one vertex
 *
 * \experimental
 *
 * \param state The state structure to search on
 * \param s The vertex index to start search with
 * \param vertices The vertex IDs of each cycle will be stored here
 * \param edges The edge IDs of each cycle will be stored here
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 * Negative for no limit.
 *
 * \return Error code.
 *
 * See also: \ref igraph_simple_cycles_search_callback_from_one_vertex()
 *
 * @see https://en.wikipedia.org/wiki/Johnson%27s_algorithm
 * @see https://stackoverflow.com/a/35922906/3909202
 * @see https://epubs.siam.org/doi/epdf/10.1137/0204007
 */
igraph_error_t igraph_simple_cycles_search_from_one_vertex(
        igraph_simple_cycle_search_state_t *state,
        igraph_integer_t s,
        igraph_vector_int_list_t *vertices,
        igraph_vector_int_list_t *edges,
        igraph_integer_t max_cycle_length) {

    igraph_i_simple_cycle_results_t result_list;
    result_list.vertices = vertices;
    result_list.edges = edges;

    igraph_simple_cycles_search_callback_from_one_vertex(
        state, s, max_cycle_length, &igraph_i_append_simple_cycle_result, &result_list);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycles_search_callback
 * \brief Search all simple cycles, use a callback function to handle found cycles
 *
 * \experimental
 *
 * This function searches for all simple cycles,
 * using Johnson's cycle detection algorithm
 * based on the original implementation in:
 * Johnson DB: Finding all the elementary circuits of a directed graph.
 * SIAM J Comput 4(1):77-84.
 * https://epubs.siam.org/doi/10.1137/0204007
 *
 * See also: \ref igraph_simple_cycles_search_callback()
 *
 * \param graph The graph to search for
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 * Negative for no limit.
 * \param cycle_handler The callback function to call when a cycle is found.
 * See \ref igraph_simple_cycle_handler_t() for details.
 * \param arg The additional argument(s) for the callback function
 *
 * \return Error code.
 */
igraph_error_t igraph_simple_cycles_search_callback(
        const igraph_t *graph,
        igraph_integer_t max_cycle_length,
        igraph_simple_cycle_handler_t *cycle_handler,
        void *arg) {

    if (max_cycle_length == 0) {
        return IGRAPH_SUCCESS;
    }

    igraph_simple_cycle_search_state_t state;

    IGRAPH_CHECK(igraph_simple_cycle_search_state_init(&state, graph));
    IGRAPH_FINALLY(igraph_simple_cycle_search_state_destroy, &state);

    // TODO: depending on the graph, it is rather unreasonable to search cycles
    // from each and every node
    for (igraph_integer_t i = 0; i < state.N; i++) {
        if (!igraph_vector_int_empty(igraph_adjlist_get(&state.AK, i))) {
            IGRAPH_CHECK(igraph_simple_cycles_search_callback_from_one_vertex(
                             &state, i, max_cycle_length, cycle_handler, arg));
            IGRAPH_ALLOW_INTERRUPTION();
        }
        if (state.stop_search) {
            state.stop_search = false;
            break;
        }
    }

    igraph_simple_cycle_search_state_destroy(&state);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_simple_cycles_search_all
 * \brief Search all simple cycles
 *
 * \experimental
 *
 * This function searches for all simple cycles,
 * using Johnson's cycle detection algorithm
 * based on the original implementation in:
 * Johnson DB: Finding all the elementary circuits of a directed graph.
 * SIAM J Comput 4(1):77-84.
 * https://epubs.siam.org/doi/10.1137/0204007
 *
 *
 * \param graph The graph to search for
 * \param vertices The vertex IDs of each cycle will be stored here
 * \param edges The edge IDs of each cycle will be stored here
 * \param max_cycle_length Limit the maximum length of cycles to search for.
 * Negative for no limit.
 *
 * \return Error code.
 */
igraph_error_t igraph_simple_cycles_search_all(
        const igraph_t *graph,
        igraph_vector_int_list_t *v_result,
        igraph_vector_int_list_t *e_result,
        igraph_integer_t max_cycle_length) {

    igraph_i_simple_cycle_results_t result_list;
    result_list.vertices = v_result;
    result_list.edges = e_result;

    igraph_simple_cycles_search_callback(graph, max_cycle_length,
                                         &igraph_i_append_simple_cycle_result,
                                         &result_list);

    return IGRAPH_SUCCESS;
}
