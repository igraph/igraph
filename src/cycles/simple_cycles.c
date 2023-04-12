/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_stack.h"
#include "igraph_structural.h"
#include "../../tests/unit/test_utilities.h"

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
typedef struct igraph_simple_cycle_search_state_t {
    /* Number of vertices in the graph */
    igraph_integer_t N;

    /* The incidence list of the graph */
    igraph_inclist_t IK;

    /* The adjacency list of the graph */
    igraph_adjlist_t AK;

    /* The set B from Johnson's paper */
    igraph_adjlist_t B;

    /* Stack in which the vertices of the current cycle are pushed */
    igraph_vector_int_t v_stack;

    /* Stack in which the edges of the current cycle are pushed */
    igraph_vector_int_t e_stack;

    /* Boolean vector indicating which vertices are blocked */
    igraph_vector_bool_t blocked;

    /* Whether the graph is directed */
    igraph_bool_t directed;
} igraph_simple_cycle_search_state_t;

/**
 * \experimental
 *
 * A custom function to prevent double results in the search when using undirected graphs
 */

static igraph_bool_t igraph_i_cycle_has_been_found_already(
    igraph_simple_cycle_search_state_t *state,
    igraph_vector_int_list_t *v_results, igraph_vector_int_list_t *e_results,
    igraph_vector_int_t *v_res, igraph_vector_int_t *e_res
) {
    if (igraph_vector_int_size(v_res) < 2) {
        return false;
    }
    // TODO: other things to compare
    return VECTOR(*v_res)[0] < VECTOR(*v_res)[1];
}

/**
 * \experimental
 *
 * The implementation of procedure UNBLOCK from Johnson's paper
 *
 */

static igraph_error_t igraph_i_simple_cycles_unblock(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t u
) {
    igraph_vector_int_t* neis;
    igraph_integer_t w;

    VECTOR(state->blocked)[u] = false;

    neis = igraph_adjlist_get(&state->B, u);
    while (!igraph_vector_int_empty(neis)) {
        w = igraph_vector_int_pop_back(neis);
        if (VECTOR(state->blocked)[w]) {
            IGRAPH_CHECK(igraph_i_simple_cycles_unblock(state, w));
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \experimental
 *
 * The implementation of procedure CIRCUIT from Johnson's paper
 */

static igraph_error_t igraph_i_simple_cycles_circuit(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t V,igraph_integer_t E,
    igraph_integer_t S, igraph_vector_int_list_t *v_results, igraph_vector_int_list_t *e_results, bool *found,
    igraph_simple_cycle_search_mode_t search_mode
) {
    bool local_found = false;
    igraph_vector_int_t* neighbors;
    igraph_vector_int_t* neighbors_e;
    igraph_integer_t num_neighbors;

    // stack v & e
    IGRAPH_CHECK(igraph_vector_int_push_back(&state->v_stack, V));
    if (E >= 0) {
        IGRAPH_CHECK(igraph_vector_int_push_back(&state->e_stack, E));
    }
    // printf("Pushing %lld to stack, stack size is %lld, result size is %lld\n", V, igraph_vector_int_size(&state->v_stack), igraph_vector_int_list_size(results));
    VECTOR(state->blocked)[V] = true;

    // L1
    neighbors = igraph_adjlist_get(&state->AK, V);
    neighbors_e = igraph_inclist_get(&state->IK, V);
    num_neighbors = igraph_vector_int_size(neighbors);
    IGRAPH_ASSERT(igraph_vector_int_size(neighbors_e) == num_neighbors);
    for (igraph_integer_t i = 0; i < num_neighbors; ++i) {
        igraph_integer_t W = VECTOR(*neighbors)[i];
        igraph_integer_t WE = VECTOR(*neighbors_e)[i];
        // NOTE: possibly dangerous fix for undirected graphs,
        // disabling finding any two-vertex-loops
        if (W == S) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->e_stack, WE));
            if ((state->directed || igraph_vector_int_size(&state->v_stack) > 2)) {
                local_found = true;
                // output circuit composed of stack
                // printf("Found cycle with size %" IGRAPH_PRId "\n", igraph_vector_int_size(&state->v_stack));

                // copy output: from stack to vector
                igraph_vector_int_t v_res;
                IGRAPH_CHECK(igraph_vector_int_init_copy(&v_res, &state->v_stack));
                IGRAPH_FINALLY(igraph_vector_int_destroy, &v_res);
                IGRAPH_CHECK(igraph_vector_int_reverse(&v_res));
                // same for edges
                igraph_vector_int_t e_res;
                IGRAPH_CHECK(igraph_vector_int_init_copy(&e_res, &state->e_stack));
                IGRAPH_FINALLY(igraph_vector_int_destroy, &e_res);
                IGRAPH_CHECK(igraph_vector_int_reverse(&e_res));
                // undirected graphs lead to every cycle being found twice.
                // this is our naïve filter for now
                igraph_bool_t persist_result = true;
                if (!state->directed && search_mode == IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE) {
                    print_vector_int(&v_res);
                    igraph_bool_t duplicate_found = igraph_i_cycle_has_been_found_already(state, v_results, e_results, &v_res, &e_res);
                    if (duplicate_found) {
                        persist_result = false;
                    }
                }
                // end filter
                if (persist_result) {
                    IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(v_results, &v_res));
                    if (e_results != NULL) {
                        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(e_results, &e_res));
                    }
                }
                igraph_vector_int_destroy(&v_res);
                igraph_vector_int_destroy(&e_res);
                IGRAPH_FINALLY_CLEAN(2);
            }
        } else if (!(VECTOR(state->blocked)[W])) {
            IGRAPH_CHECK(igraph_i_simple_cycles_circuit(state, W, WE, S, v_results, e_results, &local_found, search_mode));
        }
    }
    *found = local_found;

    // L2
    if (local_found) {
        IGRAPH_CHECK(igraph_i_simple_cycles_unblock(state, V));
    } else {
        for (igraph_integer_t i = 0; i < num_neighbors; ++i) {
            igraph_integer_t W = VECTOR(*neighbors)[i];
            igraph_integer_t pos;
            if (!igraph_vector_int_search(igraph_adjlist_get(&state->B, W), 0, V, &pos)) {
                IGRAPH_CHECK(igraph_vector_int_push_back(igraph_adjlist_get(&state->B, W), V));
            }
        }
    }

    IGRAPH_ASSERT(!igraph_vector_int_empty(&state->v_stack));

    // unstack v
    // printf("Unstacking %lld\n", V);
    igraph_vector_int_pop_back(&state->v_stack);
    if (!igraph_vector_int_empty(&state->e_stack)) {
        // can be empty for the starting point.
        // alternatively, V == S
        igraph_vector_int_pop_back(&state->e_stack);
    }

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
    igraph_simple_cycle_search_state_t *state, const igraph_t *graph
) {
    state->N = igraph_vcount(graph);

    IGRAPH_CHECK(igraph_vector_int_init(&state->v_stack, 0));
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->v_stack, 8));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &state->v_stack);
    IGRAPH_CHECK(igraph_vector_int_init(&state->e_stack, 0));
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->e_stack, 8));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &state->e_stack);
    IGRAPH_CHECK(igraph_vector_bool_init(&state->blocked, state->N));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &state->blocked);
    IGRAPH_CHECK(igraph_inclist_init(graph, &state->IK, IGRAPH_OUT, IGRAPH_LOOPS_ONCE)); // TODO: understand what we actually want to include
    IGRAPH_FINALLY(igraph_inclist_destroy, &state->IK);
    IGRAPH_CHECK(igraph_adjlist_init_from_inclist(graph, &state->AK, &state->IK));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &state->AK);
    igraph_adjlist_sort(&state->AK);
    state->directed = igraph_is_directed(graph);
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
    igraph_vector_int_destroy(&state->v_stack);
    igraph_vector_int_destroy(&state->e_stack);
    igraph_vector_bool_destroy(&state->blocked);
    igraph_adjlist_destroy(&state->AK);
    igraph_inclist_destroy(&state->IK);
    igraph_adjlist_destroy(&state->B);
}

/**
 * \function igraph_simple_cycles_search_from_one_vertex
 * \brief Search simple cycles starting from one vertex
 *
 * \experimental
 *
 * \param state The state structure to search on
 * \param s The vertex index to start search with
 * \param results The vertices of each cycle will be stored here
 * \param search_mode What to do with undirected graphs. Ignored for directed graphs.
 *                    Possible values:
 *        \clist
 *        \cli IGRAPH_UNDIRECTED_CYCLE_SEARCH_BOTH
 *          For undirected graphs, each loop will be returned twice,
 *          each once in one and once in the other direction.
 *        \cli IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE
 *          For undirected graphs, the double loops will be filtered out.
 *          This has considerable performance implications, currently,
 *          and additionally leads to the returned loops' vertices being sorted.
 *        \endclist
 * \return Error code.
 *
 * @see https://en.wikipedia.org/wiki/Johnson%27s_algorithm
 * @see https://stackoverflow.com/a/35922906/3909202
 * @see https://epubs.siam.org/doi/epdf/10.1137/0204007
 */
igraph_error_t igraph_simple_cycles_search_from_one_vertex(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t s,
    igraph_vector_int_list_t *v_results,
    igraph_vector_int_list_t *e_results,
    igraph_simple_cycle_search_mode_t search_mode
) {
    // L3:
    for (igraph_integer_t i = s; i < state->N; ++i) {
        VECTOR(state->blocked)[i] = false;
        igraph_vector_int_clear(igraph_adjlist_get(&state->B, i));
    }

    bool found = false;
    IGRAPH_CHECK(igraph_i_simple_cycles_circuit(state, s, -1, s, v_results, e_results, &found, search_mode));

    for (igraph_integer_t i = 0; i < state->N; ++i) {
        // we want to remove the vertex with value s, not at position s
        igraph_integer_t pos;
        if (igraph_vector_int_search(igraph_adjlist_get(&state->AK, i), 0, s, &pos)) {
            igraph_vector_int_remove(igraph_adjlist_get(&state->AK, i), pos);
            igraph_vector_int_remove(igraph_inclist_get(&state->IK, i), pos);
        }
    }
    igraph_vector_int_clear(igraph_adjlist_get(&state->AK, s));
    igraph_vector_int_clear(igraph_inclist_get(&state->IK, s));

    // TODO: currently, only the nodes are returned
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
 * \param v_results Initialized integer vector list; The vertices of each cycle will be stored here
 * \param e_results Initialized integer vector list; The edges of each cycle will be stored here
 * \param search_mode How search should handle undirected graphs.
 *    See \ref igraph_simple_cycles_search_from_one_vertex
 * \return Error code.
 */
igraph_error_t igraph_simple_cycles_search_all(
    const igraph_t *graph,
    igraph_vector_int_list_t *v_result,
    igraph_vector_int_list_t *e_result,
    igraph_simple_cycle_search_mode_t search_mode
) {
    igraph_simple_cycle_search_state_t state;
    igraph_integer_t i;

    IGRAPH_CHECK(igraph_simple_cycle_search_state_init(&state, graph));
    IGRAPH_FINALLY(igraph_simple_cycle_search_state_destroy, &state);

    // TODO: depending on the graph, it is rather unreasonable to search cycles from each and every node
    for (i = 0; i < state.N; i++) {
        if (!igraph_vector_int_empty(igraph_adjlist_get(&state.AK, i))) {
            IGRAPH_CHECK(igraph_simple_cycles_search_from_one_vertex(&state, i, v_result, e_result, search_mode));
        }
    }

    igraph_simple_cycle_search_state_destroy(&state);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
