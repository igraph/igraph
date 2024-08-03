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
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

#include "core/interruption.h"

// #include <math.h>

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
    igraph_vector_bool_t blocked;

    /* Whether the graph is directed */
    igraph_bool_t directed;

    /* Hashes of the vertices of found cycles (for filtering in undirected graphs)
     */
    igraph_vector_int_list_t found_cycles_vertex_hashes;

    /* Hashes of the edges of found cycles (for filtering in undirected graphs) */
    igraph_vector_int_list_t found_cycles_edge_hashes;
};

/**
 * \experimental
 *
 * A custom function to prevent double results in the search when using
 * undirected graphs
 *
 * TODO: improve performance of this by using some sort of hash table lookup
 * rather than an expensive iteration over all found cycles.
 */

static igraph_bool_t
igraph_i_cycle_has_been_found_already(igraph_simple_cycle_search_state_t *state,
                                      igraph_vector_int_t *vertex_hash,
                                      igraph_vector_int_t *edge_hash) {

    IGRAPH_ASSERT(state != NULL && vertex_hash != NULL && edge_hash != NULL);
    IGRAPH_ASSERT(
        igraph_vector_int_list_size(&state->found_cycles_vertex_hashes) ==
        igraph_vector_int_list_size(&state->found_cycles_edge_hashes));

    for (igraph_integer_t i = 0;
         i < igraph_vector_int_list_size(&state->found_cycles_vertex_hashes);
         ++i) {

        const igraph_vector_int_t *v_hashes =
            igraph_vector_int_list_get_ptr(&state->found_cycles_vertex_hashes, i);
        const igraph_vector_int_t *e_hashes =
            igraph_vector_int_list_get_ptr(&state->found_cycles_edge_hashes, i);

        if ((igraph_vector_int_size(v_hashes) != igraph_vector_int_size(vertex_hash)) ||
            (igraph_vector_int_size(e_hashes) != igraph_vector_int_size(edge_hash))) {
            continue;
        }
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(e_hashes); ++j) {
            if (VECTOR(*e_hashes)[j] != VECTOR(*edge_hash)[j]) {
                goto next_cycle;
            }
        }
        // TODO: reconsider whether it is even possible to find a cycle with the
        // same edges but different vertices
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(v_hashes); ++j) {
            if (VECTOR(*v_hashes)[j] != VECTOR(*vertex_hash)[j]) {
                goto next_cycle;
            }
        }
        return true;
next_cycle:
        continue;
    }
    return false;
}

/**
 * \experimental
 *
 * Simplify the found edge- or vertex-list to a format that
 * is irrespective of the order of vertices or edges
 */
static void igraph_i_hash_vector_int(igraph_vector_int_t *vec,
                                     igraph_vector_int_t *hash) {
    // iterate igraph_integer_t in vec, set the corresponding bit in the
    // appropriate entry in the hash vector
    igraph_integer_t max_val = igraph_vector_int_max(vec);
    igraph_integer_t n_bits_per_entries = IGRAPH_INTEGER_SIZE;
    // round up to the nearest multiple of n_bits_per_entries
    igraph_integer_t n_entries_in_hash =
        (max_val + n_bits_per_entries) / n_bits_per_entries;
    igraph_vector_int_resize(hash, n_entries_in_hash);
    igraph_vector_int_null(hash);
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(vec); ++i) {
        igraph_integer_t rel_idx = VECTOR(*vec)[i] / n_bits_per_entries;
        VECTOR(*hash)[rel_idx] =
            VECTOR(*hash)[rel_idx] | (1 << (VECTOR(*vec)[i] % n_bits_per_entries));
    }
}

/**
 * \experimental
 *
 * The implementation of procedure UNBLOCK from Johnson's paper
 *
 */

static igraph_error_t
igraph_i_simple_cycles_unblock(igraph_simple_cycle_search_state_t *state,
                               igraph_integer_t u) {
    // TODO: introduce stack for w & neis in order to reduce the number of
    // iterations.
    igraph_vector_int_t *neis;
    igraph_integer_t w;
    igraph_stack_int_t u_stack;
    IGRAPH_STACK_INT_INIT_FINALLY(&u_stack, 0);
    IGRAPH_CHECK(igraph_stack_int_push(&u_stack, u));

    while (igraph_stack_int_size(&u_stack) > 0) {
        igraph_integer_t current_u = igraph_stack_int_top(&u_stack);
        VECTOR(state->blocked)[current_u] = false;

        neis = igraph_adjlist_get(&state->B, current_u);
        igraph_bool_t recurse_deeper = false;
        while (!igraph_vector_int_empty(neis) && !recurse_deeper) {
            w = igraph_vector_int_pop_back(neis);
            if (VECTOR(state->blocked)[w]) {
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
 * - state: local state object of the search
 * - V: vertex to start the search from
 * - E: ID of edge leading into this vertex, -1 if this is the first vertex
 * - S:
 * - vertices: vector in which the results with the vertex IDs should be stored
 * - edges: vector in which the results with the edge IDs should be stored
 * - found: output argument, set to true if a cycle was found
 */

static igraph_error_t igraph_i_simple_cycles_circuit(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t V,
    igraph_integer_t E, igraph_integer_t S, igraph_vector_int_list_t *vertices,
    igraph_vector_int_list_t *edges, igraph_bool_t *found) {
    const igraph_vector_int_t *neighbors;
    const igraph_vector_int_t *incident_edges;
    igraph_integer_t num_neighbors;

    igraph_bool_t local_found = false;

    // keep track of what we were doing (rather than recursing)
    igraph_stack_int_t neigh_iteration_progress;
    IGRAPH_STACK_INT_INIT_FINALLY(&neigh_iteration_progress, 0);

    igraph_stack_int_t v_stack;
    IGRAPH_STACK_INT_INIT_FINALLY(&v_stack, 0);

    igraph_stack_int_t e_stack;
    IGRAPH_STACK_INT_INIT_FINALLY(&e_stack, 0);

    igraph_bool_t recurse_deeper = true;
    while (recurse_deeper ||
           igraph_stack_int_size(&neigh_iteration_progress) > 0) {

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
            VECTOR(state->blocked)[V] = true;
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
            // NOTE: possibly dangerous fix for undirected graphs,
            // disabling finding any two-vertex-loops
            if (W == S) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, WE));
                if (state->directed ||
                    igraph_vector_int_size(&state->vertex_stack) > 2 ||
                    E != WE) {

                    local_found = true;
                    // output circuit composed of stack
                    // printf("Found cycle with size %" IGRAPH_PRId "\n",
                    //        igraph_vector_int_size(&state->vertex_stack));

                    // copy output: from stack to vector. No need to reverse because
                    // we were putting vertices in the stack in reverse order anyway.
                    igraph_vector_int_t v_res;
                    IGRAPH_CHECK(igraph_vector_int_init_copy(&v_res, &state->vertex_stack));
                    IGRAPH_FINALLY(igraph_vector_int_destroy, &v_res);

                    // same for edges
                    igraph_vector_int_t e_res;
                    IGRAPH_CHECK(igraph_vector_int_init_copy(&e_res, &state->edge_stack));
                    IGRAPH_FINALLY(igraph_vector_int_destroy, &e_res);

                    // undirected graphs lead to some cycles being found multiple
                    // times.
                    // this is our naïve filter for now
                    igraph_bool_t persist_result = true;
                    if (!state->directed) {
                        igraph_vector_int_t vertex_hash;
                        IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_hash, 0);
                        igraph_i_hash_vector_int(&v_res, &vertex_hash);
                        igraph_vector_int_t edge_hash;
                        IGRAPH_VECTOR_INT_INIT_FINALLY(&edge_hash, 0);
                        igraph_i_hash_vector_int(&e_res, &edge_hash);
                        // printf("Has hashes %llu, %llu\n", vertex_hash, edge_hash);
                        const igraph_bool_t duplicate_found =
                            igraph_i_cycle_has_been_found_already(state,
                                                                  &vertex_hash, &edge_hash);
                        if (duplicate_found) {
                            persist_result = false;
                            igraph_vector_int_destroy(&vertex_hash);
                            igraph_vector_int_destroy(&edge_hash);
                            IGRAPH_FINALLY_CLEAN(2);
                        } else {
                            IGRAPH_CHECK(igraph_vector_int_list_push_back(
                                &state->found_cycles_vertex_hashes,
                                &vertex_hash));
                            IGRAPH_CHECK(igraph_vector_int_list_push_back(
                                &state->found_cycles_edge_hashes,
                                &edge_hash));
                            IGRAPH_FINALLY_CLEAN(2);
                        }
                    }
                    // end filter

                    if (persist_result) {
                        /* Order is important: e_res is at the top of the finally
                         * stack so we need to deal with it first */
                        if (edges != NULL) {
                            /* e_res ownership transferred to 'edges' */
                            IGRAPH_CHECK(igraph_vector_int_list_push_back(edges, &e_res));
                        } else {
                            igraph_vector_int_destroy(&e_res);
                        }
                        IGRAPH_FINALLY_CLEAN(1);

                        /* v_res ownership transferred to 'vertices' */
                        IGRAPH_CHECK(igraph_vector_int_list_push_back(vertices, &v_res));
                        IGRAPH_FINALLY_CLEAN(1);
                    } else {
                        igraph_vector_int_destroy(&v_res);
                        igraph_vector_int_destroy(&e_res);
                        IGRAPH_FINALLY_CLEAN(2);
                    }
                }
                igraph_vector_int_pop_back(&state->edge_stack);
            } else if (!(VECTOR(state->blocked)[W])) {
                // printf("Recursing deeper from %" IGRAPH_PRId " to  %" IGRAPH_PRId
                // "\n",
                //        V, W);

                recurse_deeper = true;
                IGRAPH_CHECK(igraph_stack_int_push(&neigh_iteration_progress, i + 1));
                IGRAPH_CHECK(igraph_stack_int_push(&v_stack, V));
                IGRAPH_CHECK(igraph_stack_int_push(&e_stack, E));
                V = W;
                E = WE;
                break;
            } else {
                // printf("Vertex %" IGRAPH_PRId ", neighbour of %" IGRAPH_PRId
                //        ", is blocked\n",
                //        W, V);
            }
        }

        if (!recurse_deeper) {
            // L2
            if (local_found) {
                // IGRAPH_CHECK(igraph_i_simple_cycles_unblock_recursive(state, V));
                IGRAPH_CHECK(igraph_i_simple_cycles_unblock(state, V));
            } else {
                for (igraph_integer_t i = 0; i < num_neighbors; ++i) {
                    igraph_integer_t W = VECTOR(*neighbors)[i];
                    if (!igraph_vector_int_contains(igraph_adjlist_get(&state->B, W), V)) {
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
    }

    *found = local_found;

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
igraph_error_t
igraph_simple_cycle_search_state_init(igraph_simple_cycle_search_state_t *state,
                                      const igraph_t *graph) {
    state->N = igraph_vcount(graph);
    state->NE = igraph_ecount(graph);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&state->vertex_stack, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->vertex_stack, 8));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&state->edge_stack, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->edge_stack, 8));
    IGRAPH_CHECK(igraph_vector_bool_init(&state->blocked, state->N));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &state->blocked);
    IGRAPH_CHECK(igraph_inclist_init(
                     graph, &state->IK, IGRAPH_OUT,
                     IGRAPH_LOOPS_ONCE)); // TODO: understand what we actually want to include
    IGRAPH_FINALLY(igraph_inclist_destroy, &state->IK);
    IGRAPH_CHECK(igraph_adjlist_init_from_inclist(graph, &state->AK, &state->IK));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &state->AK);
    igraph_adjlist_sort(&state->AK);
    state->directed = igraph_is_directed(graph);
    IGRAPH_CHECK(igraph_adjlist_init_empty(&state->B, state->N));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &state->B);
    if (!state->directed) {
        igraph_vector_int_list_init(&state->found_cycles_vertex_hashes, 0);
        IGRAPH_FINALLY(igraph_vector_int_list_destroy,
                       &state->found_cycles_vertex_hashes);
        igraph_vector_int_list_init(&state->found_cycles_edge_hashes, 0);
        IGRAPH_FINALLY(igraph_vector_int_list_destroy,
                       &state->found_cycles_edge_hashes);
    }

    IGRAPH_FINALLY_CLEAN(state->directed ? 6 : 8);

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
void igraph_simple_cycle_search_state_destroy(
    igraph_simple_cycle_search_state_t *state) {
    igraph_vector_int_destroy(&state->vertex_stack);
    igraph_vector_int_destroy(&state->edge_stack);
    igraph_vector_bool_destroy(&state->blocked);
    igraph_adjlist_destroy(&state->AK);
    igraph_inclist_destroy(&state->IK);
    igraph_adjlist_destroy(&state->B);
    if (!state->directed) {
        igraph_vector_int_list_destroy(&state->found_cycles_vertex_hashes);
        igraph_vector_int_list_destroy(&state->found_cycles_edge_hashes);
    }
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
 * \return Error code.
 *
 * @see https://en.wikipedia.org/wiki/Johnson%27s_algorithm
 * @see https://stackoverflow.com/a/35922906/3909202
 * @see https://epubs.siam.org/doi/epdf/10.1137/0204007
 */
igraph_error_t igraph_simple_cycles_search_from_one_vertex(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t s,
    igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges) {
    // L3:
    for (igraph_integer_t i = s; i < state->N; ++i) {
        VECTOR(state->blocked)[i] = false;
        igraph_vector_int_clear(igraph_adjlist_get(&state->B, i));
    }

    igraph_bool_t found = false;
    IGRAPH_CHECK(igraph_i_simple_cycles_circuit(state, s, -1, s, vertices, edges, &found));

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
 * \return Error code.
 */
igraph_error_t
igraph_simple_cycles_search_all(const igraph_t *graph,
                                igraph_vector_int_list_t *v_result,
                                igraph_vector_int_list_t *e_result) {

    igraph_simple_cycle_search_state_t state;

    IGRAPH_CHECK(igraph_simple_cycle_search_state_init(&state, graph));
    IGRAPH_FINALLY(igraph_simple_cycle_search_state_destroy, &state);

    // TODO: depending on the graph, it is rather unreasonable to search cycles
    // from each and every node
    for (igraph_integer_t i = 0; i < state.N; i++) {
        if (!igraph_vector_int_empty(igraph_adjlist_get(&state.AK, i))) {
            IGRAPH_CHECK(igraph_simple_cycles_search_from_one_vertex(
                             &state, i, v_result, e_result));
        }
    }

    igraph_simple_cycle_search_state_destroy(&state);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
