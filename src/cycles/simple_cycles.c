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
typedef struct igraph_simple_cycle_search_state_t
{
    /* Number of vertices in the graph */
    igraph_integer_t N;

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

    /* Hashes of the vertices of found cycles (for filtering in undirected graphs) */
    igraph_vector_int_t found_cycles_vertex_hashes;

    /* Hashes of the edges of found cycles (for filtering in undirected graphs) */
    igraph_vector_int_t found_cycles_edge_hashes;
} igraph_simple_cycle_search_state_t;

/**
 * \experimental
 *
 * A custom function to prevent double results in the search when using undirected graphs
 */

static igraph_bool_t igraph_i_cycle_has_been_found_already(
    igraph_simple_cycle_search_state_t *state,
    igraph_vector_int_list_t *v_results, igraph_vector_int_list_t *e_results,
    igraph_vector_int_t *v_res, igraph_vector_int_t *e_res,
    unsigned long long vertex_hash, unsigned long long edge_hash)
{
    if (igraph_vector_int_size(v_res) < 2)
    {
        return false;
    }

    // TODO: improve performance of this by using some sort of hash table lookup
    igraph_integer_t pos;
    if (!igraph_vector_int_search(&state->found_cycles_edge_hashes, 0, edge_hash, &pos))
    {
        return false;
    }
    // In principle, it should not be possible to have the same edges but different vertices
    // Yet this is just another safeguard against a possible (edge) hash collision
    while (VECTOR(state->found_cycles_vertex_hashes)[pos] != vertex_hash && igraph_vector_int_size(&VECTOR(*v_results)[pos]) != igraph_vector_int_size(v_res))
    {
        if (!igraph_vector_int_search(&state->found_cycles_edge_hashes, pos + 1, edge_hash, &pos))
        {
            return false;
        }
    }

    // if all of edge_hash, vertex_hash and cycle size match, we assume we got a duplicate
    return true;
}

/**
 * \experimental
 *
 * A hashing algorithm
 */
static unsigned long long igraph_i_hash_vector_int(igraph_vector_int_t *vec)
{

    unsigned long long hash_product = 1.;
    unsigned long long hash_sum = 0;
    unsigned int hash_xor = 0;

    for (igraph_integer_t i = 0; i < igraph_vector_int_size(vec); ++i)
    {
        hash_product *= VECTOR(*vec)[i];
        hash_sum += VECTOR(*vec)[i];
        hash_xor ^= VECTOR(*vec)[i];
    }

    return hash_product + hash_sum + ((unsigned long long)hash_xor << 32);
}

/**
 * \experimental
 *
 * The implementation of procedure UNBLOCK from Johnson's paper
 *
 */

static igraph_error_t igraph_i_simple_cycles_unblock(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t u)
{
    // TODO: introduce stack for w & neis in order to reduce the number of iterations.
    igraph_vector_int_t *neis;
    igraph_integer_t w;
    igraph_stack_int_t u_stack;
    igraph_stack_int_init(&u_stack, 0);
    igraph_stack_int_push(&u_stack, u);

    while (igraph_stack_int_size(&u_stack) > 0)
    {
        igraph_integer_t current_u = igraph_stack_int_top(&u_stack);
        VECTOR(state->blocked)
        [current_u] = false;

        neis = igraph_adjlist_get(&state->B, current_u);
        bool recurse_deeper = false;
        while (!igraph_vector_int_empty(neis) && !recurse_deeper)
        {
            w = igraph_vector_int_pop_back(neis);
            if (VECTOR(state->blocked)[w])
            {
                igraph_stack_int_push(&u_stack, w);
            }
        }

        if (!recurse_deeper)
        {
            igraph_stack_int_pop(&u_stack);
        }
    }

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

static igraph_error_t igraph_i_simple_cycles_circuit_recursive(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t V,
    igraph_integer_t E, igraph_integer_t S,
    igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
    igraph_bool_t *found)
{
    igraph_bool_t local_found = false;
    igraph_vector_int_t *neighbors;
    igraph_vector_int_t *incident_edges;
    igraph_integer_t num_neighbors;

    // stack v & e
    IGRAPH_CHECK(igraph_vector_int_push_back(&state->vertex_stack, V));
    if (E >= 0)
    {
        IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, E));
    }
    // printf("Pushing %lld to stack, stack size is %lld, result size is %lld\n", V, igraph_vector_int_size(&state->vertex_stack), igraph_vector_int_list_size(results));
    VECTOR(state->blocked)
    [V] = true;

    // L1
    neighbors = igraph_adjlist_get(&state->AK, V);
    incident_edges = igraph_inclist_get(&state->IK, V);
    num_neighbors = igraph_vector_int_size(neighbors);
    IGRAPH_ASSERT(igraph_vector_int_size(incident_edges) == num_neighbors);
    for (igraph_integer_t i = 0; i < num_neighbors; ++i)
    {
        igraph_integer_t W = VECTOR(*neighbors)[i];
        igraph_integer_t WE = VECTOR(*incident_edges)[i];
        // NOTE: possibly dangerous fix for undirected graphs,
        // disabling finding any two-vertex-loops
        if (W == S)
        {
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, WE));
            if (state->directed || igraph_vector_int_size(&state->vertex_stack) > 2)
            {
                local_found = true;
                // output circuit composed of stack
                // printf("Found cycle with size %" IGRAPH_PRId "\n", igraph_vector_int_size(&state->vertex_stack));

                // copy output: from stack to vector. No need to reverse because
                // we were putting vertices in the stack in reverse order anyway.
                igraph_vector_int_t v_res;
                IGRAPH_CHECK(igraph_vector_int_init_copy(&v_res, &state->vertex_stack));
                IGRAPH_FINALLY(igraph_vector_int_destroy, &v_res);
                unsigned long long vertex_hash = igraph_i_hash_vector_int(&v_res);

                // same for edges
                igraph_vector_int_t e_res;
                IGRAPH_CHECK(igraph_vector_int_init_copy(&e_res, &state->edge_stack));
                IGRAPH_FINALLY(igraph_vector_int_destroy, &e_res);
                unsigned long long edge_hash = igraph_i_hash_vector_int(&e_res);

                // undirected graphs lead to some cycles being found multiple
                // times.
                // this is our naïve filter for now
                igraph_bool_t persist_result = true;
                if (!state->directed)
                {
                    // print_vector_int(&v_res);
                    // printf("Has hashes %llu, %llu\n", vertex_hash, edge_hash);
                    igraph_bool_t duplicate_found = igraph_i_cycle_has_been_found_already(state, vertices, edges, &v_res, &e_res, vertex_hash, edge_hash);
                    if (duplicate_found)
                    {
                        persist_result = false;
                    }
                    else
                    {
                        igraph_vector_int_push_back(&state->found_cycles_vertex_hashes, vertex_hash);
                        igraph_vector_int_push_back(&state->found_cycles_edge_hashes, edge_hash);
                    }
                }
                // end filter

                if (persist_result)
                {
                    /* Order is important: e_res is at the top of the finally
                     * stack so we need to deal with it first */
                    if (edges != NULL)
                    {
                        /* e_res ownership transferred to 'edges' */
                        IGRAPH_CHECK(igraph_vector_int_list_push_back(edges, &e_res));
                    }
                    else
                    {
                        igraph_vector_int_destroy(&e_res);
                    }
                    IGRAPH_FINALLY_CLEAN(1);

                    /* v_res ownership transferred to 'vertices' */
                    IGRAPH_CHECK(igraph_vector_int_list_push_back(vertices, &v_res));
                    IGRAPH_FINALLY_CLEAN(1);
                }
                else
                {
                    igraph_vector_int_destroy(&v_res);
                    igraph_vector_int_destroy(&e_res);
                    IGRAPH_FINALLY_CLEAN(2);
                }
            }
            igraph_vector_int_pop_back(&state->edge_stack);
        }
        else if (!(VECTOR(state->blocked)[W]))
        {
            IGRAPH_CHECK(igraph_i_simple_cycles_circuit_recursive(state, W, WE, S, vertices, edges, &local_found));
        }
    }
    *found = local_found;

    // L2
    if (local_found)
    {
        IGRAPH_CHECK(igraph_i_simple_cycles_unblock(state, V));
    }
    else
    {
        for (igraph_integer_t i = 0; i < num_neighbors; ++i)
        {
            igraph_integer_t W = VECTOR(*neighbors)[i];
            igraph_integer_t pos;
            if (!igraph_vector_int_search(igraph_adjlist_get(&state->B, W), 0, V, &pos))
            {
                IGRAPH_CHECK(igraph_vector_int_push_back(igraph_adjlist_get(&state->B, W), V));
            }
        }
    }

    IGRAPH_ASSERT(!igraph_vector_int_empty(&state->vertex_stack));

    // unstack v
    // printf("Unstacking %lld\n", V);
    igraph_vector_int_pop_back(&state->vertex_stack);
    if (!igraph_vector_int_empty(&state->edge_stack))
    {
        // can be empty for the starting point.
        // alternatively, V == S
        igraph_vector_int_pop_back(&state->edge_stack);
    }

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
    igraph_simple_cycle_search_state_t *state,
    igraph_integer_t V,
    igraph_integer_t E,
    igraph_integer_t S,
    igraph_vector_int_list_t *vertices,
    igraph_vector_int_list_t *edges,
    igraph_bool_t *found)
{
    igraph_bool_t local_found = false;
    igraph_vector_int_t *neighbors;
    igraph_vector_int_t *incident_edges;
    igraph_integer_t num_neighbors;

    igraph_stack_int_t neigh_iteration_progress;
    igraph_stack_int_init(&neigh_iteration_progress, 10);

    // stack v & e
    bool recurse_deeper = true;
    while (recurse_deeper && igraph_stack_int_size(&neigh_iteration_progress))
    {
        igraph_integer_t i0 = 0;
        if (recurse_deeper)
        {
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->vertex_stack, V));
            if (E >= 0)
            {
                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, E));
            }
            // printf("Pushing %lld to stack, stack size is %lld, result size is %lld\n", V, igraph_vector_int_size(&state->vertex_stack), igraph_vector_int_list_size(results));
            VECTOR(state->blocked)
            [V] = true;
            recurse_deeper = false;
        }
        else
        {
            i0 = igraph_stack_int_pop(&neigh_iteration_progress);
        }

        // L1
        neighbors = igraph_adjlist_get(&state->AK, V);
        incident_edges = igraph_inclist_get(&state->IK, V);
        num_neighbors = igraph_vector_int_size(neighbors);
        IGRAPH_ASSERT(igraph_vector_int_size(incident_edges) == num_neighbors);
        for (igraph_integer_t i = i0; i < num_neighbors; ++i)
        {
            igraph_integer_t W = VECTOR(*neighbors)[i];
            igraph_integer_t WE = VECTOR(*incident_edges)[i];
            // NOTE: possibly dangerous fix for undirected graphs,
            // disabling finding any two-vertex-loops
            if (W == S)
            {
                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edge_stack, WE));
                if (state->directed || igraph_vector_int_size(&state->vertex_stack) > 2)
                {
                    local_found = true;
                    // output circuit composed of stack
                    // printf("Found cycle with size %" IGRAPH_PRId "\n", igraph_vector_int_size(&state->vertex_stack));

                    // copy output: from stack to vector. No need to reverse because
                    // we were putting vertices in the stack in reverse order anyway.
                    igraph_vector_int_t v_res;
                    IGRAPH_CHECK(igraph_vector_int_init_copy(&v_res, &state->vertex_stack));
                    IGRAPH_FINALLY(igraph_vector_int_destroy, &v_res);
                    unsigned long long vertex_hash = igraph_i_hash_vector_int(&v_res);

                    // same for edges
                    igraph_vector_int_t e_res;
                    IGRAPH_CHECK(igraph_vector_int_init_copy(&e_res, &state->edge_stack));
                    IGRAPH_FINALLY(igraph_vector_int_destroy, &e_res);
                    unsigned long long edge_hash = igraph_i_hash_vector_int(&e_res);

                    // undirected graphs lead to some cycles being found multiple
                    // times.
                    // this is our naïve filter for now
                    igraph_bool_t persist_result = true;
                    if (!state->directed)
                    {
                        // print_vector_int(&v_res);
                        // printf("Has hashes %llu, %llu\n", vertex_hash, edge_hash);
                        igraph_bool_t duplicate_found = igraph_i_cycle_has_been_found_already(state, vertices, edges, &v_res, &e_res, vertex_hash, edge_hash);
                        if (duplicate_found)
                        {
                            persist_result = false;
                        }
                        else
                        {
                            igraph_vector_int_push_back(&state->found_cycles_vertex_hashes, vertex_hash);
                            igraph_vector_int_push_back(&state->found_cycles_edge_hashes, edge_hash);
                        }
                    }
                    // end filter

                    if (persist_result)
                    {
                        /* Order is important: e_res is at the top of the finally
                         * stack so we need to deal with it first */
                        if (edges != NULL)
                        {
                            /* e_res ownership transferred to 'edges' */
                            IGRAPH_CHECK(igraph_vector_int_list_push_back(edges, &e_res));
                        }
                        else
                        {
                            igraph_vector_int_destroy(&e_res);
                        }
                        IGRAPH_FINALLY_CLEAN(1);

                        /* v_res ownership transferred to 'vertices' */
                        IGRAPH_CHECK(igraph_vector_int_list_push_back(vertices, &v_res));
                        IGRAPH_FINALLY_CLEAN(1);
                    }
                    else
                    {
                        igraph_vector_int_destroy(&v_res);
                        igraph_vector_int_destroy(&e_res);
                        IGRAPH_FINALLY_CLEAN(2);
                    }
                }
                igraph_vector_int_pop_back(&state->edge_stack);
            }
            else if (!(VECTOR(state->blocked)[W]))
            {
                recurse_deeper = true;
                igraph_stack_int_push(&neigh_iteration_progress, i);
                V = W;
                E = WE;
                break;
            }
        }

        if (!recurse_deeper)
        {

            // L2
            if (local_found)
            {
                IGRAPH_CHECK(igraph_i_simple_cycles_unblock(state, V));
            }
            else
            {
                for (igraph_integer_t i = 0; i < num_neighbors; ++i)
                {
                    igraph_integer_t W = VECTOR(*neighbors)[i];
                    igraph_integer_t pos;
                    if (!igraph_vector_int_search(igraph_adjlist_get(&state->B, W), 0, V, &pos))
                    {
                        IGRAPH_CHECK(igraph_vector_int_push_back(igraph_adjlist_get(&state->B, W), V));
                    }
                }
            }

            IGRAPH_ASSERT(!igraph_vector_int_empty(&state->vertex_stack));

            // unstack v
            // printf("Unstacking %lld\n", V);
            V = igraph_vector_int_pop_back(&state->vertex_stack);
            if (!igraph_vector_int_empty(&state->edge_stack))
            {
                // can be empty for the starting point.
                // alternatively, V == S
                E = igraph_vector_int_pop_back(&state->edge_stack);
            }
        }
    }
    *found = local_found;

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
    igraph_simple_cycle_search_state_t *state, const igraph_t *graph)
{
    state->N = igraph_vcount(graph);

    IGRAPH_CHECK(igraph_vector_int_init(&state->vertex_stack, 0));
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->vertex_stack, 8));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &state->vertex_stack);
    IGRAPH_CHECK(igraph_vector_int_init(&state->edge_stack, 0));
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->edge_stack, 8));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &state->edge_stack);
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
    if (!state->directed)
    {
        IGRAPH_CHECK(igraph_vector_int_init(&state->found_cycles_vertex_hashes, 0));
        IGRAPH_CHECK(igraph_vector_int_reserve(&state->found_cycles_vertex_hashes, 8));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &state->found_cycles_vertex_hashes);
        IGRAPH_CHECK(igraph_vector_int_init(&state->found_cycles_edge_hashes, 0));
        IGRAPH_CHECK(igraph_vector_int_reserve(&state->found_cycles_edge_hashes, 8));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &state->found_cycles_edge_hashes);
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
void igraph_simple_cycle_search_state_destroy(igraph_simple_cycle_search_state_t *state)
{
    igraph_vector_int_destroy(&state->vertex_stack);
    igraph_vector_int_destroy(&state->edge_stack);
    igraph_vector_bool_destroy(&state->blocked);
    igraph_adjlist_destroy(&state->AK);
    igraph_inclist_destroy(&state->IK);
    igraph_adjlist_destroy(&state->B);
    if (!state->directed)
    {
        igraph_vector_int_destroy(&state->found_cycles_vertex_hashes);
        igraph_vector_int_destroy(&state->found_cycles_edge_hashes);
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
    igraph_vector_int_list_t *vertices,
    igraph_vector_int_list_t *edges)
{
    // L3:
    for (igraph_integer_t i = s; i < state->N; ++i)
    {
        VECTOR(state->blocked)
        [i] = false;
        igraph_vector_int_clear(igraph_adjlist_get(&state->B, i));
    }

    igraph_bool_t found = false;
    IGRAPH_CHECK(igraph_i_simple_cycles_circuit(state, s, -1, s, vertices, edges, &found));

    for (igraph_integer_t i = 0; i < state->N; ++i)
    {
        // we want to remove the vertex with value s, not at position s
        igraph_integer_t pos;
        if (igraph_vector_int_search(igraph_adjlist_get(&state->AK, i), 0, s, &pos))
        {
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
igraph_error_t igraph_simple_cycles_search_all(
    const igraph_t *graph,
    igraph_vector_int_list_t *v_result,
    igraph_vector_int_list_t *e_result)
{
    igraph_simple_cycle_search_state_t state;
    igraph_integer_t i;

    IGRAPH_CHECK(igraph_simple_cycle_search_state_init(&state, graph));
    IGRAPH_FINALLY(igraph_simple_cycle_search_state_destroy, &state);

    // TODO: depending on the graph, it is rather unreasonable to search cycles from each and every node
    for (i = 0; i < state.N; i++)
    {
        if (!igraph_vector_int_empty(igraph_adjlist_get(&state.AK, i)))
        {
            IGRAPH_CHECK(igraph_simple_cycles_search_from_one_vertex(&state, i, v_result, e_result));
        }
    }

    igraph_simple_cycle_search_state_destroy(&state);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
