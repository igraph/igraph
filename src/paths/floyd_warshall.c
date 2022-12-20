/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_paths.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

/**
 * \function igraph_distances_floyd_warshall
 * \brief Weighted all-pairs shortest path lengths with the Floyd-Warshall algorithm.
 *
 * \experimental
 *
 * The Floyd-Warshall algorithm computes weighted shortest path lengths between
 * all pairs of vertices at the same time. It is useful with very dense weighted graphs,
 * as its running time is primarily determined by the vertex count, and is not sensitive
 * to the graph density. In sparse graphs, other methods such as the Dijkstra or
 * Bellman-Ford algorithms will perform significantly better.
 *
 * \param graph The graph object.
 * \param res An intialized matrix, the distances will be stored here.
 * \param weights The edge weights. If \c NULL, all weights are assumed to be 1.
 *   Negative weights are allowed, but the graph must not contain negative cycles.
 * \param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \return Error code. \c IGRAPH_ENEGLOOP is returned if a negative-weight
 *   cycle is found.
 *
 * \sa \ref igraph_distances(), \ref igraph_distances_dijkstra(),
 * \ref igraph_distances_bellman_ford(), \ref igraph_distances_johnson()
 *
 * Time complexity: O(|V|^3 + |E|) where |V| is the number of vertices
 * and |E| is the number of edges.
 */
igraph_error_t igraph_distances_floyd_warshall(
        const igraph_t *graph, igraph_matrix_t *res,
        const igraph_vector_t *weights, igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t in = false, out = false;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    switch (mode) {
    case IGRAPH_ALL:
        in = out = true;
        break;
    case IGRAPH_OUT:
        out = true;
        break;
    case IGRAPH_IN:
        in = true;
        break;
    default:
        IGRAPH_ERROR("Invalid mode.", IGRAPH_EINVAL);
    }

    if (weights && igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_fill(res, IGRAPH_INFINITY);

    for (igraph_integer_t v=0; v < no_of_nodes; v++) {
        MATRIX(*res, v, v) = 0;
    }

    for (igraph_integer_t e=0; e < no_of_edges; e++) {
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);
        igraph_real_t w = weights ? VECTOR(*weights)[e] : 1;

        if (w < 0) {
            if (mode == IGRAPH_ALL) {
                IGRAPH_ERRORF("Negative edge weight (%g) found in undirected graph "
                              "while calculating distances with Floyd-Warshall.",
                              IGRAPH_ENEGLOOP, w);
            } else if (to == from) {
                IGRAPH_ERRORF("Self-loop with negative weight (%g) found "
                              "while calculating distances with Floyd-Warshall.",
                              IGRAPH_ENEGLOOP, w);
            }
        }

        if (out && MATRIX(*res, from, to) > w) MATRIX(*res, from, to) = w;
        if (in  && MATRIX(*res, to, from) > w) MATRIX(*res, to, from) = w;
    }

    for (igraph_integer_t k=0; k < no_of_nodes; k++) {
        /* Iteration order matters for performance!
         * First j, then i, because matrices are stored as column-major. */
        for (igraph_integer_t j=0; j < no_of_nodes; j++) {
            igraph_real_t dkj = MATRIX(*res, k, j);
            if (dkj == IGRAPH_INFINITY) continue;

            for (igraph_integer_t i=0; i < no_of_nodes; i++) {
                igraph_real_t di = MATRIX(*res, i, k) + dkj;
                igraph_real_t dd = MATRIX(*res, i, j);
                if (di < dd) {
                    MATRIX(*res, i, j) = di;
                }
                if (i == j && MATRIX(*res, i, i) < 0) {
                    IGRAPH_ERROR("Negative cycle found while calculating distances with Floyd-Warshall.",
                                 IGRAPH_ENEGLOOP);
                }
            }
        }
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_distances_floyd_warshall_tree_speedup(
        const igraph_t *graph, igraph_matrix_t *res,
        const igraph_vector_t *weights, igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t in = false, out = false;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    switch (mode) {
    case IGRAPH_ALL:
        in = out = true;
        break;
    case IGRAPH_OUT:
        out = true;
        break;
    case IGRAPH_IN:
        in = true;
        break;
    default:
        IGRAPH_ERROR("Invalid mode.", IGRAPH_EINVAL);
    }

    if (weights && igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_fill(res, IGRAPH_INFINITY);

    for (igraph_integer_t v=0; v < no_of_nodes; v++) {
        MATRIX(*res, v, v) = 0;
    }

    for (igraph_integer_t e=0; e < no_of_edges; e++) {
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);
        igraph_real_t w = weights ? VECTOR(*weights)[e] : 1;

        if (w < 0) {
            if (mode == IGRAPH_ALL) {
                IGRAPH_ERRORF("Negative edge weight (%g) found in undirected graph "
                              "while calculating distances with Floyd-Warshall.",
                              IGRAPH_ENEGLOOP, w);
            } else if (to == from) {
                IGRAPH_ERRORF("Self-loop with negative weight (%g) found "
                              "while calculating distances with Floyd-Warshall.",
                              IGRAPH_ENEGLOOP, w);
            }
        }

        if (out && MATRIX(*res, from, to) > w) MATRIX(*res, from, to) = w;
        if (in && MATRIX(*res, to, from) > w) MATRIX(*res, to, from) = w;
    }

    /* predecessors[u][v] is the second but last vertex on the shortest path from u to v */
    igraph_matrix_int_t predecessors;
    /* children[u][i] is the i-th children of u in a tree of shortest distances
       rooted at k in the main loop below */
    igraph_matrix_int_t children;
    /* no_of_children[u] is the number of children of u
       rooted at k in the main loop below */
    igraph_vector_int_t no_of_children;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&no_of_children, no_of_nodes);
    igraph_stack_int_t stack;

    IGRAPH_STACK_INT_INIT_FINALLY(&stack, no_of_nodes);

    IGRAPH_CHECK(igraph_matrix_int_init(&predecessors, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_int_destroy, &predecessors);

    IGRAPH_CHECK(igraph_matrix_int_init(&children, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_int_destroy, &children);

    for (igraph_integer_t v=0; v < no_of_nodes; v++) {
        for (igraph_integer_t u=0; u < no_of_nodes; u++) {
            MATRIX(predecessors, u, v) = u;
        }
    }

    for (igraph_integer_t k=0; k < no_of_nodes; k++) {
        /* resetting no_of_children vector */
        for (igraph_integer_t parent=0; parent < no_of_nodes; parent++)
            VECTOR(no_of_children)[parent] = 0;
        /* constructing the tree out_k (as in the paper) but
           representing it as the children matrix */
        for (igraph_integer_t v=0; v < no_of_nodes; v++) {  
            if (k == v) continue;
            igraph_integer_t parent = MATRIX(predecessors, k, v);
            MATRIX(children, parent, VECTOR(no_of_children)[parent]) = v;
            VECTOR(no_of_children)[parent]++;
        }
        for (igraph_integer_t i=0; i < no_of_nodes; i++) {
            igraph_real_t dki = MATRIX(*res, k, i);
            if (dki == IGRAPH_INFINITY || i == k) continue;
            IGRAPH_CHECK(igraph_stack_int_push(&stack, k));
            while (!igraph_stack_int_empty(&stack)) {
                igraph_integer_t parent = igraph_stack_int_pop(&stack);
                for (igraph_integer_t l=0; l < VECTOR(no_of_children)[parent]; l++) {
                    igraph_integer_t j = MATRIX(children, parent, l);
                    igraph_real_t di = MATRIX(*res, j, k) + dki;
                    igraph_real_t dd = MATRIX(*res, j, i);
                    if (di < dd) {
                        MATRIX(*res, j, i) = di;
                        MATRIX(predecessors, i, j) = MATRIX(predecessors, k, j);
                        IGRAPH_CHECK(igraph_stack_int_push(&stack, j));
                    }
                    if (i == j && MATRIX(*res, i, i) < 0) {
                        IGRAPH_ERROR("Negative cycle found while calculating distances with Floyd-Warshall.",
                                    IGRAPH_ENEGLOOP);
                    }
                }
            }
        }
    }

    igraph_matrix_int_destroy(&predecessors);
    igraph_matrix_int_destroy(&children);
    igraph_stack_int_destroy(&stack);
    igraph_vector_int_destroy(&no_of_children);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}
