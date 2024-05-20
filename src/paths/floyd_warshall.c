/*
   IGraph library.
   Copyright (C) 2022-2023  The igraph development team <igraph@igraph.org>

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
#include "igraph_interface.h"
#include "igraph_stack.h"

#include "core/interruption.h"
#include "internal/utils.h"

static igraph_error_t distances_floyd_warshall_original(igraph_matrix_t *res) {

    igraph_integer_t no_of_nodes = igraph_matrix_nrow(res);

    for (igraph_integer_t k = 0; k < no_of_nodes; k++) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* Iteration order matters for performance!
         * First j, then i, because matrices are stored as column-major. */
        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            igraph_real_t dkj = MATRIX(*res, k, j);
            if (dkj == IGRAPH_INFINITY) {
                continue;
            }

            for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
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


static igraph_error_t distances_floyd_warshall_tree(igraph_matrix_t *res) {

    /* This is the "Tree" algorithm of Brodnik et al.
     * A difference from the paper is that instead of using the OUT_k tree of shortest
     * paths _starting_ in k, we use the IN_k tree of shortest paths _ending_ in k.
     * This makes it easier to iterate through matrices in column-major order,
     * i.e. storage order, thus increasing performance. */

    igraph_integer_t no_of_nodes = igraph_matrix_nrow(res);

    /* successors[v][u] is the second vertex on the shortest path from v to u,
       i.e. the parent of v in the IN_u tree. */
    igraph_matrix_int_t successors;
    IGRAPH_MATRIX_INT_INIT_FINALLY(&successors, no_of_nodes, no_of_nodes);

    /* children[children_start[u] + i] is the i-th child of u in a tree of shortest paths
       rooted at k, and ending in k, in the main loop below (IN_k). There are no_of_nodes-1
       child vertices in total, as the root vertex is excluded. This is essentially a contiguously
       stored adjacency list representation of IN_k. */
    igraph_vector_int_t children;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&children, no_of_nodes-1);

    /* children_start[u] indicates where the children of u are stored in children[].
       These are effectively the cumulative sums of no_of_children[], with the first
       element being 0. The last element, children_start[no_of_nodes], is equal to the
       total number of children in the tree, i.e. no_of_nodes-1. */
    igraph_vector_int_t children_start;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&children_start, no_of_nodes+1);

    /* no_of_children[u] is the number of children that u has in IN_k in the main loop below. */
    igraph_vector_int_t no_of_children;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&no_of_children, no_of_nodes);

    /* dfs_traversal and dfs_skip arrays for running time optimization,
       see "Practical improvement" in Section 3.1 of the paper */
    igraph_vector_int_t dfs_traversal;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dfs_traversal, no_of_nodes);
    igraph_vector_int_t dfs_skip;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dfs_skip, no_of_nodes);

    igraph_stack_int_t stack;
    IGRAPH_STACK_INT_INIT_FINALLY(&stack, no_of_nodes);

    for (igraph_integer_t u = 0; u < no_of_nodes; u++) {
        for (igraph_integer_t v = 0; v < no_of_nodes; v++) {
            MATRIX(successors, v, u) = u;
        }
    }

    for (igraph_integer_t k = 0; k < no_of_nodes; k++) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* Count the children of each node in the shortest path tree, assuming that at
           this point all elements of no_of_children[] are zeros. */
        for (igraph_integer_t v = 0; v < no_of_nodes; v++) {
            if (v == k) continue;
            igraph_integer_t parent = MATRIX(successors, v, k);
            VECTOR(no_of_children)[parent]++;
        }

        /* Note: we do not use igraph_vector_int_cumsum() here as that function produces
           an output vector of the same length as the input vector. Here we need an output
           one longer, with a 0 being prepended to what vector_cumsum() would produce. */
        igraph_integer_t cumsum = 0;
        for (igraph_integer_t v = 0; v < no_of_nodes; v++) {
            VECTOR(children_start)[v] = cumsum;
            cumsum += VECTOR(no_of_children)[v];
        }
        VECTOR(children_start)[no_of_nodes] = cumsum;

        /* Constructing the tree IN_k (as in the paper) and representing it
           as a contiguously stored adjacency list. The entries of the no_of_children
           vector as re-used as an index of where to insert child node indices.
           At the end of the calculation, all elements of no_of_children[] will be zeros,
           making this vector ready for the next iteration of the outer loop. */
        for (igraph_integer_t v = 0; v < no_of_nodes; v++) {
            if (v == k) continue;
            igraph_integer_t parent = MATRIX(successors, v, k);
            VECTOR(no_of_children)[parent]--;
            VECTOR(children)[ VECTOR(children_start)[parent] + VECTOR(no_of_children)[parent] ] = v;
        }

        /* constructing dfs-traversal and dfs-skip arrays for the IN_k tree */
        IGRAPH_CHECK(igraph_stack_int_push(&stack, k));
        igraph_integer_t counter = 0;
        while (!igraph_stack_int_empty(&stack)) {
            igraph_integer_t parent = igraph_stack_int_pop(&stack);
            if (parent >= 0) {
                VECTOR(dfs_traversal)[counter] = parent;
                counter++;
                /* a negative marker -parent - 1 that is popped right after
                   all the descendants of the parent were processed */
                IGRAPH_CHECK(igraph_stack_int_push(&stack, -parent - 1));
                for (igraph_integer_t l = VECTOR(children_start)[parent]; l < VECTOR(children_start)[parent + 1]; l++) {
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, VECTOR(children)[l]));
                }
            } else {
                VECTOR(dfs_skip)[-(parent + 1)] = counter;
            }
        }

        /* main inner loop */
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            igraph_real_t dki = MATRIX(*res, k, i);
            if (dki == IGRAPH_INFINITY || i == k) {
                continue;
            }
            igraph_integer_t counter = 1;
            while (counter < no_of_nodes) {
                igraph_integer_t j = VECTOR(dfs_traversal)[counter];
                igraph_real_t di = MATRIX(*res, j, k) + dki;
                igraph_real_t dd = MATRIX(*res, j, i);
                if (di < dd) {
                    MATRIX(*res, j, i) = di;
                    MATRIX(successors, j, i) = MATRIX(successors, j, k);
                    counter++;
                } else {
                    counter = VECTOR(dfs_skip)[j];
                }
                if (i == j && MATRIX(*res, i, i) < 0) {
                    IGRAPH_ERROR("Negative cycle found while calculating distances with Floyd-Warshall.",
                                 IGRAPH_ENEGLOOP);
                }
            }
        }
    }

    igraph_stack_int_destroy(&stack);
    igraph_vector_int_destroy(&dfs_traversal);
    igraph_vector_int_destroy(&dfs_skip);
    igraph_vector_int_destroy(&no_of_children);
    igraph_vector_int_destroy(&children_start);
    igraph_vector_int_destroy(&children);
    igraph_matrix_int_destroy(&successors);
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}

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
 * </para><para>
 * In addition to the original Floyd-Warshall algorithm, igraph contains implementations
 * of variants that offer better asymptotic complexity as well as better practical
 * running times for most instances. See the reference below for more information.
 *
 * </para><para>
 * Note that internally this function always computes the distance matrix
 * for all pairs of vertices. The \p from and \p to parameters only serve
 * to subset this matrix, but do not affect the time or memory taken by the
 * calculation.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Brodnik, A., Grgurovič, M., Požar, R.:
 * Modifications of the Floyd-Warshall algorithm with nearly quadratic expected-time,
 * Ars Mathematica Contemporanea, vol. 22, issue 1, p. #P1.01 (2021).
 * https://doi.org/10.26493/1855-3974.2467.497
 *
 * \param graph The graph object.
 * \param res An intialized matrix, the distances will be stored here.
 * \param from The source vertices.
 * \param to The target vertices.
 * \param weights The edge weights. If \c NULL, all weights are assumed to be 1.
 *   Negative weights are allowed, but the graph must not contain negative cycles.
 *   Edges with positive infinite weights are ignored.
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
 * \param method The type of the algorithm used.
 *        \clist
 *        \cli IGRAPH_FLOYD_WARSHALL_AUTOMATIC
 *          tried to select the best performing variant for the current graph;
 *          presently this option always uses the "Tree" method.
 *        \cli IGRAPH_FLOYD_WARSHALL_ORIGINAL
 *          the basic Floyd-Warshall algorithm.
 *        \cli IGRAPH_FLOYD_WARSHALL_TREE
 *          the "Tree" speedup of Brodnik et al., faster than the original algorithm
 *          in most cases.
 *        \endclist
 * \return Error code. \c IGRAPH_ENEGLOOP is returned if a negative-weight
 *   cycle is found.
 *
 * \sa \ref igraph_distances(), \ref igraph_distances_dijkstra(),
 * \ref igraph_distances_bellman_ford(), \ref igraph_distances_johnson()
 *
 * Time complexity:
 * The original variant has complexity O(|V|^3 + |E|).
 * The "Tree" variant has expected-case complexity of O(|V|^2 log^2 |V|)
 * according to Brodnik et al., while its worst-time complexity remains O(|V|^3).
 * Here |V| denotes the number of vertices and |E| is the number of edges.
 */
igraph_error_t igraph_distances_floyd_warshall(
        const igraph_t *graph, igraph_matrix_t *res,
        igraph_vs_t from, igraph_vs_t to,
        const igraph_vector_t *weights, igraph_neimode_t mode,
        const igraph_floyd_warshall_algorithm_t method) {

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
        IGRAPH_ERROR("Invalid mode for Floyd-Warshall shortest path calculation.", IGRAPH_EINVMODE);
    }

    if (weights && igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_fill(res, IGRAPH_INFINITY);

    for (igraph_integer_t v = 0; v < no_of_nodes; v++) {
        MATRIX(*res, v, v) = 0;
    }

    for (igraph_integer_t e = 0; e < no_of_edges; e++) {
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
        } else if (w == IGRAPH_INFINITY) {
            /* Ignore edges with infinite weight */
            continue;
        }

        if (out && MATRIX(*res, from, to) > w) {
            MATRIX(*res, from, to) = w;
        }
        if (in  && MATRIX(*res, to, from) > w) {
            MATRIX(*res, to, from) = w;
        }
    }

    /* If there are zero or one vertices, nothing needs to be done.
     * This is special-cased so that at later stages we can rely on no_of_nodes - 1 >= 0. */
    if (no_of_nodes <= 1) {
        return IGRAPH_SUCCESS;
    }

    switch (method) {
    case IGRAPH_FLOYD_WARSHALL_ORIGINAL:
        IGRAPH_CHECK(distances_floyd_warshall_original(res));
        break;
    case IGRAPH_FLOYD_WARSHALL_AUTOMATIC:
    case IGRAPH_FLOYD_WARSHALL_TREE:
        IGRAPH_CHECK(distances_floyd_warshall_tree(res));
        break;
    default:
        IGRAPH_ERROR("Invalid method.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_matrix_subset_vertices(res, graph, from, to));

    return IGRAPH_SUCCESS;
}
