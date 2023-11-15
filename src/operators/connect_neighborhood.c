/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_adjlist.h"
#include "igraph_error.h"
#include "igraph_operators.h"

#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "graph/attributes.h"

/**
 * \function igraph_connect_neighborhood
 * \brief Connects each vertex to its neighborhood.
 *
 * This function adds new edges to the input graph. Each vertex is connected
 * to all vertices reachable by at most \p order steps from it
 * (unless a connection already existed).
 *
 * </para><para>
 * Note that the input graph is modified in place, no
 * new graph is created. Call \ref igraph_copy() if you want to keep
 * the original graph as well.
 *
 * </para><para>
 * For undirected graphs reachability is always
 * symmetric: if vertex A can be reached from vertex B in at
 * most \p order steps, then the opposite is also true. Only one
 * undirected (A,B) edge will be added in this case.
 *
 * \param graph The input graph. It will be modified in-place.
 * \param order Integer constant, it gives the distance within which
 *    the vertices will be connected to the source vertex.
 * \param mode Constant, it specifies how the neighborhood search is
 *    performed for directed graphs. If \c IGRAPH_OUT then vertices
 *    reachable from the source vertex will be connected, \c IGRAPH_IN
 *    is the opposite. If \c IGRAPH_ALL then the directed graph is
 *    considered as an undirected one.
 * \return Error code.
 *
 * \sa \ref igraph_graph_power() to compute the kth power of a graph;
 * \ref igraph_square_lattice() uses this function to connect the
 * neighborhood of the vertices.
 *
 * Time complexity: O(|V|*d^k), |V| is the number of vertices in the
 * graph, d is the average degree and k is the \p order argument.
 */
igraph_error_t igraph_connect_neighborhood(igraph_t *graph, igraph_integer_t order,
                                igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vector_int_t edges;
    igraph_integer_t i, j, in;
    igraph_integer_t *added;
    igraph_vector_int_t neis;

    if (order < 0) {
        IGRAPH_ERRORF("Order must not be negative, found %" IGRAPH_PRId ".",
                IGRAPH_EINVAL, order);
    }

    if (order < 2) {
        IGRAPH_WARNING("Order smaller than two, graph will be unchanged.");
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(added, "Cannot connect neighborhood.");

    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    for (i = 0; i < no_of_nodes; i++) {
        added[i] = i + 1;
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, mode));
        in = igraph_vector_int_size(&neis);
        if (order > 1) {
            for (j = 0; j < in; j++) {
                igraph_integer_t nei = VECTOR(neis)[j];
                added[nei] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, 1));
            }
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            igraph_integer_t n;
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
            n = igraph_vector_int_size(&neis);

            if (actdist < order - 1) {
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                        if (mode != IGRAPH_ALL || i < nei) {
                            if (mode == IGRAPH_IN) {
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, nei));
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                            } else {
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, nei));
                            }
                        }
                    }
                }
            } else {
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (mode != IGRAPH_ALL || i < nei) {
                            if (mode == IGRAPH_IN) {
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, nei));
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                            } else {
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, nei));
                            }
                        }
                    }
                }
            }

        } /* while q not empty */
    } /* for i < no_of_nodes */

    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&q);
    igraph_free(added);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_add_edges(graph, &edges, NULL));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_graph_power
 * \brief The kth power of a graph.
 *
 * \experimental
 *
 * The kth power of a graph G is a simple graph where vertex \c u is connected to
 * \c v by a single edge if \c v is reachable from \c u in G within at most k steps.
 * By convention, the zeroth power of a graph has no edges. The first power is
 * identical to the original graph, except that multiple edges and self-loops
 * are removed.
 *
 * </para><para>
 * Graph power is usually defined only for undirected graphs. igraph extends the concept
 * to directed graphs. To ignore edge directions in the input, set the \p directed
 * parameter to \c false. In this case, the result will be an undirected graph.
 *
 * </para><para>
 * Graph and vertex attributes are preserved, but edge attributes are discarded.
 *
 * \param graph The input graph.
 * \param res The graph power of the given \p order.
 * \param order Non-negative integer, the power to raise the graph to.
 *    In other words, vertices within a distance \p order will be connected.
 * \param directed Logical, whether to take edge directions into account.
 * \return Error code.
 *
 * \sa \ref igraph_connect_neighborhood() to connect each vertex to its
 * neighborhood, modifying a graph in-place.
 *
 * Time complexity: O(|V|*d^k), |V| is the number of vertices in the
 * graph, d is the average degree and k is the \p order argument.
 */
igraph_error_t igraph_graph_power(const igraph_t *graph, igraph_t *res,
                                  igraph_integer_t order, igraph_bool_t directed) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t edges;
    igraph_adjlist_t al;
    igraph_bool_t dir = igraph_is_directed(graph) && directed;
    igraph_neimode_t mode = dir ? IGRAPH_OUT : IGRAPH_ALL;

    if (order < 0) {
        IGRAPH_ERRORF("Order must not be negative, found %" IGRAPH_PRId ".",
                IGRAPH_EINVAL, order);
    }

    IGRAPH_CHECK(igraph_empty(res, no_of_nodes, dir));
    IGRAPH_I_ATTRIBUTE_DESTROY(res);
    IGRAPH_I_ATTRIBUTE_COPY(res, graph, /* graph */ true, /* vertex */ true, /* edge */ false);
    if (order == 0) {
        return IGRAPH_SUCCESS;
    }

    /* Initialize res with a copy of the graph, but with multi-edges and self-loops removed.
     * Also convert the graph to undirected if this is requested. */
    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    /* Reserve initial space for no_of_edges. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges);
    igraph_vector_int_clear(&edges);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *tmp = igraph_adjlist_get(&al, i);
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(tmp); j++) {
            if (dir || i < VECTOR(*tmp)[j]) {
                igraph_vector_int_push_back(&edges, i);
                igraph_vector_int_push_back(&edges, VECTOR(*tmp)[j]);
            }
        }
    }

    if (order > 1) {
        /* order > 1, so add more edges. */

        igraph_integer_t d_i, d_actnode;
        igraph_integer_t *added;
        const igraph_vector_int_t *neis;
        igraph_dqueue_int_t q;

        added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
        IGRAPH_CHECK_OOM(added, "Insufficient memory for graph power.");
        IGRAPH_FINALLY(igraph_free, added);

        IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            added[i] = i + 1;
            neis = igraph_adjlist_get(&al, i);
            d_i = igraph_vector_int_size(neis);

            for (igraph_integer_t j = 0; j < d_i; j++) {
                igraph_integer_t nei = VECTOR(*neis)[j];
                added[nei] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, 1));
            }

            while (!igraph_dqueue_int_empty(&q)) {
                igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
                igraph_integer_t actdist = igraph_dqueue_int_pop(&q);

                neis = igraph_adjlist_get(&al, actnode);
                d_actnode = igraph_vector_int_size(neis);

                if (actdist < order - 1) {
                    for (igraph_integer_t j = 0; j < d_actnode; j++) {
                        igraph_integer_t nei = VECTOR(*neis)[j];
                        if (added[nei] != i + 1) {
                            added[nei] = i + 1;
                            IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                            IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                            if (dir || i < nei) {
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, nei));
                            }
                        }
                    }
                } else {
                    for (igraph_integer_t j = 0; j < d_actnode; j++) {
                        igraph_integer_t nei = VECTOR(*neis)[j];
                        if (added[nei] != i + 1) {
                            added[nei] = i + 1;
                            if (dir || i < nei) {
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, nei));
                            }
                        }
                    }
                }

            } /* while q not empty */
        } /* for i < no_of_nodes */

        igraph_dqueue_int_destroy(&q);
        igraph_free(added);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_add_edges(res, &edges, 0));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
