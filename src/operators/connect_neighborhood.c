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

#include "igraph_operators.h"

#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

/**
 * \function igraph_connect_neighborhood
 * \brief Graph power: connect each vertex to its neighborhood.
 *
 * This function adds new edges to the input graph. Each vertex is connected
 * to all vertices reachable by at most \p order steps from it
 * (unless a connection already existed).  In other words, the \p order power of
 * the graph is computed.
 *
 * </para><para> Note that the input graph is modified in place, no
 * new graph is created. Call \ref igraph_copy() if you want to keep
 * the original graph as well.
 *
 * </para><para> For undirected graphs reachability is always
 * symmetric: if vertex A can be reached from vertex B in at
 * most \p order steps, then the opposite is also true. Only one
 * undirected (A,B) edge will be added in this case.
 * \param graph The input graph, this is the output graph as well.
 * \param order Integer constant, it gives the distance within which
 *    the vertices will be connected to the source vertex.
 * \param mode Constant, it specifies how the neighborhood search is
 *    performed for directed graphs. If \c IGRAPH_OUT then vertices
 *    reachable from the source vertex will be connected, \c IGRAPH_IN
 *    is the opposite. If \c IGRAPH_ALL then the directed graph is
 *    considered as an undirected one.
 * \return Error code.
 *
 * \sa \ref igraph_square_lattice() uses this function to connect the
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
        IGRAPH_ERRORF("Order can not be negative, found %" IGRAPH_PRId ".",
                IGRAPH_EINVAL, order);
    }

    if (order < 2) {
        IGRAPH_WARNING("Order smaller than two, graph will be unchanged");
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (added == 0) {
        IGRAPH_ERROR("Cannot connect neighborhood", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
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

    IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
