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
 * \brief Connects every vertex to its neighborhood
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
 * \sa \ref igraph_lattice() uses this function to connect the
 * neighborhood of the vertices.
 *
 * Time complexity: O(|V|*d^k), |V| is the number of vertices in the
 * graph, d is the average degree and k is the \p order argument.
 */
int igraph_connect_neighborhood(igraph_t *graph, igraph_integer_t order,
                                igraph_neimode_t mode) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vector_t edges;
    long int i, j, in;
    long int *added;
    igraph_vector_t neis;

    if (order < 0) {
        IGRAPH_ERROR("Negative order, cannot connect neighborhood", IGRAPH_EINVAL);
    }

    if (order < 2) {
        IGRAPH_WARNING("Order smaller than two, graph will be unchanged");
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    added = IGRAPH_CALLOC(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot connect neighborhood", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

    for (i = 0; i < no_of_nodes; i++) {
        added[i] = i + 1;
        igraph_neighbors(graph, &neis, (igraph_integer_t) i, mode);
        in = igraph_vector_size(&neis);
        if (order > 1) {
            for (j = 0; j < in; j++) {
                long int nei = (long int) VECTOR(neis)[j];
                added[nei] = i + 1;
                igraph_dqueue_push(&q, nei);
                igraph_dqueue_push(&q, 1);
            }
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (mode != IGRAPH_ALL || i < nei) {
                            if (mode == IGRAPH_IN) {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                            } else {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                            }
                        }
                    }
                }
            } else {
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (mode != IGRAPH_ALL || i < nei) {
                            if (mode == IGRAPH_IN) {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                            } else {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                            }
                        }
                    }
                }
            }

        } /* while q not empty */
    } /* for i < no_of_nodes */

    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&q);
    igraph_free(added);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
