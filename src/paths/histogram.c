/*
   igraph library.
   Copyright (C) 2005-2024  The igraph development team <igraph@igraph.org>

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

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_progress.h"

#include "core/interruption.h"

/**
 * \function igraph_path_length_hist
 * \brief Create a histogram of all shortest path lengths.
 *
 * This function calculates a histogram by calculating the shortest path
 * length between all pairs of vertices. In directed graphs, both directions
 * are considered, meaning that each vertex pair appears twice in the histogram.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored here. The
 *     first (i.e. index 0) element contains the number of shortest paths of
 *     length 1, the second of length 2, etc. The supplied vector is resized
 *     as needed.
 * \param unconnected Pointer to a real number, the number of vertex
 *     pairs for which the second vertex is not reachable from the
 *     first is stored here.
 * \param directed Whether to consider directed paths in a directed
 *     graph. This argument is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V||E|), the number of vertices times the number
 * of edges.
 *
 * \sa \ref igraph_average_path_length() and \ref igraph_distances()
 */

igraph_error_t igraph_path_length_hist(const igraph_t *graph, igraph_vector_t *res,
                            igraph_real_t *unconnected, igraph_bool_t directed) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t already_added;
    igraph_int_t nodes_reached;
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_adjlist_t allneis;
    igraph_real_t unconn = 0;
    igraph_int_t ressize;

    if (! igraph_is_directed(graph)) {
        directed = false;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&already_added, no_of_nodes);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis,
                                     directed ? IGRAPH_OUT : IGRAPH_ALL,
                                     IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    igraph_vector_clear(res);
    ressize = 0;

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        nodes_reached = 1;      /* itself */
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        VECTOR(already_added)[i] = i + 1;

        IGRAPH_PROGRESS("Path length histogram: ", 100.0 * i / no_of_nodes, NULL);

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_int_t actnode = igraph_dqueue_int_pop(&q);
            igraph_int_t actdist = igraph_dqueue_int_pop(&q);

            neis = igraph_adjlist_get(&allneis, actnode);
            const igraph_int_t n = igraph_vector_int_size(neis);
            for (igraph_int_t j = 0; j < n; j++) {
                igraph_int_t neighbor = VECTOR(*neis)[j];
                if (VECTOR(already_added)[neighbor] == i + 1) {
                    continue;
                }
                VECTOR(already_added)[neighbor] = i + 1;
                nodes_reached++;
                if (actdist + 1 > ressize) {
                    IGRAPH_CHECK(igraph_vector_resize(res, actdist + 1));
                    for (; ressize < actdist + 1; ressize++) {
                        VECTOR(*res)[ressize] = 0;
                    }
                }
                VECTOR(*res)[actdist] += 1;

                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
            }
        } /* while !igraph_dqueue_int_empty */

        unconn += (no_of_nodes - nodes_reached);

    } /* for i<no_of_nodes */

    IGRAPH_PROGRESS("Path length histogram: ", 100.0, NULL);

    /* count every pair only once for an undirected graph */
    if (!directed) {
        for (igraph_int_t i = 0; i < ressize; i++) {
            VECTOR(*res)[i] /= 2;
        }
        unconn /= 2;
    }

    igraph_vector_int_destroy(&already_added);
    igraph_dqueue_int_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(3);

    if (unconnected) {
        *unconnected = unconn;
    }

    return IGRAPH_SUCCESS;
}
