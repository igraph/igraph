/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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

#include "igraph_paths.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_progress.h"

#include "core/interruption.h"

/**
 * \function igraph_path_length_hist
 * Create a histogram of all shortest path lengths.
 *
 * This function calculates a histogram, by calculating the
 * shortest path length between each pair of vertices. For directed
 * graphs both directions might be considered and then every pair of vertices
 * appears twice in the histogram.
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *     here. The first (i.e. zeroth) element contains the number of
 *     shortest paths of length 1, etc. The supplied vector is resized
 *     as needed.
 * \param unconnected Pointer to a real number, the number of
 *     pairs for which the second vertex is not reachable from the
 *     first is stored here.
 * \param directed Whether to consider directed paths in a directed
 *     graph (if not zero). This argument is ignored for undirected
 *     graphs.
 * \return Error code.
 *
 * Time complexity: O(|V||E|), the number of vertices times the number
 * of edges.
 *
 * \sa \ref igraph_average_path_length() and \ref igraph_distances()
 */

igraph_error_t igraph_path_length_hist(const igraph_t *graph, igraph_vector_t *res,
                            igraph_real_t *unconnected, igraph_bool_t directed) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i, j, n;
    igraph_vector_int_t already_added;
    igraph_integer_t nodes_reached;

    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_neimode_t dirmode;
    igraph_adjlist_t allneis;
    igraph_real_t unconn = 0;
    igraph_integer_t ressize;

    if (directed) {
        dirmode = IGRAPH_OUT;
    } else {
        dirmode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&already_added, no_of_nodes);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, dirmode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    igraph_vector_clear(res);
    ressize = 0;

    for (i = 0; i < no_of_nodes; i++) {
        nodes_reached = 1;      /* itself */
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        VECTOR(already_added)[i] = i + 1;

        IGRAPH_PROGRESS("Path length histogram: ", 100.0 * i / no_of_nodes, NULL);

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);

            neis = igraph_adjlist_get(&allneis, actnode);
            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                igraph_integer_t neighbor = VECTOR(*neis)[j];
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
    if (!directed || !igraph_is_directed(graph)) {
        for (i = 0; i < ressize; i++) {
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
