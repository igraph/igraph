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

#include "igraph_centrality.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "core/interruption.h"

#include <string.h>

/**
 * \function igraph_convergence_degree
 * \brief Calculates the convergence degree of each edge in a graph.
 *
 * Let us define the input set of an edge (i, j) as the set of vertices where
 * the shortest paths passing through (i, j) originate, and similarly, let us
 * defined the output set of an edge (i, j) as the set of vertices where the
 * shortest paths passing through (i, j) terminate. The convergence degree of
 * an edge is defined as the normalized value of the difference between the
 * size of the input set and the output set, i.e. the difference of them
 * divided by the sum of them. Convergence degrees are in the range (-1, 1); a
 * positive value indicates that the edge is \em convergent since the shortest
 * paths passing through it originate from a larger set and terminate in a
 * smaller set, while a negative value indicates that the edge is \em divergent
 * since the paths originate from a small set and terminate in a larger set.
 *
 * </para><para>
 * Note that the convergence degree as defined above does not make sense in
 * undirected graphs as there is no distinction between the input and output
 * set. Therefore, for undirected graphs, the input and output sets of an edge
 * are determined by orienting the edge arbitrarily while keeping the remaining
 * edges undirected, and then taking the absolute value of the convergence
 * degree.
 *
 * \param graph The input graph, it can be either directed or undirected.
 * \param result Pointer to an initialized vector; the convergence degrees of
 *   each edge will be stored here. May be \c NULL if we are not interested in
 *   the exact convergence degrees.
 * \param ins Pointer to an initialized vector; the size of the input set of
 *   each edge will be stored here. May be \c NULL if we are not interested in
 *   the sizes of the input sets.
 * \param outs Pointer to an initialized vector; the size of the output set of
 *   each edge will be stored here. May be \c NULL if we are not interested in
 *   the sizes of the output sets.
 * \return Error code.
 *
 * Time complexity: O(|V||E|), the number of vertices times the number of edges.
 */
int igraph_convergence_degree(const igraph_t *graph, igraph_vector_t *result,
                              igraph_vector_t *ins, igraph_vector_t *outs) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int i, j, k, n;
    long int *geodist;
    igraph_vector_int_t *eids;
    igraph_vector_t *ins_p, *outs_p, ins_v, outs_v;
    igraph_dqueue_t q;
    igraph_inclist_t inclist;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (result != 0) {
        IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
    }
    IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &q);

    if (ins == 0) {
        ins_p = &ins_v;
        IGRAPH_VECTOR_INIT_FINALLY(ins_p, no_of_edges);
    } else {
        ins_p = ins;
        IGRAPH_CHECK(igraph_vector_resize(ins_p, no_of_edges));
        igraph_vector_null(ins_p);
    }

    if (outs == 0) {
        outs_p = &outs_v;
        IGRAPH_VECTOR_INIT_FINALLY(outs_p, no_of_edges);
    } else {
        outs_p = outs;
        IGRAPH_CHECK(igraph_vector_resize(outs_p, no_of_edges));
        igraph_vector_null(outs_p);
    }

    geodist = IGRAPH_CALLOC(no_of_nodes, long int);
    if (geodist == 0) {
        IGRAPH_ERROR("Cannot calculate convergence degrees", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, geodist);

    /* Collect shortest paths originating from/to every node to correctly
     * determine input and output field sizes */
    for (k = 0; k < (directed ? 2 : 1); k++) {
        igraph_neimode_t neimode = (k == 0) ? IGRAPH_OUT : IGRAPH_IN;
        igraph_real_t *vec;
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, neimode, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
        vec = (k == 0) ? VECTOR(*ins_p) : VECTOR(*outs_p);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_dqueue_clear(&q);
            memset(geodist, 0, sizeof(long int) * (size_t) no_of_nodes);
            geodist[i] = 1;
            IGRAPH_CHECK(igraph_dqueue_push(&q, i));
            IGRAPH_CHECK(igraph_dqueue_push(&q, 0.0));
            while (!igraph_dqueue_empty(&q)) {
                long int actnode = (long int) igraph_dqueue_pop(&q);
                long int actdist = (long int) igraph_dqueue_pop(&q);
                IGRAPH_ALLOW_INTERRUPTION();
                eids = igraph_inclist_get(&inclist, actnode);
                n = igraph_vector_int_size(eids);
                for (j = 0; j < n; j++) {
                    long int neighbor = IGRAPH_OTHER(graph, VECTOR(*eids)[j], actnode);
                    if (geodist[neighbor] != 0) {
                        /* we've already seen this node, another shortest path? */
                        if (geodist[neighbor] - 1 == actdist + 1) {
                            /* Since this edge is in the BFS tree rooted at i, we must
                             * increase either the size of the infield or the outfield */
                            if (!directed) {
                                if (actnode < neighbor) {
                                    VECTOR(*ins_p)[(long int)VECTOR(*eids)[j]] += 1;
                                } else {
                                    VECTOR(*outs_p)[(long int)VECTOR(*eids)[j]] += 1;
                                }
                            } else {
                                vec[(long int)VECTOR(*eids)[j]] += 1;
                            }
                        } else if (geodist[neighbor] - 1 < actdist + 1) {
                            continue;
                        }
                    } else {
                        /* we haven't seen this node yet */
                        IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        /* Since this edge is in the BFS tree rooted at i, we must
                         * increase either the size of the infield or the outfield */
                        if (!directed) {
                            if (actnode < neighbor) {
                                VECTOR(*ins_p)[(long int)VECTOR(*eids)[j]] += 1;
                            } else {
                                VECTOR(*outs_p)[(long int)VECTOR(*eids)[j]] += 1;
                            }
                        } else {
                            vec[(long int)VECTOR(*eids)[j]] += 1;
                        }
                        geodist[neighbor] = actdist + 2;
                    }
                }
            }
        }

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (result != 0) {
        for (i = 0; i < no_of_edges; i++) {
            VECTOR(*result)[i] = (VECTOR(*ins_p)[i] - VECTOR(*outs_p)[i]) /
                                 (VECTOR(*ins_p)[i] + VECTOR(*outs_p)[i]);
        }

        if (!directed) {
            for (i = 0; i < no_of_edges; i++) {
                if (VECTOR(*result)[i] < 0) {
                    VECTOR(*result)[i] = -VECTOR(*result)[i];
                }
            }
        }
    }

    if (ins == 0) {
        igraph_vector_destroy(ins_p);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (outs == 0) {
        igraph_vector_destroy(outs_p);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_FREE(geodist);
    igraph_dqueue_destroy(&q);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}
