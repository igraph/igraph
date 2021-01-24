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

/**
 * \function igraph_shortest_paths_bellman_ford
 * Weighted shortest paths from some sources allowing negative weights.
 *
 * This function is the Bellman-Ford algorithm to find the weighted
 * shortest paths to all vertices from a single source. (It is run
 * independently for the given sources.). If there are no negative
 * weights, you are better off with \ref igraph_shortest_paths_dijkstra() .
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here, the matrix will be resized if needed.
 *    Each row contains the distances from a single source, to all
 *    vertices in the graph, in the order of vertex ids. For unreachable
 *    vertices the matrix contains \c IGRAPH_INFINITY.
 * \param from The source vertices.
 * \param weights The edge weights. There mustn't be any closed loop in
 *    the graph that has a negative total weight (since this would allow
 *    us to decrease the weight of any path containing at least a single
 *    vertex of this loop infinitely). Additionally, no edge weight may
 *    be NaN. If either case does not hold, an error is returned. If this
 *    is a null pointer, then the unweighted version,
 *    \ref igraph_shortest_paths() is called.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(s*|E|*|V|), where |V| is the number of
 * vertices, |E| the number of edges and s the number of sources.
 *
 * \sa \ref igraph_shortest_paths() for a faster unweighted version
 * or \ref igraph_shortest_paths_dijkstra() if you do not have negative
 * edge weights.
 *
 * \example examples/simple/bellman_ford.c
 */
igraph_error_t igraph_shortest_paths_bellman_ford(const igraph_t *graph,
                                       igraph_matrix_t *res,
                                       const igraph_vs_t from,
                                       const igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode) {
    igraph_long_t no_of_nodes = igraph_vcount(graph);
    igraph_long_t no_of_edges = igraph_ecount(graph);
    igraph_lazy_inclist_t inclist;
    igraph_long_t i, j, k;
    igraph_long_t no_of_from, no_of_to;
    igraph_dqueue_t Q;
    igraph_vector_t clean_vertices;
    igraph_vector_t num_queued;
    igraph_vit_t fromvit, tovit;
    igraph_real_t my_infinity = IGRAPH_INFINITY;
    igraph_bool_t all_to;
    igraph_vector_t dist;

    /*
       - speedup: a vertex is marked clean if its distance from the source
         did not change during the last phase. Neighbors of a clean vertex
         are not relaxed again, since it would mean no change in the
         shortest path values. Dirty vertices are queued. Negative loops can
         be detected by checking whether a vertex has been queued at least
         n times.
    */
    if (!weights) {
        return igraph_shortest_paths(graph, res, from, to, mode);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
    }
    if (no_of_edges > 0 && igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);
    no_of_from = IGRAPH_VIT_SIZE(fromvit);

    IGRAPH_DQUEUE_INIT_FINALLY(&Q, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&clean_vertices, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&num_queued, no_of_nodes);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    all_to = igraph_vs_is_all(&to);
    if (all_to) {
        no_of_to = no_of_nodes;
    } else {
        IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
        IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
        no_of_to = IGRAPH_VIT_SIZE(tovit);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));

    for (IGRAPH_VIT_RESET(fromvit), i = 0;
         !IGRAPH_VIT_END(fromvit);
         IGRAPH_VIT_NEXT(fromvit), i++) {
        igraph_long_t source = IGRAPH_VIT_GET(fromvit);

        igraph_vector_fill(&dist, my_infinity);
        VECTOR(dist)[source] = 0;
        igraph_vector_null(&clean_vertices);
        igraph_vector_null(&num_queued);

        /* Fill the queue with vertices to be checked */
        for (j = 0; j < no_of_nodes; j++) {
            IGRAPH_CHECK(igraph_dqueue_push(&Q, j));
        }

        while (!igraph_dqueue_empty(&Q)) {
            igraph_vector_t *neis;
            igraph_long_t nlen;

            j = (igraph_long_t) igraph_dqueue_pop(&Q);
            VECTOR(clean_vertices)[j] = 1;
            VECTOR(num_queued)[j] += 1;
            if (VECTOR(num_queued)[j] > no_of_nodes) {
                IGRAPH_ERROR("cannot run Bellman-Ford algorithm", IGRAPH_ENEGLOOP);
            }

            /* If we cannot get to j in finite time yet, there is no need to relax
             * its edges */
            if (!IGRAPH_FINITE(VECTOR(dist)[j])) {
                continue;
            }

            neis = igraph_lazy_inclist_get(&inclist, (igraph_long_t) j);
            nlen = igraph_vector_size(neis);

            for (k = 0; k < nlen; k++) {
                igraph_long_t nei = (igraph_long_t) VECTOR(*neis)[k];
                igraph_long_t target = IGRAPH_OTHER(graph, nei, j);
                if (VECTOR(dist)[target] > VECTOR(dist)[j] + VECTOR(*weights)[nei]) {
                    /* relax the edge */
                    VECTOR(dist)[target] = VECTOR(dist)[j] + VECTOR(*weights)[nei];
                    if (VECTOR(clean_vertices)[target]) {
                        VECTOR(clean_vertices)[target] = 0;
                        IGRAPH_CHECK(igraph_dqueue_push(&Q, target));
                    }
                }
            }
        }

        /* Copy it to the result */
        if (all_to) {
            igraph_matrix_set_row(res, &dist, i);
        } else {
            for (IGRAPH_VIT_RESET(tovit), j = 0; !IGRAPH_VIT_END(tovit);
                 IGRAPH_VIT_NEXT(tovit), j++) {
                igraph_long_t v = IGRAPH_VIT_GET(tovit);
                MATRIX(*res, i, j) = VECTOR(dist)[v];
            }
        }
    }

    igraph_vector_destroy(&dist);
    IGRAPH_FINALLY_CLEAN(1);

    if (!all_to) {
        igraph_vit_destroy(&tovit);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vit_destroy(&fromvit);
    igraph_dqueue_destroy(&Q);
    igraph_vector_destroy(&clean_vertices);
    igraph_vector_destroy(&num_queued);
    igraph_lazy_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}
