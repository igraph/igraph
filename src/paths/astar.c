/*
   IGraph library.
   Copyright (C) 2022 The igraph development team

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
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_memory.h"

#include "core/indheap.h"
#include "core/interruption.h"

static igraph_error_t null_heuristic(
    igraph_real_t *result, igraph_integer_t from, igraph_integer_t to,
	void *extra
) {
    IGRAPH_UNUSED(from);
    IGRAPH_UNUSED(to);
    IGRAPH_UNUSED(extra);

    *result = 0;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_shortest_path_astar
 * \brief A* gives the shortest path from one vertex to another, with heuristic.
 *
 * \experimental
 *
 * Calculates a shortest path from a single source vertex to a single
 * target, using the A* algorithm. A* tries to find a shortest path by
 * starting at \p from and moving to vertices that lie on a path with
 * the lowest estimated length. This length estimate is the sum of two
 * numbers: the distance from the source (\p from) to the intermediate vertex,
 * and the value returned by the heuristic function. The heuristic function
 * provides an estimate the distance between intermediate candidate
 * vertices and the target vertex \p to. The A* algorithm is guaranteed
 * to give the correct shortest path (if one exists) only if the heuristic
 * does not overestimate distances, i.e. if the heuristic function is
 * \em admissible.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param vertices Pointer to an initialized vector or the \c NULL
 *        pointer. If not \c NULL, then the vertex IDs along
 *        the path are stored here, including the source and target
 *        vertices.
 * \param edges Pointer to an initialized vector or the \c NULL
 *        pointer. If not \c NULL, then the edge IDs along the
 *        path are stored here.
 * \param from The ID of the source vertex.
 * \param to The ID of the target vertex.
 * \param weights Optional edge weights. Supply \c NULL for unweighted graphs.
 *        All edge weights must be non-negative. Additionally, no
 *        edge weight may be NaN. If either case does not hold, an error
 *        is returned. Edges with positive infinite weights are ignored.
 * \param mode A constant specifying how edge directions are
 *        considered in directed graphs. \c IGRAPH_OUT follows edge
 *        directions, \c IGRAPH_IN follows the opposite directions,
 *        and \c IGRAPH_ALL ignores edge directions. This argument is
 *        ignored for undirected graphs.
 * \param heuristic A function that provides distance estimates to the
 *        target vertex. See \ref igraph_astar_heuristic_func_t for
 *        more information.
 * \param extra This is passed on to the heuristic function.
 * \return Error code.
 *
 * Time complexity: In the worst case, O(|E|log|V|+|V|), where
 * |V| is the number of vertices and
 * |E| is the number of edges in the graph.
 * The running time depends on the accuracy of the distance estimates
 * returned by the heuristic function. Assuming that the heuristic
 * is admissible, the better the estimates, the shortert the running
 * time.
 */

igraph_error_t igraph_get_shortest_path_astar(const igraph_t *graph,
                                      igraph_vector_int_t *vertices,
                                      igraph_vector_int_t *edges,
                                      igraph_integer_t from,
                                      igraph_integer_t to,
                                      const igraph_vector_t *weights,
                                      igraph_neimode_t mode,
                                      igraph_astar_heuristic_func_t *heuristic,
                                      void *extra)
{
    igraph_real_t heur_res;

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t dists;
    igraph_integer_t *parent_eids;

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("Starting vertex out of range.", IGRAPH_EINVVID);
    }

    if (to < 0 || to >= no_of_nodes) {
        IGRAPH_ERROR("End vertex out of range.", IGRAPH_EINVVID);
    }

    if (!heuristic) {
        heuristic = null_heuristic;
    }

    if (weights) { /* If there are no weights, they are treated as 1. */
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                          IGRAPH_EINVAL, igraph_vector_size(weights), no_of_edges);
        }
        if (no_of_edges > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERRORF("Weight vector must be non-negative, found weight of %g.", IGRAPH_EINVAL, min);
            }
            else if (isnan(min)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    /* dists[v] is the length of the shortest path found so far between 'from' and 'v'. */
    IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
    igraph_vector_fill(&dists, IGRAPH_INFINITY);

    /* parent_eids[v] is the 1 + the ID of v's inbound edge in the shortest path tree.
     * A value of 0 indicates unreachable vertices. */
    parent_eids = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(parent_eids, "Insufficient memory for shortest paths with A* algorithm.");
    IGRAPH_FINALLY(igraph_free, parent_eids);

    VECTOR(dists)[from] = 0.0;
    IGRAPH_CHECK(heuristic(&heur_res, from, to, extra));
    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, from, -heur_res));

    igraph_bool_t found = false;
    while (!igraph_2wheap_empty(&Q)) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* The from -> u -> to distance estimate is the sum of the
         * from -> u distance and the u -> to distance estimate
         * obtained from the heuristic.
         *
         * We use an indexed heap to process 'u' vertices in order
         * of the smallest from -> u -> to distance estimate. Since
         * we only have a maximum heap available, we store negated values
         * in order to obtain smallest values first. The value taken off
         * the heap is ignored, we just want the index of 'u'. */

        igraph_integer_t u;
        igraph_2wheap_delete_max_index(&Q, &u);

        /* Reached the target vertex, the search can be stopped. */
        if (u == to) {
            found = true;
            break;
        }

        /* Now we check all neighbors 'u' for a path with a shorter actual (not estimated)
         * length than what was found so far. */

        igraph_vector_int_t *neis = igraph_lazy_inclist_get(&inclist, u);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");

        igraph_integer_t nlen = igraph_vector_int_size(neis);
        for (igraph_integer_t i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t v = IGRAPH_OTHER(graph, edge, u);
            igraph_real_t altdist; /* candidate from -> v distance */
            if (weights) {
                igraph_real_t weight = VECTOR(*weights)[edge];
                if (weight == IGRAPH_INFINITY) {
                    continue;
                }
                altdist = VECTOR(dists)[u] + weight;
            } else {
                altdist = VECTOR(dists)[u] + 1;
            }
            igraph_real_t curdist = VECTOR(dists)[v];
            if (curdist == IGRAPH_INFINITY) {
                /* This is the first finite from -> v distance we found.
                 * Here we rely on infinite weight edges having been skipped, see TODO above. */
                VECTOR(dists)[v] = altdist;
                parent_eids[v] = edge + 1;
                IGRAPH_CHECK(heuristic(&heur_res, v, to, extra));
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, v, -(altdist + heur_res)));
            } else if (altdist < curdist) {
                /* This is a shorter from -> v path than what was found before. */
                VECTOR(dists)[v] = altdist;
                parent_eids[v] = edge + 1;
                IGRAPH_CHECK(heuristic(&heur_res, v, to, extra));
                igraph_2wheap_modify(&Q, v, -(altdist + heur_res));
            }
        }
    } /* !igraph_2wheap_empty(&Q) */

    if (!found) {
        IGRAPH_WARNING("Couldn't reach the target vertex.");
    }

    /* Reconstruct the shortest paths based on vertex and/or edge IDs */
    if (vertices || edges) {
        igraph_integer_t size, act, edge;

        if (vertices) {
            igraph_vector_int_clear(vertices);
        }
        if (edges) {
            igraph_vector_int_clear(edges);
        }

        IGRAPH_ALLOW_INTERRUPTION();

        size = 0;
        act = to;
        while (parent_eids[act]) {
            size++;
            edge = parent_eids[act] - 1;
            act = IGRAPH_OTHER(graph, edge, act);
        }
        if (vertices && (size > 0 || to == from)) {
            IGRAPH_CHECK(igraph_vector_int_resize(vertices, size + 1));
            VECTOR(*vertices)[size] = to;
        }
        if (edges) {
            IGRAPH_CHECK(igraph_vector_int_resize(edges, size));
        }
        act = to;
        while (parent_eids[act]) {
            edge = parent_eids[act] - 1;
            act = IGRAPH_OTHER(graph, edge, act);
            size--;
            if (vertices) {
                VECTOR(*vertices)[size] = act;
            }
            if (edges) {
                VECTOR(*edges)[size] = edge;
            }
        }
    }


    IGRAPH_FREE(parent_eids);
    igraph_vector_destroy(&dists);
    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}
