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
#include "igraph_vector_list.h"
#include "igraph_adjlist.h"
#include "igraph_memory.h"

#include "core/indheap.h"
#include "core/interruption.h"

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths_astar
 * \brief Shortest paths from a vertex with heuristic.
 *
 * </para><para>
 * Finds the shortest path between on starting vertex and one or several enpoints.
 * Compared to Dijkstra's algorithm, A* uses a heuristic to decide the direction to move to.
 *
 * </para><para>
 * If there is more than one path with the smallest weight between two vertices, this
 * function gives only one of them.
 * \param graph The graph object.
 * \param vertices The result, the IDs of the vertices along the paths.
 *        This is a list of integer vectors where each element is an
 *        \ref igraph_vector_int_t object. The list will be resized as needed.
 *        Supply a null pointer here if you don't need these vectors.
 * \param edges The result, the IDs of the edges along the paths.
 *        This is a list of integer vectors where each element is an
 *        \ref igraph_vector_int_t object. The list will be resized as needed.
 *        Supply a null pointer here if you don't need these vectors.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the IDs of the vertices to/from which the
 *        shortest paths will be calculated. A vertex might be given multiple
 *        times.
 * \param weights Optional edge weights. Supply \c NULL for unweighted graphs.
 *        all edge weights must be non-negative. Additionally, no
 *        edge weight may be NaN. If either case does not hold, an error
 *        is returned.
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
 * \param parents A pointer to an initialized igraph vector or null.
 *        If not null, a vector containing the parent of each vertex in
 *        the single source shortest path tree is returned here. The
 *        parent of vertex i in the tree is the vertex from which vertex i
 *        was reached. The parent of the start vertex (in the \c from
 *        argument) is -1. If the parent is -2, it means
 *        that the given vertex was not reached from the source during the
 *        search. Note that the search terminates if all the vertices in
 *        \c to are reached.
 * \param inbound_edges A pointer to an initialized igraph vector or null.
 *        If not null, a vector containing the inbound edge of each vertex in
 *        the single source shortest path tree is returned here. The
 *        inbound edge of vertex i in the tree is the edge via which vertex i
 *        was reached. The start vertex and vertices that were not reached
 *        during the search will have -1 in the corresponding entry of the
 *        vector. Note that the search terminates if all the vertices in
 *        \c to are reached.
 * \param heuristic A function that returns an estimate of the distance as
 *        \c igraph_real_t in its first argument. The second parameter is
 *        passed the vertex id as an \c igraph_integer_t, and the third
 *         parameter is passed \p extra.
 * \param extra This is passed on to the heuristic.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p from is invalid vertex ID
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(|E|log|V|+|V|), where |V| is the number of
 * vertices and |E| is the number of edges
 *
 */
igraph_error_t igraph_get_shortest_paths_astar(const igraph_t *graph, 
                                       igraph_vector_int_list_t *vertices,
                                       igraph_vector_int_list_t *edges,
                                       igraph_integer_t from,
                                       igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode,
                                       igraph_vector_int_t *parents,
                                       igraph_vector_int_t *inbound_edges,
                                       igraph_astar_heuristic_t *heuristic,
                                       void *extra
                                      ) {
    igraph_real_t heur_res;

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t dists;
    igraph_integer_t *parent_eids;
    igraph_bool_t *is_target;
    igraph_integer_t i, to_reach;

    if (weights) { //If there are no weights, they are treated as 1.
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERROR("Weight vector must be non-negative.", IGRAPH_EINVAL);
            }
            else if (isnan(min)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values", IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    if (vertices) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(vertices, IGRAPH_VIT_SIZE(vit)));
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(edges, IGRAPH_VIT_SIZE(vit)));
    }


    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
    igraph_vector_fill(&dists, IGRAPH_INFINITY);

    parent_eids = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (parent_eids == 0) {
        IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, parent_eids);
    is_target = IGRAPH_CALLOC(no_of_nodes, igraph_bool_t);
    if (is_target == 0) {
        IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, is_target);

    /* Mark the vertices we need to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (!is_target[ IGRAPH_VIT_GET(vit) ]) {
            is_target[ IGRAPH_VIT_GET(vit) ] = 1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    VECTOR(dists)[from] = 0.0;  /* zero distance */
    IGRAPH_CHECK(heuristic(&heur_res, from, extra));
    igraph_2wheap_push_with_index(&Q, from, -heur_res);

    while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
        igraph_integer_t nlen, minnei;
        /*The sum of the distance and the heuristic is the estimate.
         *The estimates should be on the heap,
         *because the minimum estimate should always be handled next.
         *The value taken off the heap is ignored, we just want the index.
         */
        igraph_2wheap_delete_max_index(&Q, &minnei);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        if (is_target[minnei]) {
            is_target[minnei] = 0;
            to_reach--;
        }

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_lazy_inclist_get(&inclist, minnei);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist;
            if (weights) {
                igraph_real_t weight = VECTOR(*weights)[edge];
                if (weight == IGRAPH_INFINITY) continue;
                altdist = VECTOR(dists)[minnei] + weight;
            } else {
                altdist = VECTOR(dists)[minnei] + 1;
            }
            igraph_real_t curdist = VECTOR(dists)[tto];
            if (curdist == IGRAPH_INFINITY) {
                /* This is the first finite distance */
                VECTOR(dists)[tto] = altdist;
                parent_eids[tto] = edge + 1;
                IGRAPH_CHECK(heuristic(&heur_res, tto, extra));
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist -heur_res));
            } else if (altdist < curdist) {
                /* This is a shorter path */
                VECTOR(dists)[tto] = altdist;
                parent_eids[tto] = edge + 1;
                IGRAPH_CHECK(heuristic(&heur_res, tto, extra));
                igraph_2wheap_modify(&Q, tto, -altdist -heur_res);
            }
        }
    } /* !igraph_2wheap_empty(&Q) */

    if (to_reach > 0) {
        IGRAPH_WARNING("Couldn't reach some vertices");
    }

    /* Create `parents' if needed */
    if (parents) {
        IGRAPH_CHECK(igraph_vector_int_resize(parents, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (i == from) {
                /* i is the start vertex */
                VECTOR(*parents)[i] = -1;
            } else if (parent_eids[i] <= 0) {
                /* i was not reached */
                VECTOR(*parents)[i] = -2;
            } else {
                /* i was reached via the edge with ID = parent_eids[i] - 1 */
                VECTOR(*parents)[i] = IGRAPH_OTHER(graph, parent_eids[i] - 1, i);
            }
        }
    }

    /* Create `inbound_edges' if needed */
    if (inbound_edges) {
        IGRAPH_CHECK(igraph_vector_int_resize(inbound_edges, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (parent_eids[i] <= 0) {
                /* i was not reached */
                VECTOR(*inbound_edges)[i] = -1;
            } else {
                /* i was reached via the edge with ID = parent_eids[i] - 1 */
                VECTOR(*inbound_edges)[i] = parent_eids[i] - 1;
            }
        }
    }

    /* Reconstruct the shortest paths based on vertex and/or edge IDs */
    if (vertices || edges) {
        for (IGRAPH_VIT_RESET(vit), i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
            igraph_integer_t node = IGRAPH_VIT_GET(vit);
            igraph_integer_t size, act, edge;
            igraph_vector_int_t *vvec = 0, *evec = 0;
            if (vertices) {
                vvec = igraph_vector_int_list_get_ptr(vertices, i);
                igraph_vector_int_clear(vvec);
            }
            if (edges) {
                evec = igraph_vector_int_list_get_ptr(edges, i);
                igraph_vector_int_clear(evec);
            }

            IGRAPH_ALLOW_INTERRUPTION();

            size = 0;
            act = node;
            while (parent_eids[act]) {
                size++;
                edge = parent_eids[act] - 1;
                act = IGRAPH_OTHER(graph, edge, act);
            }
            if (vvec && (size > 0 || node == from)) {
                IGRAPH_CHECK(igraph_vector_int_resize(vvec, size + 1));
                VECTOR(*vvec)[size] = node;
            }
            if (evec) {
                IGRAPH_CHECK(igraph_vector_int_resize(evec, size));
            }
            act = node;
            while (parent_eids[act]) {
                edge = parent_eids[act] - 1;
                act = IGRAPH_OTHER(graph, edge, act);
                size--;
                if (vvec) {
                    VECTOR(*vvec)[size] = act;
                }
                if (evec) {
                    VECTOR(*evec)[size] = edge;
                }
            }
        }
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    igraph_vector_destroy(&dists);
    IGRAPH_FREE(is_target);
    IGRAPH_FREE(parent_eids);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(6);
    return IGRAPH_SUCCESS;
}
