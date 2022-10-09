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

#include "igraph_interface.h"
/**
 * \function igraph_get_shortest_path_dijkstra
 * \brief Weighted shortest path from one vertex to another one.
 *
 * Calculates a single (positively) weighted shortest path from
 * a single vertex to another one, using Dijkstra's algorithm.
 *
 * </para><para>This function is a special case (and a wrapper) to
 * \ref igraph_get_shortest_paths_dijkstra().
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param vertices Pointer to an initialized vector or a null
 *        pointer. If not a null pointer, then the vertex IDs along
 *        the path are stored here, including the source and target
 *        vertices.
 * \param edges Pointer to an initialized vector or a null
 *        pointer. If not a null pointer, then the edge IDs along the
 *        path are stored here.
 * \param from The id of the source vertex.
 * \param to The id of the target vertex.
 * \param weights The edge weights. All edge weights must be
 *       non-negative for Dijkstra's algorithm to work. Additionally, no
 *       edge weight may be NaN. If either case does not hold, an error
 *       is returned. If this is a null pointer, then the unweighted
 *       version, \ref igraph_get_shortest_paths() is called.
 * \param mode A constant specifying how edge directions are
 *        considered in directed graphs. \c IGRAPH_OUT follows edge
 *        directions, \c IGRAPH_IN follows the opposite directions,
 *        and \c IGRAPH_ALL ignores edge directions. This argument is
 *        ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|E|log|V|+|V|), |V| is the number of vertices,
 * |E| is the number of edges in the graph.
 *
 * \sa \ref igraph_get_shortest_paths_dijkstra() for the version with
 * more target vertices.
 */

//this can hold any static information in extra, and vertex_id seems like all the dynamic information we need
typedef igraph_error_t igraph_astar_heuristic_t(igraph_real_t *result, igraph_integer_t vertex_id, void *extra);


igraph_error_t igraph_get_shortest_path_astar(const igraph_t *graph,
                                      igraph_vector_int_t *vertices,
                                      igraph_vector_int_t *edges,
                                      igraph_integer_t from,
                                      igraph_integer_t to,
                                      const igraph_vector_t *weights,
                                      igraph_neimode_t mode,
                                      igraph_astar_heuristic_t *heuristic,
                                      //extra is needed so users can pass in information to the heuristic
                                      //maybe the whole graph? and also the weights? or a location matrix?
                                      void *extra
                                      ) {
    //mostly copied from igraph_get_shortest_path_dijkstra
    //and igraph_get_shortest_paths_dijkstra
 
    igraph_real_t heur_res;
    igraph_bool_t found = false;

    igraph_vector_int_list_t vertices2, *vp = &vertices2;
    igraph_vector_int_list_t edges2, *ep = &edges2;

    if (vertices) {
        IGRAPH_CHECK(igraph_vector_int_list_init(&vertices2, 1));
        IGRAPH_FINALLY(igraph_vector_int_list_destroy, &vertices2);
    } else {
        vp = NULL;
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_int_list_init(&edges2, 1));
        IGRAPH_FINALLY(igraph_vector_int_list_destroy, &edges2);
    } else {
        ep = NULL;
    }

    /*
    IGRAPH_CHECK(igraph_get_shortest_paths_dijkstra(graph, vp, ep,
                 from, igraph_vss_1(to),
                 weights, mode, NULL, NULL));
                 */

    /*just copying stuff from get_shortest_paths that seems to be useful.
      this needs to be fixed using #2221 */
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t dists;
    igraph_integer_t *parent_eids;
    igraph_bool_t *is_target;
    igraph_integer_t i, to_reach;

    if (weights) { //for now if there are now weights, they are treated as 1
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERROR("Weight vector must be non-negative", IGRAPH_EINVAL);
            }
            else if (igraph_is_nan(min)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values", IGRAPH_EINVAL);
            }
        }
    }

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

    //maybe we also want an estimates vector here? or we can calculate them when they're needed from the heuristic
    IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
    igraph_vector_fill(&dists, -1.0);

    //we don't need to mark any vertices we need to reach, because there's only one, so we're skipping some code here
    
    VECTOR(dists)[from] = 0.0;  /* zero distance */
    //igraph_2wheap_push_with_index(&Q, from, 0);
    IGRAPH_CHECK(heuristic(&heur_res, from, extra));
    igraph_2wheap_push_with_index(&Q, from, -heur_res);

    while (!igraph_2wheap_empty(&Q)) {
        igraph_integer_t nlen, minnei = igraph_2wheap_max_index(&Q);
        //so there's the distance, the heuristic, and the sum of the distance and the heuristic, which is the estimate. The estimate should be on the heap, because the minimum estimate should always be handled next
        igraph_real_t minest = -igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        if (minnei == to) {
            found = true;
            break;
        }

        /* Now check all neighbors of 'minnei' for a shorter path */
        //distance updates should not use the heuristic, so the old code should be fine, except
        // of course for the heap updates
        neis = igraph_lazy_inclist_get(&inclist, minnei);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist;
            if (weights) {
                altdist = mindist + VECTOR(*weights)[edge];
            } else {
                altdist = mindist + 1;
            }
            igraph_real_t curdist = VECTOR(dists)[tto];
            if (curdist < 0) {
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

    if (!found) {
        IGRAPH_WARNING("Couldn't reach the vertex");
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
    if (edges) {
        IGRAPH_CHECK(igraph_vector_int_update(edges, igraph_vector_int_list_get_ptr(&edges2, 0)));
        igraph_vector_int_list_destroy(&edges2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (vertices) {
        IGRAPH_CHECK(igraph_vector_int_update(vertices, igraph_vector_int_list_get_ptr(&vertices2, 0)));
        igraph_vector_int_list_destroy(&vertices2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
