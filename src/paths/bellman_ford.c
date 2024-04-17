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
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_paths.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "core/interruption.h"

/**
 * \function igraph_distances_bellman_ford
 * \brief Weighted shortest path lengths between vertices, allowing negative weights.
 *
 * This function implements the Bellman-Ford algorithm to find the weighted
 * shortest paths to all vertices from a single source, allowing negative weights.
 * It is run independently for the given sources. If there are no negative
 * weights, you are better off with \ref igraph_distances_dijkstra() .
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here, the matrix will be resized if needed.
 *    Each row contains the distances from a single source, to all
 *    vertices in the graph, in the order of vertex IDs. For unreachable
 *    vertices the matrix contains \c IGRAPH_INFINITY.
 * \param from The source vertices.
 * \param to The target vertices.
 * \param weights The edge weights. There must not be any closed loop in
 *    the graph that has a negative total weight (since this would allow
 *    us to decrease the weight of any path containing at least a single
 *    vertex of this loop infinitely). Additionally, no edge weight may
 *    be NaN. If either case does not hold, an error is returned. If this
 *    is a null pointer, then the unweighted version,
 *    \ref igraph_distances() is called.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(s*|E|*|V|), where |V| is the number of
 * vertices, |E| the number of edges and s the number of sources.
 *
 * \sa \ref igraph_distances() for a faster unweighted version
 * or \ref igraph_distances_dijkstra() if you do not have negative
 * edge weights.
 *
 * \example examples/simple/bellman_ford.c
 */
igraph_error_t igraph_distances_bellman_ford(const igraph_t *graph,
                                       igraph_matrix_t *res,
                                       const igraph_vs_t from,
                                       const igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_lazy_inclist_t inclist;
    igraph_integer_t i;
    igraph_integer_t no_of_from, no_of_to;
    igraph_dqueue_int_t Q;
    igraph_vector_bool_t clean_vertices;
    igraph_vector_int_t num_queued;
    igraph_vit_t fromvit, tovit;
    igraph_bool_t all_to;
    igraph_vector_t dist;
    int counter = 0;

    /*
       - speedup: a vertex is marked clean if its distance from the source
         did not change during the last phase. Neighbors of a clean vertex
         are not relaxed again, since it would mean no change in the
         shortest path values. Dirty vertices are queued. Negative loops can
         be detected by checking whether a vertex has been queued at least
         n times.
    */
    if (!weights) {
        return igraph_distances(graph, res, from, to, mode);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }
    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);
    no_of_from = IGRAPH_VIT_SIZE(fromvit);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&Q, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&clean_vertices, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&num_queued, no_of_nodes);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    all_to = igraph_vs_is_all(&to);
    if (all_to) {
        no_of_to = no_of_nodes;
    } else {
        IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
        IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
        no_of_to = IGRAPH_VIT_SIZE(tovit);

        /* No need to check here whether the vertices in 'to' are unique because
         * the loop below uses a temporary distance vector that is then copied
         * into the result matrix at the end of the outer loop iteration, and
         * this is safe even if 'to' contains the same vertex multiple times */
    }

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));

    for (IGRAPH_VIT_RESET(fromvit), i = 0;
         !IGRAPH_VIT_END(fromvit);
         IGRAPH_VIT_NEXT(fromvit), i++) {
        igraph_integer_t source = IGRAPH_VIT_GET(fromvit);

        igraph_vector_fill(&dist, IGRAPH_INFINITY);
        VECTOR(dist)[source] = 0;
        igraph_vector_bool_null(&clean_vertices);
        igraph_vector_int_null(&num_queued);

        /* Fill the queue with vertices to be checked */
        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&Q, j));
        }

        while (!igraph_dqueue_int_empty(&Q)) {
            if (++counter >= 10000) {
                counter = 0;
                IGRAPH_ALLOW_INTERRUPTION();
            }

            igraph_integer_t j = igraph_dqueue_int_pop(&Q);
            VECTOR(clean_vertices)[j] = true;
            VECTOR(num_queued)[j] += 1;
            if (VECTOR(num_queued)[j] > no_of_nodes) {
                IGRAPH_ERROR("Negative loop in graph while calculating distances with Bellman-Ford algorithm.",
                             IGRAPH_ENEGLOOP);
            }

            /* If we cannot get to j in finite time yet, there is no need to relax
             * its edges */
            if (VECTOR(dist)[j] == IGRAPH_INFINITY) {
                continue;
            }

            igraph_vector_int_t *neis = igraph_lazy_inclist_get(&inclist, j);
            IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");

            igraph_integer_t nlen = igraph_vector_int_size(neis);
            for (igraph_integer_t k = 0; k < nlen; k++) {
                igraph_integer_t nei = VECTOR(*neis)[k];
                igraph_integer_t target = IGRAPH_OTHER(graph, nei, j);
                igraph_real_t altdist = VECTOR(dist)[j] + VECTOR(*weights)[nei];
                if (VECTOR(dist)[target] > altdist) {
                    /* relax the edge */
                    VECTOR(dist)[target] = altdist;
                    if (VECTOR(clean_vertices)[target]) {
                        VECTOR(clean_vertices)[target] = false;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, target));
                    }
                }
            }
        }

        /* Copy it to the result */
        if (all_to) {
            igraph_matrix_set_row(res, &dist, i);
        } else {
            igraph_integer_t j;
            for (IGRAPH_VIT_RESET(tovit), j = 0; !IGRAPH_VIT_END(tovit);
                 IGRAPH_VIT_NEXT(tovit), j++) {
                igraph_integer_t v = IGRAPH_VIT_GET(tovit);
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
    igraph_dqueue_int_destroy(&Q);
    igraph_vector_bool_destroy(&clean_vertices);
    igraph_vector_int_destroy(&num_queued);
    igraph_lazy_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_shortest_paths_bellman_ford
 * \brief Weighted shortest path lengths between vertices, allowing negative weights (deprecated).
 *
 * \deprecated-by igraph_distances_bellman_ford 0.10.0
 */
igraph_error_t igraph_shortest_paths_bellman_ford(const igraph_t *graph,
                                       igraph_matrix_t *res,
                                       const igraph_vs_t from,
                                       const igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode) {
    return igraph_distances_bellman_ford(graph, res, from, to, weights, mode);
}

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths_bellman_ford
 * \brief Weighted shortest paths from a vertex, allowing negative weights.
 *
 * This function calculates weighted shortest paths from or to a single vertex
 * using the Bellman-Ford algorithm, whihc can handle negative weights. When
 * there is more than one shortest path between two vertices, only one of them
 * is returned. When there are no negative weights,
 * \ref igraph_get_shortest_paths_dijkstra() is likely to be faster.
 *
 * \param graph The input graph, can be directed.
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
 * \param weights The edge weights. There must not be any closed loop in
 *    the graph that has a negative total weight (since this would allow
 *    us to decrease the weight of any path containing at least a single
 *    vertex of this loop infinitely). If this is a null pointer, then the
 *    unweighted version, \ref igraph_get_shortest_paths() is called.
 *    Edges with positive infinite weights are ignored.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
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
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_ENOMEM
 *           Not enough memory for temporary data.
 *         \cli IGRAPH_EINVAL
 *           The weight vector doesn't math the number of edges.
 *         \cli IGRAPH_EINVVID
 *           \p from is invalid vertex ID
 *         \cli IGRAPH_ENEGLOOP
 *           Bellman-ford algorithm encounted a negative loop.
 *         \endclist
 *
 * Time complexity: O(|E|*|V|), where |V| is the number of
 * vertices, |E| the number of edges.
 *
 * \sa \ref igraph_distances_bellman_ford() to compute only shortest path
 * lengths, but not the paths themselves; \ref igraph_get_shortest_paths() for
 * a faster unweighted version or \ref igraph_get_shortest_paths_dijkstra()
 * if you do not have negative edge weights.
 */

igraph_error_t igraph_get_shortest_paths_bellman_ford(const igraph_t *graph,
                                        igraph_vector_int_list_t *vertices,
                                        igraph_vector_int_list_t *edges,
                                        igraph_integer_t from,
                                        igraph_vs_t to,
                                        const igraph_vector_t *weights,
                                        igraph_neimode_t mode,
                                        igraph_vector_int_t *parents,
                                        igraph_vector_int_t *inbound_edges) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t *parent_eids;
    igraph_lazy_inclist_t inclist;
    igraph_integer_t i, j, k;
    igraph_dqueue_int_t Q;
    igraph_vector_bool_t clean_vertices;
    igraph_vector_int_t num_queued;
    igraph_vit_t tovit;
    igraph_vector_t dist;
    int counter = 0;

    if (!weights) {
        return igraph_get_shortest_paths(graph, vertices, edges, from, to, mode,
                                         parents, inbound_edges);
    }

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("Index of source vertex is out of range.", IGRAPH_EINVVID);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length must match number of edges.", IGRAPH_EINVAL);
    }

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&Q, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&clean_vertices, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&num_queued, no_of_nodes);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
    IGRAPH_FINALLY(igraph_vit_destroy, &tovit);

    if (vertices) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(vertices, IGRAPH_VIT_SIZE(tovit)));
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(edges, IGRAPH_VIT_SIZE(tovit)));
    }

    parent_eids = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(parent_eids, "Insufficient memory for shortest paths with Bellman-Ford.");
    IGRAPH_FINALLY(igraph_free, parent_eids);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    igraph_vector_fill(&dist, IGRAPH_INFINITY);
    VECTOR(dist)[from] = 0;

    /* Fill the queue with vertices to be checked */
    for (j = 0; j < no_of_nodes; j++) {
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, j));
    }

    while (!igraph_dqueue_int_empty(&Q)) {
        if (++counter >= 10000) {
            counter = 0;
            IGRAPH_ALLOW_INTERRUPTION();
        }

        j = igraph_dqueue_int_pop(&Q);
        VECTOR(clean_vertices)[j] = true;
        VECTOR(num_queued)[j] += 1;
        if (VECTOR(num_queued)[j] > no_of_nodes) {
            IGRAPH_ERROR("Negative loop in graph while calculating distances with Bellman-Ford algorithm.",
                         IGRAPH_ENEGLOOP);
        }

        /* If we cannot get to j in finite time yet, there is no need to relax its edges */
        if (VECTOR(dist)[j] == IGRAPH_INFINITY) {
            continue;
        }

        igraph_vector_int_t *neis = igraph_lazy_inclist_get(&inclist, j);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");

        igraph_integer_t nlen = igraph_vector_int_size(neis);
        for (k = 0; k < nlen; k++) {
            igraph_integer_t nei = VECTOR(*neis)[k];
            igraph_integer_t target = IGRAPH_OTHER(graph, nei, j);
            igraph_real_t weight = VECTOR(*weights)[nei];
            igraph_real_t altdist = VECTOR(dist)[j] + weight;

            if (isnan(weight)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            }

            /* infinite weights are handled correctly here; if an edge has
             * infinite weight, altdist will also be infinite so the condition
             * will never be true as if the edge was ignored */

            if (VECTOR(dist)[target] > altdist) {
                /* relax the edge */
                VECTOR(dist)[target] = altdist;
                parent_eids[target] = nei + 1;
                if (VECTOR(clean_vertices)[target]) {
                    VECTOR(clean_vertices)[target] = false;
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, target));
                }
            }
        }
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
        for (IGRAPH_VIT_RESET(tovit), i = 0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit), i++) {
            igraph_integer_t node = IGRAPH_VIT_GET(tovit);
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

    igraph_vector_destroy(&dist);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vit_destroy(&tovit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FREE(parent_eids);
    igraph_dqueue_int_destroy(&Q);
    igraph_vector_bool_destroy(&clean_vertices);
    igraph_vector_int_destroy(&num_queued);
    igraph_lazy_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_shortest_path_bellman_ford
 * \brief Weighted shortest path from one vertex to another one (Bellman-Ford).
 *
 * Finds a weighted shortest path from a single source vertex to
 * a single target using the Bellman-Ford algorithm.
 *
 * </para><para>
 * This function is a special case (and a wrapper) to
 * \ref igraph_get_shortest_paths_bellman_ford().
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param vertices Pointer to an initialized vector or a null
 *        pointer. If not a null pointer, then the vertex IDs along
 *        the path are stored here, including the source and target
 *        vertices.
 * \param edges Pointer to an initialized vector or a null
 *        pointer. If not a null pointer, then the edge IDs along the
 *        path are stored here.
 * \param from The ID of the source vertex.
 * \param to The ID of the target vertex.
 * \param weights The edge weights. There must not be any closed loop in
 *        the graph that has a negative total weight (since this would allow
 *        us to decrease the weight of any path containing at least a single
 *        vertex of this loop infinitely). If this is a null pointer, then the
 *        unweighted version is called.
 * \param mode A constant specifying how edge directions are
 *        considered in directed graphs. \c IGRAPH_OUT follows edge
 *        directions, \c IGRAPH_IN follows the opposite directions,
 *        and \c IGRAPH_ALL ignores edge directions. This argument is
 *        ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|E|log|E|+|V|), |V| is the number of vertices,
 * |E| is the number of edges in the graph.
 *
 * \sa \ref igraph_get_shortest_paths_bellman_ford() for the version with
 * more target vertices.
 */

igraph_error_t igraph_get_shortest_path_bellman_ford(const igraph_t *graph,
                                          igraph_vector_int_t *vertices,
                                          igraph_vector_int_t *edges,
                                          igraph_integer_t from,
                                          igraph_integer_t to,
                                          const igraph_vector_t *weights,
                                          igraph_neimode_t mode) {

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

    IGRAPH_CHECK(igraph_get_shortest_paths_bellman_ford(graph, vp, ep,
                                                        from, igraph_vss_1(to),
                                                        weights, mode, NULL, NULL));

    /* We use the constant time vector_swap() instead of the linear-time vector_update() to move the
       result to the output parameter. */
    if (edges) {
        IGRAPH_CHECK(igraph_vector_int_swap(edges, igraph_vector_int_list_get_ptr(&edges2, 0)));
        igraph_vector_int_list_destroy(&edges2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (vertices) {
        IGRAPH_CHECK(igraph_vector_int_swap(vertices, igraph_vector_int_list_get_ptr(&vertices2, 0)));
        igraph_vector_int_list_destroy(&vertices2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
