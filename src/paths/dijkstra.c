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
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_nongraph.h"
#include "igraph_stack.h"
#include "igraph_vector_ptr.h"

#include "core/indheap.h"
#include "core/interruption.h"

#include <string.h>   /* memset */

/**
 * \function igraph_distances_dijkstra_cutoff
 * \brief Weighted shortest path lengths between vertices, with cutoff.
 *
 * \experimental
 *
 * This function is similar to \ref igraph_distances_dijkstra(), but
 * paths longer than \p cutoff will not be considered.
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here. The matrix will be resized as needed.
 *    Each row contains the distances from a single source, to the
 *    vertices given in the \p to argument.
 *    Vertices that are not reachable within distance \p cutoff will
 *    be assigned distance \c IGRAPH_INFINITY.
 * \param from The source vertices.
 * \param to The target vertices. It is not allowed to include a
 *    vertex twice or more.
 * \param weights The edge weights. All edge weights must be
 *    non-negative for Dijkstra's algorithm to work. Additionally, no
 *    edge weight may be NaN. If either case does not hold, an error
 *    is returned. If this is a null pointer, then the unweighted
 *    version, \ref igraph_distances() is called. Edges with positive infinite
 *    weights are ignored.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \param cutoff The maximal length of paths that will be considered.
 *    When the distance of two vertices is greater than this value,
 *    it will be returned as \c IGRAPH_INFINITY. Negative cutoffs are
 *    treated as infinity.
 * \return Error code.
 *
 * Time complexity: at most O(s |E| log|V| + |V|), where |V| is the number of
 * vertices, |E| the number of edges and s the number of sources. The
 * \p cutoff parameter will limit the number of edges traversed from each
 * source vertex, which reduces the computation time.
 *
 * \sa \ref igraph_distances_cutoff() for a (slightly) faster unweighted
 * version.
 *
 * \example examples/simple/distances.c
 */
igraph_error_t igraph_distances_dijkstra_cutoff(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode,
                                   igraph_real_t cutoff) {

    /* Implementation details. This is the basic Dijkstra algorithm,
       with a binary heap. The heap is indexed, i.e. it stores not only
       the distances, but also which vertex they belong to.

       From now on we use a 2-way heap, so the distances can be queried
       directly from the heap.

       Tricks:
       - The opposite of the distance is stored in the heap, as it is a
         maximum heap and we need a minimum heap.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_vit_t fromvit, tovit;
    igraph_integer_t no_of_from, no_of_to;
    igraph_lazy_inclist_t inclist;
    igraph_integer_t i, j;
    igraph_bool_t all_to;
    igraph_vector_int_t indexv;

    if (!weights) {
        return igraph_distances_cutoff(graph, res, from, to, mode, cutoff);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weights must not be negative, got %g.", IGRAPH_EINVAL, min);
        } else if (isnan(min)) {
            IGRAPH_ERROR("Weights must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);
    no_of_from = IGRAPH_VIT_SIZE(fromvit);

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    all_to = igraph_vs_is_all(&to);
    if (all_to) {
        no_of_to = no_of_nodes;
    } else {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&indexv, no_of_nodes);
        IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
        IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
        no_of_to = IGRAPH_VIT_SIZE(tovit);

        /* We need to check whether the vertices in 'tovit' are unique; this is
         * because the inner while loop of the main algorithm updates the
         * distance matrix whenever a shorter path is encountered from the
         * source vertex 'i' to a target vertex, and we need to be able to
         * map a target vertex to its column in the distance matrix. The mapping
         * is constructed by the loop below */
        for (i = 0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit)) {
            igraph_integer_t v = IGRAPH_VIT_GET(tovit);
            if (VECTOR(indexv)[v]) {
                IGRAPH_ERROR("Target vertex list must not have any duplicates.",
                             IGRAPH_EINVAL);
            }
            VECTOR(indexv)[v] = ++i;
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));
    igraph_matrix_fill(res, IGRAPH_INFINITY);

    for (IGRAPH_VIT_RESET(fromvit), i = 0;
         !IGRAPH_VIT_END(fromvit);
         IGRAPH_VIT_NEXT(fromvit), i++) {

        igraph_integer_t reached = 0;
        igraph_integer_t source = IGRAPH_VIT_GET(fromvit);

        igraph_2wheap_clear(&Q);

        /* Many systems distinguish between +0.0 and -0.0.
         * Since we store negative distances in the heap,
         * we must insert -0.0 in order to get +0.0 as the
         * final distance result. */
        igraph_2wheap_push_with_index(&Q, source, -0.0);

        while (!igraph_2wheap_empty(&Q)) {
            igraph_integer_t minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            igraph_integer_t nlen;

            if (cutoff >= 0 && mindist > cutoff) {
                continue;
            }

            if (all_to) {
                MATRIX(*res, i, minnei) = mindist;
            } else {
                if (VECTOR(indexv)[minnei]) {
                    MATRIX(*res, i, VECTOR(indexv)[minnei] - 1) = mindist;
                    reached++;
                    if (reached == no_of_to) {
                        igraph_2wheap_clear(&Q);
                        break;
                    }
                }
            }

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_lazy_inclist_get(&inclist, minnei);
            IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_integer_t edge = VECTOR(*neis)[j];
                igraph_real_t weight = VECTOR(*weights)[edge];

                /* Optimization: do not follow infinite-weight edges. */
                if (weight == IGRAPH_INFINITY) {
                    continue;
                }

                igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + weight;

                if (! igraph_2wheap_has_elem(&Q, tto)) {
                    /* This is the first non-infinite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
                } else if (igraph_2wheap_has_active(&Q, tto)) {
                    igraph_real_t curdist = -igraph_2wheap_get(&Q, tto);
                    if (altdist < curdist) {
                        /* This is a shorter path */
                        igraph_2wheap_modify(&Q, tto, -altdist);
                    }
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

    } /* !IGRAPH_VIT_END(fromvit) */

    if (!all_to) {
        igraph_vit_destroy(&tovit);
        igraph_vector_int_destroy(&indexv);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    igraph_vit_destroy(&fromvit);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_distances_dijkstra
 * \brief Weighted shortest path lengths between vertices.
 *
 * This function implements Dijkstra's algorithm, which can find
 * the weighted shortest path lengths from a source vertex to all
 * other vertices. This function allows specifying a set of source
 * and target vertices. The algorithm is run independently for each
 * source and the results are retained only for the specified targets.
 * This implementation uses a binary heap for efficiency.
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here. The matrix will be resized as needed.
 *    Each row contains the distances from a single source, to the
 *    vertices given in the \p to argument.
 *    Unreachable vertices have distance \c IGRAPH_INFINITY.
 * \param from The source vertices.
 * \param to The target vertices. It is not allowed to include a
 *    vertex twice or more.
 * \param weights The edge weights. All edge weights must be
 *    non-negative for Dijkstra's algorithm to work. Additionally, no
 *    edge weight may be NaN. If either case does not hold, an error
 *    is returned. If this is a null pointer, then the unweighted
 *    version, \ref igraph_distances() is called.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(s*|E|log|V|+|V|), where |V| is the number of
 * vertices, |E| the number of edges and s the number of sources.
 *
 * \sa \ref igraph_distances() for a (slightly) faster unweighted
 * version or \ref igraph_distances_bellman_ford() for a weighted
 * variant that works in the presence of negative edge weights (but no
 * negative loops)
 *
 * \example examples/simple/distances.c
 */
igraph_error_t igraph_distances_dijkstra(const igraph_t *graph,
                                         igraph_matrix_t *res,
                                         const igraph_vs_t from,
                                         const igraph_vs_t to,
                                         const igraph_vector_t *weights,
                                         igraph_neimode_t mode) {
    return igraph_distances_dijkstra_cutoff(graph, res, from, to, weights, mode, -1);
}

/**
 * \function igraph_shortest_paths_dijkstra
 * \brief Weighted shortest path lengths between vertices (deprecated).
 *
 * \deprecated-by igraph_distances_dijkstra 0.10.0
 */
igraph_error_t igraph_shortest_paths_dijkstra(const igraph_t *graph,
                                       igraph_matrix_t *res,
                                       const igraph_vs_t from,
                                       const igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode) {
    return igraph_distances_dijkstra(graph, res, from, to, weights, mode);
}

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths_dijkstra
 * \brief Weighted shortest paths from a vertex.
 *
 * Finds weighted shortest paths from a single source vertex to the specified
 * sets of target vertices using Dijkstra's algorithm. If there is more than
 * one path with the smallest weight between two vertices, this function gives
 * only one of them. To find all such paths, use
 * \ref igraph_get_all_shortest_paths_dijkstra().
 *
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
* \param weights The edge weights. All edge weights must be
 *       non-negative for Dijkstra's algorithm to work. Additionally, no
 *       edge weight may be NaN. If either case does not hold, an error
 *       is returned. If this is a null pointer, then the unweighted
 *       version, \ref igraph_get_shortest_paths() is called.
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
 * \sa \ref igraph_distances_dijkstra() if you only need the path length but
 * not the paths themselves; \ref igraph_get_shortest_paths() if all edge
 * weights are equal; \ref igraph_get_all_shortest_paths() to find all
 * shortest paths between (source, target) pairs;
 * \ref igraph_get_shortest_paths_bellman_ford() if some edge weighted are
 * negative.
 *
 * \example examples/simple/igraph_get_shortest_paths_dijkstra.c
 */
igraph_error_t igraph_get_shortest_paths_dijkstra(const igraph_t *graph,
                                       igraph_vector_int_list_t *vertices,
                                       igraph_vector_int_list_t *edges,
                                       igraph_integer_t from,
                                       igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode,
                                       igraph_vector_int_t *parents,
                                       igraph_vector_int_t *inbound_edges) {
    /* Implementation details. This is the basic Dijkstra algorithm,
       with a binary heap. The heap is indexed, i.e. it stores not only
       the distances, but also which vertex they belong to. The other
       mapping, i.e. getting the distance for a vertex is not in the
       heap (that would by the double-indexed heap), but in the result
       matrix.

       Dirty tricks:
       - the opposite of the distance is stored in the heap, as it is a
         maximum heap and we need a minimum heap.
       - we don't use IGRAPH_INFINITY in the distance vector during the
         computation, as isfinite() might involve a function call
         and we want to spare that. So we store distance+1.0 instead of
         distance, and zero denotes infinity.
       - `parent_eids' assigns the inbound edge IDs of all vertices in the
         shortest path tree to the vertices. In this implementation, the
         edge ID + 1 is stored, zero means unreachable vertices.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t dists;
    igraph_integer_t *parent_eids;
    igraph_bool_t *is_target;
    igraph_integer_t i, to_reach;

    if (!weights) {
        return igraph_get_shortest_paths(graph, vertices, edges, from, to, mode,
                                         parents, inbound_edges);
    }

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("Index of source vertex is out of range.", IGRAPH_EINVVID);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length does not match number of edges.", IGRAPH_EINVAL);
    }
    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weights must not be negative, got %g.", IGRAPH_EINVAL, min);
        }
        else if (isnan(min)) {
            IGRAPH_ERROR("Weights must not contain NaN values.", IGRAPH_EINVAL);
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
    igraph_vector_fill(&dists, -1.0);

    parent_eids = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(parent_eids, "Insufficient memory for shortest paths with Dijkstra's algorithm.");
    IGRAPH_FINALLY(igraph_free, parent_eids);

    is_target = IGRAPH_CALLOC(no_of_nodes, igraph_bool_t);
    IGRAPH_CHECK_OOM(is_target, "Insufficient memory for shortest paths with Dijkstra's algorithm.");
    IGRAPH_FINALLY(igraph_free, is_target);

    /* Mark the vertices we need to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (!is_target[ IGRAPH_VIT_GET(vit) ]) {
            is_target[ IGRAPH_VIT_GET(vit) ] = true;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    VECTOR(dists)[from] = 0.0;  /* zero distance */
    parent_eids[from] = 0;
    igraph_2wheap_push_with_index(&Q, from, 0);

    while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
        igraph_integer_t nlen, minnei = igraph_2wheap_max_index(&Q);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        if (is_target[minnei]) {
            is_target[minnei] = false;
            to_reach--;
        }

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_lazy_inclist_get(&inclist, minnei);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
            igraph_real_t curdist = VECTOR(dists)[tto];
            if (curdist < 0) {
                /* This is the first finite distance */
                VECTOR(dists)[tto] = altdist;
                parent_eids[tto] = edge + 1;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
            } else if (altdist < curdist) {
                /* This is a shorter path */
                VECTOR(dists)[tto] = altdist;
                parent_eids[tto] = edge + 1;
                igraph_2wheap_modify(&Q, tto, -altdist);
            }
        }
    } /* !igraph_2wheap_empty(&Q) */

    if (to_reach > 0) {
        IGRAPH_WARNING("Couldn't reach some vertices.");
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

/**
 * \function igraph_get_shortest_path_dijkstra
 * \brief Weighted shortest path from one vertex to another one (Dijkstra).
 *
 * Finds a weighted shortest path from a single source vertex to
 * a single target, using Dijkstra's algorithm. If more than one
 * shortest path exists, an arbitrary one is returned.
 *
 * </para><para>
 * This function is a special case (and a wrapper) to
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
 * \param from The ID of the source vertex.
 * \param to The ID of the target vertex.
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

igraph_error_t igraph_get_shortest_path_dijkstra(const igraph_t *graph,
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

    IGRAPH_CHECK(igraph_get_shortest_paths_dijkstra(graph, vp, ep,
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

/**
 * \ingroup structural
 * \function igraph_get_all_shortest_paths_dijkstra
 * \brief All weighted shortest paths (geodesics) from a vertex.
 *
 * \param graph The graph object.
 * \param vertices Pointer to an initialized integer vector list or NULL.
 *   If not NULL, then each vector object contains the vertices along a
 *   shortest path from \p from to another vertex. The vectors are
 *   ordered according to their target vertex: first the shortest
 *   paths to vertex 0, then to vertex 1, etc. No data is included
 *   for unreachable vertices.
 * \param edges Pointer to an initialized integer vector list or NULL. If
 *   not NULL, then each vector object contains the edges along a
 *   shortest path from \p from to another vertex. The vectors are
 *   ordered according to their target vertex: first the shortest
 *   paths to vertex 0, then to vertex 1, etc. No data is included for
 *   unreachable vertices.
 * \param nrgeo Pointer to an initialized igraph_vector_int_t object or
 *   NULL. If not NULL the number of shortest paths from \p from are
 *   stored here for every vertex in the graph. Note that the values
 *   will be accurate only for those vertices that are in the target
 *   vertex sequence (see \p to), since the search terminates as soon
 *   as all the target vertices have been found.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the IDs of the vertices to/from which the
 *        shortest paths will be calculated. A vertex might be given multiple
 *        times.
 * \param weights The edge weights. All edge weights must be
 *       non-negative for Dijkstra's algorithm to work. Additionally, no
 *       edge weight may be NaN. If either case does not hold, an error
 *       is returned. If this is a null pointer, then the unweighted
 *       version, \ref igraph_get_all_shortest_paths() is called.
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
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p from is an invalid vertex ID
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(|E|log|V|+|V|), where |V| is the number of
 * vertices and |E| is the number of edges
 *
 * \sa \ref igraph_distances_dijkstra() if you only need the path
 * length but not the paths themselves, \ref igraph_get_all_shortest_paths()
 * if all edge weights are equal.
 *
 * \example examples/simple/igraph_get_all_shortest_paths_dijkstra.c
 */
igraph_error_t igraph_get_all_shortest_paths_dijkstra(const igraph_t *graph,
        igraph_vector_int_list_t *vertices,
        igraph_vector_int_list_t *edges,
        igraph_vector_int_t *nrgeo,
        igraph_integer_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode) {
    /* Implementation details: see igraph_get_shortest_paths_dijkstra,
       it's basically the same.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t dists;
    igraph_vector_int_t index;
    igraph_vector_int_t order;
    igraph_vector_ptr_t parents, parents_edge;

    unsigned char *is_target; /* uses more than two discrete values, can't be 'bool' */
    igraph_integer_t i, n, to_reach;
    igraph_bool_t free_vertices = false;
    int cmp_result;
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;

    if (!weights) {
        return igraph_get_all_shortest_paths(graph, vertices, edges, nrgeo, from, to, mode);
    }

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("Index of source vertex is out of range.", IGRAPH_EINVVID);
    }

    if (vertices == NULL && nrgeo == NULL && edges == NULL) {
        return IGRAPH_SUCCESS;
    }
    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length does not match number of edges.", IGRAPH_EINVAL);
    }
    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Edge weights must not be negative, got %g.", IGRAPH_EINVAL, min);
        }
        else if (isnan(min)) {
            IGRAPH_ERROR("Weights must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    /* parents stores a vector for each vertex, listing the parent vertices
     * of each vertex in the traversal. Right now we do not use an
     * igraph_vector_int_list_t because that would pre-initialize vectors
     * for all the nodes even if the traversal would involve only a small part
     * of the graph */
    IGRAPH_CHECK(igraph_vector_ptr_init(&parents, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &parents);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&parents, igraph_vector_destroy);

    /* parents_edge stores a vector for each vertex, listing the parent edges
     * of each vertex in the traversal */
    IGRAPH_CHECK(igraph_vector_ptr_init(&parents_edge, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &parents_edge);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&parents_edge, igraph_vector_destroy);

    for (i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *parent_vec, *parent_edge_vec;

        parent_vec = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(parent_vec, "Insufficient memory for all shortest paths with Dijkstra's algorithm.");
        IGRAPH_FINALLY(igraph_free, parent_vec);
        IGRAPH_CHECK(igraph_vector_int_init(parent_vec, 0));
        VECTOR(parents)[i] = parent_vec;
        IGRAPH_FINALLY_CLEAN(1);

        parent_edge_vec = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(parent_edge_vec, "Insufficient memory for all shortest paths with Dijkstra's algorithm.");
        IGRAPH_FINALLY(igraph_free, parent_edge_vec);
        IGRAPH_CHECK(igraph_vector_int_init(parent_edge_vec, 0));
        VECTOR(parents_edge)[i] = parent_edge_vec;
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* distance of each vertex from the root */
    IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
    igraph_vector_fill(&dists, -1.0);

    /* order lists the order of vertices in which they were found during
     * the traversal */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, 0);

    /* boolean array to mark whether a given vertex is a target or not */
    is_target = IGRAPH_CALLOC(no_of_nodes, unsigned char);
    IGRAPH_CHECK_OOM(is_target, "Insufficient memory for all shortest paths with Dijkstra's algorithm.");
    IGRAPH_FINALLY(igraph_free, is_target);

    /* two-way heap storing vertices and distances */
    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);

    /* lazy adjacency edge list to query neighbours efficiently */
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    /* Mark the vertices we need to reach */
    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (!is_target[ IGRAPH_VIT_GET(vit) ]) {
            is_target[ IGRAPH_VIT_GET(vit) ] = 1;
        } else {
            to_reach--; /* this node was given multiple times */
        }
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    VECTOR(dists)[from] = 0.0; /* zero distance */
    igraph_2wheap_push_with_index(&Q, from, 0.0);

    while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
        igraph_integer_t nlen, minnei = igraph_2wheap_max_index(&Q);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        if (is_target[minnei]) {
            is_target[minnei] = 0;
            to_reach--;
        }

        /* Mark that we have reached this vertex */
        IGRAPH_CHECK(igraph_vector_int_push_back(&order, minnei));

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_lazy_inclist_get(&inclist, minnei);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
            igraph_real_t curdist = VECTOR(dists)[tto];
            igraph_vector_int_t *parent_vec, *parent_edge_vec;

            cmp_result = igraph_cmp_epsilon(curdist, altdist, eps);
            if (curdist < 0) {
                /* This is the first non-infinite distance */
                VECTOR(dists)[tto] = altdist;

                parent_vec = (igraph_vector_int_t*)VECTOR(parents)[tto];
                IGRAPH_CHECK(igraph_vector_int_push_back(parent_vec, minnei));
                parent_edge_vec = (igraph_vector_int_t*)VECTOR(parents_edge)[tto];
                IGRAPH_CHECK(igraph_vector_int_push_back(parent_edge_vec, edge));

                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
            } else if (cmp_result == 0 /* altdist == curdist */ && VECTOR(*weights)[edge] > 0) {
                /* This is an alternative path with exactly the same length.
                 * Note that we consider this case only if the edge via which we
                 * reached the node has a nonzero weight; otherwise we could create
                 * infinite loops in undirected graphs by traversing zero-weight edges
                 * back-and-forth */
                parent_vec = (igraph_vector_int_t*) VECTOR(parents)[tto];
                IGRAPH_CHECK(igraph_vector_int_push_back(parent_vec, minnei));
                parent_edge_vec = (igraph_vector_int_t*) VECTOR(parents_edge)[tto];
                IGRAPH_CHECK(igraph_vector_int_push_back(parent_edge_vec, edge));
            } else if (cmp_result > 0 /* altdist < curdist */) {
                /* This is a shorter path */
                VECTOR(dists)[tto] = altdist;

                parent_vec = (igraph_vector_int_t*)VECTOR(parents)[tto];
                igraph_vector_int_clear(parent_vec);
                IGRAPH_CHECK(igraph_vector_int_push_back(parent_vec, minnei));
                parent_edge_vec = (igraph_vector_int_t*)VECTOR(parents_edge)[tto];
                igraph_vector_int_clear(parent_edge_vec);
                IGRAPH_CHECK(igraph_vector_int_push_back(parent_edge_vec, edge));

                igraph_2wheap_modify(&Q, tto, -altdist);
            }
        }
    } /* !igraph_2wheap_empty(&Q) */

    if (to_reach > 0) {
        IGRAPH_WARNING("Couldn't reach some of the requested target vertices.");
    }

    /* we don't need these anymore */
    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    /*
    printf("Order:\n");
    igraph_vector_int_print(&order);

    printf("Parent vertices:\n");
    for (i = 0; i < no_of_nodes; i++) {
      if (igraph_vector_int_size(VECTOR(parents)[i]) > 0) {
        printf("[%ld]: ", i);
        igraph_vector_int_print(VECTOR(parents)[i]);
      }
    }
    */

    if (nrgeo) {
        IGRAPH_CHECK(igraph_vector_int_resize(nrgeo, no_of_nodes));
        igraph_vector_int_null(nrgeo);

        /* Theoretically, we could calculate nrgeo in parallel with the traversal.
         * However, that way we would have to check whether nrgeo is null or not
         * every time we want to update some element in nrgeo. Since we need the
         * order vector anyway for building the final result, we could just as well
         * build nrgeo here.
         */
        VECTOR(*nrgeo)[from] = 1;
        n = igraph_vector_int_size(&order);
        for (i = 1; i < n; i++) {
            igraph_integer_t node, j, k;
            igraph_vector_int_t *parent_vec;

            node = VECTOR(order)[i];
            /* now, take the parent vertices */
            parent_vec = (igraph_vector_int_t*)VECTOR(parents)[node];
            k = igraph_vector_int_size(parent_vec);
            for (j = 0; j < k; j++) {
                VECTOR(*nrgeo)[node] += VECTOR(*nrgeo)[VECTOR(*parent_vec)[j]];
            }
        }
    }

    if (vertices || edges) {
        igraph_vector_int_t *path, *parent_vec, *parent_edge_vec;
        igraph_vector_t *paths_index;
        igraph_stack_int_t stack;
        igraph_integer_t j, node;

        /* a shortest path from the starting vertex to vertex i can be
         * obtained by calculating the shortest paths from the "parents"
         * of vertex i in the traversal. Knowing which of the vertices
         * are "targets" (see is_target), we can collect for which other
         * vertices do we need to calculate the shortest paths. We reuse
         * is_target for that; is_target = 0 means that we don't need the
         * vertex, is_target = 1 means that the vertex is a target (hence
         * we need it), is_target = 2 means that the vertex is not a target
         * but it stands between a shortest path between the root and one
         * of the targets
         */
        if (igraph_vs_is_all(&to)) {
            memset(is_target, 1, sizeof(unsigned char) * (size_t) no_of_nodes);
        } else {
            memset(is_target, 0, sizeof(unsigned char) * (size_t) no_of_nodes);

            IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
            IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);

            /* Add the target vertices to the queue */
            IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
            IGRAPH_FINALLY(igraph_vit_destroy, &vit);
            for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
                i = IGRAPH_VIT_GET(vit);
                if (!is_target[i]) {
                    is_target[i] = 1;
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, i));
                }
            }
            igraph_vit_destroy(&vit);
            IGRAPH_FINALLY_CLEAN(1);

            while (!igraph_stack_int_empty(&stack)) {
                /* For each parent of node i, get its parents */
                igraph_integer_t el = igraph_stack_int_pop(&stack);
                parent_vec = (igraph_vector_int_t*) VECTOR(parents)[el];
                i = igraph_vector_int_size(parent_vec);

                for (j = 0; j < i; j++) {
                    /* For each parent, check if it's already in the stack.
                     * If not, push it and mark it in is_target */
                    n = VECTOR(*parent_vec)[j];
                    if (!is_target[n]) {
                        is_target[n] = 2;
                        IGRAPH_CHECK(igraph_stack_int_push(&stack, n));
                    }
                }
            }
            igraph_stack_int_destroy(&stack);
            IGRAPH_FINALLY_CLEAN(1);
        }

        /* now, reconstruct the shortest paths from the parent list in the
         * order we've found the nodes during the traversal.
         * dists is being re-used as a vector where element i tells the
         * index in vertices where the shortest paths leading to vertex i
         * start, plus one (so that zero means that there are no paths
         * for a given vertex).
         */
        paths_index = &dists;
        n = igraph_vector_int_size(&order);
        igraph_vector_null(paths_index);

        if (edges) {
            igraph_vector_int_list_clear(edges);
        }

        if (vertices) {
            igraph_vector_int_list_clear(vertices);
        } else {
            /* If the 'vertices' vector doesn't exist, then create one, in order
             * for the algorithm to work. */
            vertices = IGRAPH_CALLOC(1, igraph_vector_int_list_t);
            IGRAPH_CHECK_OOM(vertices, "Insufficient memory for all shortest paths with Dijkstra's algorithm.");
            IGRAPH_FINALLY(igraph_free, vertices);
            IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(vertices, 0);
            free_vertices = true;
        }

        /* by definition, the shortest path leading to the starting vertex
         * consists of the vertex itself only */
        IGRAPH_CHECK(igraph_vector_int_list_push_back_new(vertices, &path));
        IGRAPH_CHECK(igraph_vector_int_push_back(path, from));

        if (edges) {
            /* the shortest path from the source to itself is empty */
            IGRAPH_CHECK(igraph_vector_int_list_push_back_new(edges, &path));
        }
        VECTOR(*paths_index)[from] = 1;

        for (i = 1; i < n; i++) {
            igraph_integer_t m, path_count;
            igraph_vector_int_t *parent_path, *parent_path_edge;

            node = VECTOR(order)[i];

            /* if we don't need the shortest paths for this node (because
             * it is not standing in a shortest path between the source
             * node and any of the target nodes), skip it */
            if (!is_target[node]) {
                continue;
            }

            IGRAPH_ALLOW_INTERRUPTION();

            /* we are calculating the shortest paths of node now. */
            /* first, we update the paths_index */
            path_count = igraph_vector_int_list_size(vertices);
            VECTOR(*paths_index)[node] = path_count + 1;

            /* now, take the parent vertices */
            parent_vec = (igraph_vector_int_t*) VECTOR(parents)[node];
            parent_edge_vec = (igraph_vector_int_t*) VECTOR(parents_edge)[node];
            m = igraph_vector_int_size(parent_vec);

            /*
            printf("Calculating shortest paths to vertex %ld\n", node);
            printf("Parents are: ");
            igraph_vector_print(parent_vec);
            */

            for (j = 0; j < m; j++) {
                /* for each parent, copy the shortest paths leading to that parent
                 * and add the current vertex in the end */
                igraph_integer_t parent_node = VECTOR(*parent_vec)[j];
                igraph_integer_t parent_edge = VECTOR(*parent_edge_vec)[j];
                igraph_integer_t parent_path_idx = VECTOR(*paths_index)[parent_node] - 1;
                /*
                printf("  Considering parent: %ld\n", parent_node);
                printf("  Paths to parent start at index %ld in vertices\n", parent_path_idx);
                */
                IGRAPH_ASSERT(parent_path_idx >= 0);
                for (; parent_path_idx < path_count; parent_path_idx++) {
                    parent_path = igraph_vector_int_list_get_ptr(vertices, parent_path_idx);
                    if (igraph_vector_int_tail(parent_path) != parent_node) {
                        break;
                    }

                    IGRAPH_CHECK(igraph_vector_int_list_push_back_new(vertices, &path));

                    /* We need to re-read parent_path because the previous push_back_new()
                     * call might have reallocated the entire vector list */
                    parent_path = igraph_vector_int_list_get_ptr(vertices, parent_path_idx);
                    IGRAPH_CHECK(igraph_vector_int_update(path, parent_path));
                    IGRAPH_CHECK(igraph_vector_int_push_back(path, node));

                    if (edges) {
                        IGRAPH_CHECK(igraph_vector_int_list_push_back_new(edges, &path));
                        if (parent_node != from) {
                            parent_path_edge = igraph_vector_int_list_get_ptr(edges, parent_path_idx);
                            IGRAPH_CHECK(igraph_vector_int_update(path, parent_path_edge));
                        }
                        IGRAPH_CHECK(igraph_vector_int_push_back(path, parent_edge));
                    }
                }
            }
        }

        /* free those paths from the result vector that we won't need */
        n = igraph_vector_int_list_size(vertices);
        i = 0;
        while (i < n) {
            igraph_integer_t tmp;
            path = igraph_vector_int_list_get_ptr(vertices, i);
            tmp = igraph_vector_int_tail(path);
            if (is_target[tmp] == 1) {
                /* we need this path, keep it */
                i++;
            } else {
                /* we don't need this path, free it */
                igraph_vector_int_list_discard_fast(vertices, i);
                if (edges) {
                    igraph_vector_int_list_discard_fast(edges, i);
                }
                n--;
            }
        }

        /* sort the remaining paths by the target vertices */
        IGRAPH_VECTOR_INT_INIT_FINALLY(&index, 0);
        igraph_vector_int_list_sort_ind(vertices, &index, igraph_vector_int_colex_cmp);
        IGRAPH_CHECK(igraph_vector_int_list_permute(vertices, &index));
        if (edges) {
            IGRAPH_CHECK(igraph_vector_int_list_permute(edges, &index));
        }
        igraph_vector_int_destroy(&index);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* free the allocated memory */
    if (free_vertices) {
        igraph_vector_int_list_destroy(vertices);
        IGRAPH_FREE(vertices);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_int_destroy(&order);
    IGRAPH_FREE(is_target);
    igraph_vector_destroy(&dists);
    igraph_vector_ptr_destroy_all(&parents);
    igraph_vector_ptr_destroy_all(&parents_edge);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}
