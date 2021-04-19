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
#include "igraph_memory.h"
#include "igraph_stack.h"

#include "core/indheap.h"
#include "core/interruption.h"

#include <string.h>   /* memset */

/**
 * \function igraph_shortest_paths_dijkstra
 * \brief Weighted shortest path lengths between vertices.
 *
 * This function implements Dijkstra's algorithm to find the weighted
 * shortest path lengths to all vertices from a single source. It is run
 * independently for the given sources. It uses a binary heap for
 * efficient implementation.
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here. The matrix will be resized as needed.
 *    Each row contains the distances from a single source, to the
 *    vertices given in the \c to argument.
 *    Unreachable vertices has distance
 *    \c IGRAPH_INFINITY.
 * \param from The source vertices.
 * \param to The target vertices. It is not allowed to include a
 *    vertex twice or more.
 * \param weights The edge weights. All edge weights must be
 *    non-negative for Dijkstra's algorithm to work. Additionally, no
 *    edge weight may be NaN. If either case does not hold, an error
 *    is returned. If this is a null pointer, then the unweighted
 *    version, \ref igraph_shortest_paths() is called.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(s*|E|log|E|+|V|), where |V| is the number of
 * vertices, |E| the number of edges and s the number of sources.
 *
 * \sa \ref igraph_shortest_paths() for a (slightly) faster unweighted
 * version or \ref igraph_shortest_paths_bellman_ford() for a weighted
 * variant that works in the presence of negative edge weights (but no
 * negative loops).
 *
 * \example examples/simple/dijkstra.c
 */
int igraph_shortest_paths_dijkstra(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    /* Implementation details. This is the basic Dijkstra algorithm,
       with a binary heap. The heap is indexed, i.e. it stores not only
       the distances, but also which vertex they belong to.

       From now on we use a 2-way heap, so the distances can be queried
       directly from the heap.

       Dirty tricks:
       - the opposite of the distance is stored in the heap, as it is a
         maximum heap and we need a minimum heap.
       - we don't use IGRAPH_INFINITY in the res matrix during the
         computation, as IGRAPH_FINITE() might involve a function call
         and we want to spare that. -1 will denote infinity instead.
    */

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_vit_t fromvit, tovit;
    long int no_of_from, no_of_to;
    igraph_lazy_inclist_t inclist;
    long int i, j;
    igraph_real_t my_infinity = IGRAPH_INFINITY;
    igraph_bool_t all_to;
    igraph_vector_t indexv;

    if (!weights) {
        return igraph_shortest_paths(graph, res, from, to, mode);
    }

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
        IGRAPH_VECTOR_INIT_FINALLY(&indexv, no_of_nodes);
        IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
        IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
        no_of_to = IGRAPH_VIT_SIZE(tovit);
        for (i = 0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit)) {
            long int v = IGRAPH_VIT_GET(tovit);
            if (VECTOR(indexv)[v]) {
                IGRAPH_ERROR("Duplicate vertices in `to', this is not allowed",
                             IGRAPH_EINVAL);
            }
            VECTOR(indexv)[v] = ++i;
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));
    igraph_matrix_fill(res, my_infinity);

    for (IGRAPH_VIT_RESET(fromvit), i = 0;
         !IGRAPH_VIT_END(fromvit);
         IGRAPH_VIT_NEXT(fromvit), i++) {

        long int reached = 0;
        long int source = IGRAPH_VIT_GET(fromvit);
        igraph_2wheap_clear(&Q);
        igraph_2wheap_push_with_index(&Q, source, -1.0);

        while (!igraph_2wheap_empty(&Q)) {
            long int minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            long int nlen;

            if (all_to) {
                MATRIX(*res, i, minnei) = mindist - 1.0;
            } else {
                if (VECTOR(indexv)[minnei]) {
                    MATRIX(*res, i, (long int)(VECTOR(indexv)[minnei] - 1)) = mindist - 1.0;
                    reached++;
                    if (reached == no_of_to) {
                        igraph_2wheap_clear(&Q);
                        break;
                    }
                }
            }

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_lazy_inclist_get(&inclist, (igraph_integer_t) minnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                long int tto = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
                igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
                igraph_real_t curdist = active ? -igraph_2wheap_get(&Q, tto) : 0.0;
                if (!has) {
                    /* This is the first non-infinite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
                } else if (altdist < curdist) {
                    /* This is a shorter path */
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, -altdist));
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

    } /* !IGRAPH_VIT_END(fromvit) */

    if (!all_to) {
        igraph_vit_destroy(&tovit);
        igraph_vector_destroy(&indexv);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    igraph_vit_destroy(&fromvit);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths_dijkstra
 * \brief Weighted shortest paths from a vertex.
 *
 * </para><para>
 * If there is more than one path with the smallest weight between two vertices, this
 * function gives only one of them.
 * \param graph The graph object.
 * \param vertices The result, the ids of the vertices along the paths.
 *        This is a pointer vector, each element points to a vector
 *        object. These should be initialized before passing them to
 *        the function, which will properly clear and/or resize them
 *        and fill the ids of the vertices along the geodesics from/to
 *        the vertices. Supply a null pointer here if you don't need
 *        these vectors. Normally, either this argument, or the \c
 *        edges should be non-null, but no error or warning is given
 *        if they are both null pointers.
 * \param edges The result, the ids of the edges along the paths.
 *        This is a pointer vector, each element points to a vector
 *        object. These should be initialized before passing them to
 *        the function, which will properly clear and/or resize them
 *        and fill the ids of the vertices along the geodesics from/to
 *        the vertices. Supply a null pointer here if you don't need
 *        these vectors. Normally, either this argument, or the \c
 *        vertices should be non-null, but no error or warning is given
 *        if they are both null pointers.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the ids of the vertices to/from which the
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
 * \param predecessors A pointer to an initialized igraph vector or null.
 *        If not null, a vector containing the predecessor of each vertex in
 *        the single source shortest path tree is returned here. The
 *        predecessor of vertex i in the tree is the vertex from which vertex i
 *        was reached. The predecessor of the start vertex (in the \c from
 *        argument) is itself by definition. If the predecessor is -1, it means
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
 *           \p from is invalid vertex id, or the length of \p to is
 *           not the same as the length of \p res.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(|E|log|E|+|V|), where |V| is the number of
 * vertices and |E| is the number of edges
 *
 * \sa \ref igraph_shortest_paths_dijkstra() if you only need the path length but
 * not the paths themselves, \ref igraph_get_shortest_paths() if all edge
 * weights are equal.
 *
 * \example examples/simple/igraph_get_shortest_paths_dijkstra.c
 */
int igraph_get_shortest_paths_dijkstra(const igraph_t *graph,
                                       igraph_vector_ptr_t *vertices,
                                       igraph_vector_ptr_t *edges,
                                       igraph_integer_t from,
                                       igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode,
                                       igraph_vector_long_t *predecessors,
                                       igraph_vector_long_t *inbound_edges) {
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
         computation, as IGRAPH_FINITE() might involve a function call
         and we want to spare that. So we store distance+1.0 instead of
         distance, and zero denotes infinity.
       - `parents' assigns the inbound edge IDs of all vertices in the
         shortest path tree to the vertices. In this implementation, the
         edge ID + 1 is stored, zero means unreachable vertices.
    */

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t dists;
    long int *parents;
    igraph_bool_t *is_target;
    long int i, to_reach;

    if (!weights) {
        return igraph_get_shortest_paths(graph, vertices, edges, from, to, mode,
                                         predecessors, inbound_edges);
    }

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

    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    if (vertices && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(vertices)) {
        IGRAPH_ERROR("Size of `vertices' and `to' should match", IGRAPH_EINVAL);
    }
    if (edges && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(edges)) {
        IGRAPH_ERROR("Size of `edges' and `to' should match", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
    igraph_vector_fill(&dists, -1.0);

    parents = IGRAPH_CALLOC(no_of_nodes, long int);
    if (parents == 0) {
        IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, parents);
    is_target = IGRAPH_CALLOC(no_of_nodes, igraph_bool_t);
    if (is_target == 0) {
        IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, is_target);

    /* Mark the vertices we need to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (!is_target[ (long int) IGRAPH_VIT_GET(vit) ]) {
            is_target[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    VECTOR(dists)[(long int)from] = 0.0;  /* zero distance */
    parents[(long int)from] = 0;
    igraph_2wheap_push_with_index(&Q, from, 0);

    while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
        long int nlen, minnei = igraph_2wheap_max_index(&Q);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        if (is_target[minnei]) {
            is_target[minnei] = 0;
            to_reach--;
        }

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_lazy_inclist_get(&inclist, (igraph_integer_t) minnei);
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            long int edge = (long int) VECTOR(*neis)[i];
            long int tto = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
            igraph_real_t curdist = VECTOR(dists)[tto];
            if (curdist < 0) {
                /* This is the first finite distance */
                VECTOR(dists)[tto] = altdist;
                parents[tto] = edge + 1;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
            } else if (altdist < curdist) {
                /* This is a shorter path */
                VECTOR(dists)[tto] = altdist;
                parents[tto] = edge + 1;
                IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, -altdist));
            }
        }
    } /* !igraph_2wheap_empty(&Q) */

    if (to_reach > 0) {
        IGRAPH_WARNING("Couldn't reach some vertices");
    }

    /* Create `predecessors' if needed */
    if (predecessors) {
        IGRAPH_CHECK(igraph_vector_long_resize(predecessors, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (i == from) {
                /* i is the start vertex */
                VECTOR(*predecessors)[i] = i;
            } else if (parents[i] <= 0) {
                /* i was not reached */
                VECTOR(*predecessors)[i] = -1;
            } else {
                /* i was reached via the edge with ID = parents[i] - 1 */
                VECTOR(*predecessors)[i] = IGRAPH_OTHER(graph, parents[i] - 1, i);
            }
        }
    }

    /* Create `inbound_edges' if needed */
    if (inbound_edges) {
        IGRAPH_CHECK(igraph_vector_long_resize(inbound_edges, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (parents[i] <= 0) {
                /* i was not reached */
                VECTOR(*inbound_edges)[i] = -1;
            } else {
                /* i was reached via the edge with ID = parents[i] - 1 */
                VECTOR(*inbound_edges)[i] = parents[i] - 1;
            }
        }
    }

    /* Reconstruct the shortest paths based on vertex and/or edge IDs */
    if (vertices || edges) {
        for (IGRAPH_VIT_RESET(vit), i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
            long int node = IGRAPH_VIT_GET(vit);
            long int size, act, edge;
            igraph_vector_t *vvec = 0, *evec = 0;
            if (vertices) {
                vvec = VECTOR(*vertices)[i];
                igraph_vector_clear(vvec);
            }
            if (edges) {
                evec = VECTOR(*edges)[i];
                igraph_vector_clear(evec);
            }

            IGRAPH_ALLOW_INTERRUPTION();

            size = 0;
            act = node;
            while (parents[act]) {
                size++;
                edge = parents[act] - 1;
                act = IGRAPH_OTHER(graph, edge, act);
            }
            if (vvec && (size > 0 || node == from)) {
                IGRAPH_CHECK(igraph_vector_resize(vvec, size + 1));
                VECTOR(*vvec)[size] = node;
            }
            if (evec) {
                IGRAPH_CHECK(igraph_vector_resize(evec, size));
            }
            act = node;
            while (parents[act]) {
                edge = parents[act] - 1;
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
    IGRAPH_FREE(parents);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(6);

    return 0;
}

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
 *        pointer. If not a null pointer, then the vertex ids along
 *        the path are stored here, including the source and target
 *        vertices.
 * \param edges Pointer to an uninitialized vector or a null
 *        pointer. If not a null pointer, then the edge ids along the
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
 * Time complexity: O(|E|log|E|+|V|), |V| is the number of vertices,
 * |E| is the number of edges in the graph.
 *
 * \sa \ref igraph_get_shortest_paths_dijkstra() for the version with
 * more target vertices.
 */

int igraph_get_shortest_path_dijkstra(const igraph_t *graph,
                                      igraph_vector_t *vertices,
                                      igraph_vector_t *edges,
                                      igraph_integer_t from,
                                      igraph_integer_t to,
                                      const igraph_vector_t *weights,
                                      igraph_neimode_t mode) {

    igraph_vector_ptr_t vertices2, *vp = &vertices2;
    igraph_vector_ptr_t edges2, *ep = &edges2;

    if (vertices) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&vertices2, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &vertices2);
        VECTOR(vertices2)[0] = vertices;
    } else {
        vp = 0;
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&edges2, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &edges2);
        VECTOR(edges2)[0] = edges;
    } else {
        ep = 0;
    }

    IGRAPH_CHECK(igraph_get_shortest_paths_dijkstra(graph, vp, ep,
                 from, igraph_vss_1(to),
                 weights, mode, 0, 0));

    if (edges) {
        igraph_vector_ptr_destroy(&edges2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (vertices) {
        igraph_vector_ptr_destroy(&vertices2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/* Compares two paths based on their last elements. Required by
 * igraph_get_all_shortest_paths_dijkstra to put the final result
 * in order. Assumes that both paths are pointers to igraph_vector_t
 * objects and that they are not empty
 */
static int igraph_i_vector_tail_cmp(const void* path1, const void* path2) {
    return (int) (igraph_vector_tail(*(const igraph_vector_t**)path1) -
                  igraph_vector_tail(*(const igraph_vector_t**)path2));
}

/**
 * \ingroup structural
 * \function igraph_get_all_shortest_paths_dijkstra
 * \brief All weighted shortest paths (geodesics) from a vertex.
 *
 * \param graph The graph object.
 * \param res Pointer to an initialized pointer vector, the result
 *   will be stored here in igraph_vector_t objects. Each vector
 *   object contains the vertices along a shortest path from \p from
 *   to another vertex. The vectors are ordered according to their
 *   target vertex: first the shortest paths to vertex 0, then to
 *   vertex 1, etc. No data is included for unreachable vertices.
 * \param nrgeo Pointer to an initialized igraph_vector_t object or
 *   NULL. If not NULL the number of shortest paths from \p from are
 *   stored here for every vertex in the graph. Note that the values
 *   will be accurate only for those vertices that are in the target
 *   vertex sequence (see \p to), since the search terminates as soon
 *   as all the target vertices have been found.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the ids of the vertices to/from which the
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
 *           \p from is invalid vertex id, or the length of \p to is
 *           not the same as the length of \p res.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(|E|log|E|+|V|), where |V| is the number of
 * vertices and |E| is the number of edges
 *
 * \sa \ref igraph_shortest_paths_dijkstra() if you only need the path
 * length but not the paths themselves, \ref igraph_get_all_shortest_paths()
 * if all edge weights are equal.
 *
 * \example examples/simple/igraph_get_all_shortest_paths_dijkstra.c
 */
int igraph_get_all_shortest_paths_dijkstra(const igraph_t *graph,
        igraph_vector_ptr_t *res,
        igraph_vector_t *nrgeo,
        igraph_integer_t from, igraph_vs_t to,
        const igraph_vector_t *weights,
        igraph_neimode_t mode) {
    /* Implementation details: see igraph_get_shortest_paths_dijkstra,
       it's basically the same.
    */

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t dists, order;
    igraph_vector_ptr_t parents;
    igraph_finally_func_t *res_item_destructor;
    unsigned char *is_target;
    long int i, n, to_reach;

    if (!weights) {
        return igraph_get_all_shortest_paths(graph, res, nrgeo, from, to, mode);
    }

    if (res == 0 && nrgeo == 0) {
        return IGRAPH_SUCCESS;
    }

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

    /* parents stores a vector for each vertex, listing the parent vertices
     * of each vertex in the traversal */
    IGRAPH_CHECK(igraph_vector_ptr_init(&parents, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &parents);
    igraph_vector_ptr_set_item_destructor(&parents, (igraph_finally_func_t*)igraph_vector_destroy);
    for (i = 0; i < no_of_nodes; i++) {
        igraph_vector_t* parent_vec;
        parent_vec = IGRAPH_CALLOC(1, igraph_vector_t);
        if (parent_vec == 0) {
            IGRAPH_ERROR("cannot run igraph_get_all_shortest_paths", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_init(parent_vec, 0));
        VECTOR(parents)[i] = parent_vec;
    }

    /* distance of each vertex from the root */
    IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
    igraph_vector_fill(&dists, -1.0);

    /* order lists the order of vertices in which they were found during
     * the traversal */
    IGRAPH_VECTOR_INIT_FINALLY(&order, 0);

    /* boolean array to mark whether a given vertex is a target or not */
    is_target = IGRAPH_CALLOC(no_of_nodes, unsigned char);
    if (is_target == 0) {
        IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM);
    }
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
        if (!is_target[ (long int) IGRAPH_VIT_GET(vit) ]) {
            is_target[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    VECTOR(dists)[(long int)from] = 0.0;  /* zero distance */
    igraph_2wheap_push_with_index(&Q, from, 0);

    while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
        long int nlen, minnei = igraph_2wheap_max_index(&Q);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        /*
        printf("Reached vertex %ld, is_target[%ld] = %d, %ld to go\n",
            minnei, minnei, (int)is_target[minnei], to_reach - is_target[minnei]);
        */

        if (is_target[minnei]) {
            is_target[minnei] = 0;
            to_reach--;
        }

        /* Mark that we have reached this vertex */
        IGRAPH_CHECK(igraph_vector_push_back(&order, minnei));

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_lazy_inclist_get(&inclist, (igraph_integer_t) minnei);
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            long int edge = (long int) VECTOR(*neis)[i];
            long int tto = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
            igraph_real_t curdist = VECTOR(dists)[tto];
            igraph_vector_t *parent_vec;

            if (curdist < 0) {
                /* This is the first non-infinite distance */
                VECTOR(dists)[tto] = altdist;
                parent_vec = (igraph_vector_t*)VECTOR(parents)[tto];
                IGRAPH_CHECK(igraph_vector_push_back(parent_vec, minnei));
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
            } else if (altdist == curdist && VECTOR(*weights)[edge] > 0) {
                /* This is an alternative path with exactly the same length.
                     * Note that we consider this case only if the edge via which we
                     * reached the node has a nonzero weight; otherwise we could create
                     * infinite loops in undirected graphs by traversing zero-weight edges
                     * back-and-forth */
                parent_vec = (igraph_vector_t*)VECTOR(parents)[tto];
                IGRAPH_CHECK(igraph_vector_push_back(parent_vec, minnei));
            } else if (altdist < curdist) {
                /* This is a shorter path */
                VECTOR(dists)[tto] = altdist;
                parent_vec = (igraph_vector_t*)VECTOR(parents)[tto];
                igraph_vector_clear(parent_vec);
                IGRAPH_CHECK(igraph_vector_push_back(parent_vec, minnei));
                IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, -altdist));
            }
        }
    } /* !igraph_2wheap_empty(&Q) */

    if (to_reach > 0) {
        IGRAPH_WARNING("Couldn't reach some vertices");
    }

    /* we don't need these anymore */
    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    /*
    printf("Order:\n");
    igraph_vector_print(&order);

    printf("Parent vertices:\n");
    for (i = 0; i < no_of_nodes; i++) {
      if (igraph_vector_size(VECTOR(parents)[i]) > 0) {
        printf("[%ld]: ", (long int)i);
        igraph_vector_print(VECTOR(parents)[i]);
      }
    }
    */

    if (nrgeo) {
        IGRAPH_CHECK(igraph_vector_resize(nrgeo, no_of_nodes));
        igraph_vector_null(nrgeo);

        /* Theoretically, we could calculate nrgeo in parallel with the traversal.
         * However, that way we would have to check whether nrgeo is null or not
         * every time we want to update some element in nrgeo. Since we need the
         * order vector anyway for building the final result, we could just as well
         * build nrgeo here.
         */
        VECTOR(*nrgeo)[(long int)from] = 1;
        n = igraph_vector_size(&order);
        for (i = 1; i < n; i++) {
            long int node, j, k;
            igraph_vector_t *parent_vec;

            node = (long int)VECTOR(order)[i];
            /* now, take the parent vertices */
            parent_vec = (igraph_vector_t*)VECTOR(parents)[node];
            k = igraph_vector_size(parent_vec);
            for (j = 0; j < k; j++) {
                VECTOR(*nrgeo)[node] += VECTOR(*nrgeo)[(long int)VECTOR(*parent_vec)[j]];
            }
        }
    }

    if (res) {
        igraph_vector_t *path, *paths_index, *parent_vec;
        igraph_stack_t stack;
        long int j, node;

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

            IGRAPH_CHECK(igraph_stack_init(&stack, 0));
            IGRAPH_FINALLY(igraph_stack_destroy, &stack);

            /* Add the target vertices to the queue */
            IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
            IGRAPH_FINALLY(igraph_vit_destroy, &vit);
            for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
                i = (long int) IGRAPH_VIT_GET(vit);
                if (!is_target[i]) {
                    is_target[i] = 1;
                    IGRAPH_CHECK(igraph_stack_push(&stack, i));
                }
            }
            igraph_vit_destroy(&vit);
            IGRAPH_FINALLY_CLEAN(1);

            while (!igraph_stack_empty(&stack)) {
                /* For each parent of node i, get its parents */
                igraph_real_t el = igraph_stack_pop(&stack);
                parent_vec = (igraph_vector_t*)VECTOR(parents)[(long int) el];
                i = igraph_vector_size(parent_vec);

                for (j = 0; j < i; j++) {
                    /* For each parent, check if it's already in the stack.
                     * If not, push it and mark it in is_target */
                    n = (long int) VECTOR(*parent_vec)[j];
                    if (!is_target[n]) {
                        is_target[n] = 2;
                        IGRAPH_CHECK(igraph_stack_push(&stack, n));
                    }
                }
            }
            igraph_stack_destroy(&stack);
            IGRAPH_FINALLY_CLEAN(1);
        }

        /* now, reconstruct the shortest paths from the parent list in the
         * order we've found the nodes during the traversal.
         * dists is being re-used as a vector where element i tells the
         * index in res where the shortest paths leading to vertex i
         * start, plus one (so that zero means that there are no paths
         * for a given vertex).
         */
        paths_index = &dists;
        n = igraph_vector_size(&order);
        igraph_vector_null(paths_index);

        /* clear the paths vector */
        igraph_vector_ptr_clear(res);
        res_item_destructor = igraph_vector_ptr_get_item_destructor(res);
        igraph_vector_ptr_set_item_destructor(res,
                                              (igraph_finally_func_t*)igraph_vector_destroy);

        /* by definition, the shortest path leading to the starting vertex
         * consists of the vertex itself only */
        path = IGRAPH_CALLOC(1, igraph_vector_t);
        if (path == 0)
            IGRAPH_ERROR("cannot run igraph_get_all_shortest_paths_dijkstra",
                         IGRAPH_ENOMEM);
        IGRAPH_FINALLY(igraph_free, path);
        IGRAPH_CHECK(igraph_vector_init(path, 1));
        IGRAPH_CHECK(igraph_vector_ptr_push_back(res, path));
        IGRAPH_FINALLY_CLEAN(1);  /* ownership of path passed to res */
        VECTOR(*path)[0] = from;
        VECTOR(*paths_index)[(long int)from] = 1;

        for (i = 1; i < n; i++) {
            long int m, path_count;
            igraph_vector_t *parent_path;

            node = (long int) VECTOR(order)[i];

            /* if we don't need the shortest paths for this node (because
             * it is not standing in a shortest path between the source
             * node and any of the target nodes), skip it */
            if (!is_target[node]) {
                continue;
            }

            IGRAPH_ALLOW_INTERRUPTION();

            /* we are calculating the shortest paths of node now. */
            /* first, we update the paths_index */
            path_count = igraph_vector_ptr_size(res);
            VECTOR(*paths_index)[node] = path_count + 1;
            /* res_end = (igraph_vector_t*)&(VECTOR(*res)[path_count]); */

            /* now, take the parent vertices */
            parent_vec = (igraph_vector_t*)VECTOR(parents)[node];
            m = igraph_vector_size(parent_vec);

            /*
            printf("Calculating shortest paths to vertex %ld\n", node);
            printf("Parents are: ");
            igraph_vector_print(parent_vec);
            */

            for (j = 0; j < m; j++) {
                /* for each parent, copy the shortest paths leading to that parent
                 * and add the current vertex in the end */
                long int parent_node = (long int) VECTOR(*parent_vec)[j];
                long int parent_path_idx = (long int) VECTOR(*paths_index)[parent_node] - 1;
                /*
                printf("  Considering parent: %ld\n", parent_node);
                printf("  Paths to parent start at index %ld in res\n", parent_path_idx);
                */
                IGRAPH_ASSERT(parent_path_idx >= 0);
                for (; parent_path_idx < path_count; parent_path_idx++) {
                    parent_path = (igraph_vector_t*)VECTOR(*res)[parent_path_idx];
                    if (igraph_vector_tail(parent_path) != parent_node) {
                        break;
                    }

                    path = IGRAPH_CALLOC(1, igraph_vector_t);
                    if (path == 0)
                        IGRAPH_ERROR("cannot run igraph_get_all_shortest_paths_dijkstra",
                                     IGRAPH_ENOMEM);
                    IGRAPH_FINALLY(igraph_free, path);
                    IGRAPH_CHECK(igraph_vector_copy(path, parent_path));
                    IGRAPH_CHECK(igraph_vector_ptr_push_back(res, path));
                    IGRAPH_FINALLY_CLEAN(1);  /* ownership of path passed to res */
                    IGRAPH_CHECK(igraph_vector_push_back(path, node));
                }
            }
        }

        /* remove the path vector's original item destructor */
        igraph_vector_ptr_set_item_destructor(res, res_item_destructor);

        /* free those paths from the result vector which we won't need */
        n = igraph_vector_ptr_size(res);
        j = 0;
        for (i = 0; i < n; i++) {
            igraph_real_t tmp;
            path = (igraph_vector_t*)VECTOR(*res)[i];
            tmp = igraph_vector_tail(path);
            if (is_target[(long int)tmp] == 1) {
                /* we need this path, keep it */
                VECTOR(*res)[j] = path;
                j++;
            } else {
                /* we don't need this path, free it */
                igraph_vector_destroy(path); free(path);
            }
        }
        IGRAPH_CHECK(igraph_vector_ptr_resize(res, j));

        /* sort the paths by the target vertices */
        igraph_vector_ptr_sort(res, igraph_i_vector_tail_cmp);
    }

    /* free the allocated memory */
    igraph_vector_destroy(&order);
    IGRAPH_FREE(is_target);
    igraph_vector_destroy(&dists);
    igraph_vector_ptr_destroy_all(&parents);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}
