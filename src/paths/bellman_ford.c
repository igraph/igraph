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
#include "igraph_stack.h"

#include "core/indheap.h"
#include "core/interruption.h"

#include <string.h>

/**
 * \function igraph_shortest_paths_bellman_ford
 * \brief Weighted shortest path lengths between vertices, allowing negative weights.
 *
 * This function implements the Bellman-Ford algorithm to find the weighted
 * shortest paths to all vertices from a single source, allowing negative weights.
 * It is run independently for the given sources. If there are no negative
 * weights, you are better off with \ref igraph_shortest_paths_dijkstra() .
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here, the matrix will be resized if needed.
 *    Each row contains the distances from a single source, to all
 *    vertices in the graph, in the order of vertex ids. For unreachable
 *    vertices the matrix contains \c IGRAPH_INFINITY.
 * \param from The source vertices.
 * \param to The target vertices. It is not allowed to include a
 *    vertex twice or more.
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
int igraph_shortest_paths_bellman_ford(const igraph_t *graph,
                                       igraph_matrix_t *res,
                                       const igraph_vs_t from,
                                       const igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_lazy_inclist_t inclist;
    long int i, j, k;
    long int no_of_from, no_of_to;
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
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
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
        long int source = IGRAPH_VIT_GET(fromvit);

        igraph_vector_fill(&dist, my_infinity);
        VECTOR(dist)[source] = 0;
        igraph_vector_null(&clean_vertices);
        igraph_vector_null(&num_queued);

        /* Fill the queue with vertices to be checked */
        for (j = 0; j < no_of_nodes; j++) {
            IGRAPH_CHECK(igraph_dqueue_push(&Q, j));
        }

        while (!igraph_dqueue_empty(&Q)) {
            igraph_vector_int_t *neis;
            long int nlen;

            j = (long int) igraph_dqueue_pop(&Q);
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

            neis = igraph_lazy_inclist_get(&inclist, (igraph_integer_t) j);
            nlen = igraph_vector_int_size(neis);

            for (k = 0; k < nlen; k++) {
                long int nei = (long int) VECTOR(*neis)[k];
                long int target = IGRAPH_OTHER(graph, nei, j);
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
                long int v = IGRAPH_VIT_GET(tovit);
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


/**
 * \ingroup structural
 * \function igraph_get_shortest_paths_bellman_ford
 * \brief Weighted shortest paths from a vertex, allowing negative weights.
 *
 * This function calculates weighted shortest paths from or to a single vertex,
 * and allows negative weights. When there is more than one shortest path between
 * two vertices, only one of them is returned.
 *
 * If there are no negative weights, you are better off with 
 * \ref igraph_get_shortest_paths_dijkstra() .
 *
 * \param graph The input graph, can be directed.
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
 * \param weights The edge weights. There mustn't be any closed loop in
 *    the graph that has a negative total weight (since this would allow
 *    us to decrease the weight of any path containing at least a single
 *    vertex of this loop infinitely). If this is a null pointer, then the
 *    unweighted version, \ref igraph_shortest_paths() is called.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
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
 *         \clist
 *         \cli IGRAPH_ENOMEM
 *           Not enough memory for temporary data.
 *         \cli IGRAPH_EINVAL
 *           The weight vector doesn't math the number of edges.
 *         \cli IGRAPH_EINVVID
 *           \p from is invalid vertex id, or the length of \p to is
 *           not the same as the length of \p vertices or \p edges.
 *         \cli IGRAPH_ENEGLOOP
 *           Bellman-ford algorithm encounted a negative loop.
 *         \endclist
 *
 * Time complexity: O(|E|*|V|), where |V| is the number of
 * vertices, |E| the number of edges.
 *
 * \sa \ref igraph_shortest_paths() for a faster unweighted version
 * or \ref igraph_shortest_paths_dijkstra() if you do not have negative
 * edge weights.
 */

int igraph_get_shortest_paths_bellman_ford(const igraph_t *graph,
                                        igraph_vector_ptr_t *vertices,
                                        igraph_vector_ptr_t *edges,
                                        igraph_integer_t from,
                                        igraph_vs_t to,
                                        const igraph_vector_t *weights,
                                        igraph_neimode_t mode,
                                        igraph_vector_long_t *predecessors,
                                        igraph_vector_long_t *inbound_edges) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int *parents;
    igraph_lazy_inclist_t inclist;
    long int i, j, k;
    igraph_dqueue_t Q;
    igraph_vector_t clean_vertices;
    igraph_vector_t num_queued;
    igraph_vit_t tovit;
    igraph_real_t my_infinity = IGRAPH_INFINITY;
    igraph_vector_t dist;

    if (!weights) {
        return  igraph_get_shortest_paths(graph, vertices, edges, from, to, mode,
                                         predecessors, inbound_edges);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length must match number of edges.", IGRAPH_EINVAL);
    }

    IGRAPH_DQUEUE_INIT_FINALLY(&Q, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&clean_vertices, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&num_queued, no_of_nodes);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
    IGRAPH_FINALLY(igraph_vit_destroy, &tovit);

    if (vertices && IGRAPH_VIT_SIZE(tovit) != igraph_vector_ptr_size(vertices)) {
        IGRAPH_ERROR("Size of `vertices' and `to' should match.", IGRAPH_EINVAL);
    }
    if (edges && IGRAPH_VIT_SIZE(tovit) != igraph_vector_ptr_size(edges)) {
        IGRAPH_ERROR("Size of `edges' and `to' should match.", IGRAPH_EINVAL);
    }

    parents = IGRAPH_CALLOC(no_of_nodes, long int);
    if (parents == 0) {
        IGRAPH_ERROR("Insufficient memory for shortest paths with Bellman-Ford.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, parents);
    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    igraph_vector_fill(&dist, my_infinity);
    VECTOR(dist)[from] = 0;
    igraph_vector_null(&clean_vertices);
    igraph_vector_null(&num_queued);

    /* Fill the queue with vertices to be checked */
    for (j = 0; j < no_of_nodes; j++) {
        IGRAPH_CHECK(igraph_dqueue_push(&Q, j));
    }

    while (!igraph_dqueue_empty(&Q)) {
        igraph_vector_int_t *neis;
        long int nlen;

        j = (long int) igraph_dqueue_pop(&Q);
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

        neis = igraph_lazy_inclist_get(&inclist, (igraph_integer_t) j);
        nlen = igraph_vector_int_size(neis);

        for (k = 0; k < nlen; k++) {
            long int nei = (long int) VECTOR(*neis)[k];
            long int target = IGRAPH_OTHER(graph, nei, j);
            if (VECTOR(dist)[target] > VECTOR(dist)[j] + VECTOR(*weights)[nei]) {
                /* relax the edge */
                VECTOR(dist)[target] = VECTOR(dist)[j] + VECTOR(*weights)[nei];
                parents[target] = nei + 1;
                if (VECTOR(clean_vertices)[target]) {
                    VECTOR(clean_vertices)[target] = 0;
                    IGRAPH_CHECK(igraph_dqueue_push(&Q, target));
                }
            }
        }
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
        for (IGRAPH_VIT_RESET(tovit), i = 0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit), i++) {
            long int node = IGRAPH_VIT_GET(tovit);
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
            if (vvec) {
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

    igraph_vector_destroy(&dist);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vit_destroy(&tovit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FREE(parents);
    igraph_dqueue_destroy(&Q);
    igraph_vector_destroy(&clean_vertices);
    igraph_vector_destroy(&num_queued);
    igraph_lazy_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_shortest_path_bellman_ford
 * \brief Weighted shortest path from one vertex to another one.
 *
 * Calculates a single (positively) weighted shortest path from
 * a single vertex to another one, using Bellman-Ford algorithm.
 *
 * </para><para>
 * This function is a special case (and a wrapper) to
 * \ref igraph_get_shortest_paths_bellman_ford().
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
 * \param weights The edge weights. There mustn't be any closed loop in
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

int igraph_get_shortest_path_bellman_ford(const igraph_t *graph,
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
        vp = NULL;
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&edges2, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &edges2);
        VECTOR(edges2)[0] = edges;
    } else {
        ep = NULL;
    }

    IGRAPH_CHECK(igraph_get_shortest_paths_bellman_ford(graph, vp, ep,
                                                        from, igraph_vss_1(to),
                                                        weights, mode, NULL, NULL));

    if (edges) {
        igraph_vector_ptr_destroy(&edges2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (vertices) {
        igraph_vector_ptr_destroy(&vertices2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
