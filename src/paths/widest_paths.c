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

#include "core/indheap.h"
#include "core/interruption.h"
#include "internal/utils.h"

/**
 * \function igraph_get_widest_paths
 * \brief Widest paths from a single vertex.
 *
 * Calculates the widest paths from a single node to all other specified nodes,
 * using a modified Dijkstra's algorithm.  If there is more than one path with
 * the largest width between two vertices, this function gives only one of them.
 * \param graph The graph object.
 * \param vertices The result, the IDs of the vertices along the paths.
 *        This is a list of integer vectors where each element is an
 *        \ref igraph_vector_int_t object. The list will be resized as needed.
 *        Supply a null pointer here if you don't need these vectors.
 * \param edges The result, the IDs of the edges along the paths.
 *        This is a list of integer vectors where each element is an
 *        \ref igraph_vector_int_t object. The list will be resized as needed.
 *        Supply a null pointer here if you don't need these vectors.
 * \param from The id of the vertex from/to which the widest paths are
 *        calculated.
 * \param to Vertex sequence with the IDs of the vertices to/from which the
 *        widest paths will be calculated. A vertex might be given multiple
 *        times.
 * \param weights The edge weights. Edge weights can be negative. If this
 *        is a null pointer or if any edge weight is NaN, then an error
 *        is returned. Edges with positive infinite weight are ignored.
 * \param mode The type of widest paths to be used for the
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
 *        the single source widest path tree is returned here. The
 *        parent of vertex i in the tree is the vertex from which vertex i
 *        was reached. The parent of the start vertex (in the \c from
 *        argument) is -1. If the parent is -2, it means
 *        that the given vertex was not reached from the source during the
 *        search. Note that the search terminates if all the vertices in
 *        \c to are reached.
 * \param inbound_edges A pointer to an initialized igraph vector or null.
 *        If not null, a vector containing the inbound edge of each vertex in
 *        the single source widest path tree is returned here. The
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
 * Time complexity: O(|E|log|E|+|V|), where |V| is the number of
 * vertices in the graph and |E| is the number of edges
 *
 * \sa \ref igraph_widest_path_widths_dijkstra() or
 * \ref igraph_widest_path_widths_floyd_warshall() if you only need the
 * widths of the paths but not the paths themselves.
 */
igraph_error_t igraph_get_widest_paths(const igraph_t *graph,
                                       igraph_vector_int_list_t *vertices,
                                       igraph_vector_int_list_t *edges,
                                       igraph_integer_t from,
                                       igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode,
                                       igraph_vector_int_t *parents,
                                       igraph_vector_int_t *inbound_edges) {

    /* Implementation details: This is a Dijkstra algorithm with a
    binary heap, modified to support widest paths. The heap is indexed,
    so it stores both the widest path to a node, as well as it's index. We
    use a 2 way heap so that we can query indexes directly in the heap.

    To adapt a Dijkstra to handle widest path, instead of prioritising candidate
    nodes with the minimum distance, we prioritise those with the maximum
    width instead. When adding a node into our set of 'completed' nodes, we
    update all neighbouring nodes with a width that is equal to the min of the
    width to the current node and the width of the edge.

    We denote the widest path from a node to itself as infinity, and the widest
    path from a node to a node it cannot reach as negative infinity.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t widths;
    igraph_integer_t *parent_eids;
    bool *is_target;
    igraph_integer_t i, to_reach;

    if (!weights) {
        IGRAPH_ERROR("Weight vector is required.", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
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

    IGRAPH_VECTOR_INIT_FINALLY(&widths, no_of_nodes);
    igraph_vector_fill(&widths, IGRAPH_NEGINFINITY);

    parent_eids = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(parent_eids, "Insufficient memory for widest paths.");
    IGRAPH_FINALLY(igraph_free, parent_eids);

    is_target = IGRAPH_CALLOC(no_of_nodes, bool);
    IGRAPH_CHECK_OOM(is_target, "Insufficient memory for widest paths.");
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

    VECTOR(widths)[from] = IGRAPH_POSINFINITY;
    parent_eids[from] = 0;
    igraph_2wheap_push_with_index(&Q, from, IGRAPH_POSINFINITY);

    while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
        igraph_integer_t nlen, maxnei = igraph_2wheap_max_index(&Q);
        igraph_real_t maxwidth = igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        if (is_target[maxnei]) {
            is_target[maxnei] = false;
            to_reach--;
        }

        /* Now check all neighbors of 'maxnei' for a wider path */
        neis = igraph_lazy_inclist_get(&inclist, maxnei);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, maxnei);
            igraph_real_t edgewidth = VECTOR(*weights)[edge];
            igraph_real_t altwidth = maxwidth < edgewidth ? maxwidth : edgewidth;
            igraph_real_t curwidth = VECTOR(widths)[tto];
            if (edgewidth == IGRAPH_INFINITY) {
                /* Ignore edges with infinite weight */
            } else if (curwidth < 0) {
                /* This is the first assigning a width to this vertex */
                VECTOR(widths)[tto] = altwidth;
                parent_eids[tto] = edge + 1;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, altwidth));
            } else if (altwidth > curwidth) {
                /* This is a wider path */
                VECTOR(widths)[tto] = altwidth;
                parent_eids[tto] = edge + 1;
                igraph_2wheap_modify(&Q, tto, altwidth);
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
    /* Reconstruct the widest paths based on vertex and/or edge IDs */
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
    igraph_vector_destroy(&widths);
    IGRAPH_FREE(is_target);
    IGRAPH_FREE(parent_eids);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_widest_path
 * \brief Widest path from one vertex to another one.
 *
 * Calculates a single widest path from a single vertex to another
 * one, using Dijkstra's algorithm.
 *
 * </para><para>This function is a special case (and a wrapper) to
 * \ref igraph_get_widest_paths().
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
 * \param weights The edge weights. Edge weights can be negative. If this
 *        is a null pointer or if any edge weight is NaN, then an error
 *        is returned. Edges with positive infinite weight are ignored.
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
 * \sa \ref igraph_get_widest_paths() for the version with
 * more target vertices.
 */
igraph_error_t igraph_get_widest_path(const igraph_t *graph,
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

    IGRAPH_CHECK(igraph_get_widest_paths(graph, vp, ep,
                 from, igraph_vss_1(to),
                 weights, mode, 0, 0));

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

/**
 * \function igraph_widest_path_widths_floyd_warshall
 * \brief Widths of widest paths between vertices.
 *
 * This function implements a modified Floyd-Warshall algorithm,
 * to find the widest path widths between a set of source and target
 * vertices. It is primarily useful for all-pairs path widths in very dense
 * graphs, as its running time is manily determined by the vertex count,
 * and is not sensitive to the graph density. In sparse graphs, other methods
 * such as the Dijkstra algorithm, will perform better.
 *
 * </para><para>
 * Note that internally this function always computes the path width matrix
 * for all pairs of vertices. The \p from and \p to parameters only serve
 * to subset this matrix, but do not affect the time taken by the
 * calculation.
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here. The matrix will be resized as needed.
 *    Each row contains the widths from a single source, to the
 *    vertices given in the \c to argument.
 *    Unreachable vertices have width \c IGRAPH_NEGINFINITY, and vertices
 *    have a width of \c IGRAPH_POSINFINITY to themselves.
 * \param from The source vertices.
 * \param to The target vertices.
 * \param weights The edge weights. Edge weights can be negative. If this
 *        is a null pointer or if any edge weight is NaN, then an error
 *        is returned. Edges with positive infinite weight are ignored.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V|^3), where |V| is the number of vertices in the graph.
 *
 * \sa \ref igraph_widest_path_widths_dijkstra() for a variant that runs faster
 * on sparse graphs.
 */
igraph_error_t igraph_widest_path_widths_floyd_warshall(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    /* Implementation Details: This is a modified Floyd Warshall algorithm
    which computes the widest path between every pair of nodes. The key
    difference between this and the regular Floyd Warshall is that instead
    of updating the distance between two nodes to be the minimum of itself
    and the distance through an intermediate node, we instead set the width
    to be the maximum of itself and the width through the intermediate node.

    We denote the widest path from a node to itself as infinity, and the widest
    path from a node to a node it cannot reach as negative infinity.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t in = false, out = false;

    if (! weights) {
        IGRAPH_ERROR("Weight vector is required.", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    switch (mode) {
    case IGRAPH_ALL:
        in = out = true;
        break;
    case IGRAPH_OUT:
        out = true;
        break;
    case IGRAPH_IN:
        in = true;
        break;
    default:
        IGRAPH_ERROR("Invalid mode for Floyd-Warshall shortest path calculation.", IGRAPH_EINVMODE);
    }

    /* Fill out adjacency matrix */
    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_fill(res, IGRAPH_NEGINFINITY);
    for (igraph_integer_t i=0; i < no_of_nodes; i++) {
        MATRIX(*res, i, i) = IGRAPH_POSINFINITY;
    }

    for (igraph_integer_t edge=0; edge < no_of_edges; edge++) {
        igraph_integer_t from = IGRAPH_FROM(graph, edge);
        igraph_integer_t to = IGRAPH_TO(graph, edge);
        igraph_real_t w = VECTOR(*weights)[edge];

        if (w == IGRAPH_INFINITY) {
            /* Ignore edges with infinite weight */
            continue;
        }

        if (out && MATRIX(*res, from, to) < w) MATRIX(*res, from, to) = w;
        if (in  && MATRIX(*res, to, from) < w) MATRIX(*res, to, from) = w;
    }

    /* Run modified Floyd Warshall */
    for (igraph_integer_t k = 0; k < no_of_nodes; k++) {
        /* Iterate in column-major order for better performance */
        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            igraph_real_t width_kj = MATRIX(*res, k, j);
            if (j == k || width_kj == IGRAPH_NEGINFINITY) continue;

            IGRAPH_ALLOW_INTERRUPTION();

            for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
                if (i == j || i == k) continue;

                /* alternative_width := min(A(i,k), A(k,j))
                   A(i,j) := max(A(i,j), alternative_width) */

                igraph_real_t altwidth = MATRIX(*res, i, k);
                if (width_kj < altwidth) {
                    altwidth = width_kj;
                }
                if (altwidth > MATRIX(*res, i, j)) {
                    MATRIX(*res, i, j) = altwidth;
                }
            }
        }
    }

    IGRAPH_CHECK(igraph_i_matrix_subset_vertices(res, graph, from, to));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_widest_path_widths_dijkstra
 * \brief Widths of widest paths between vertices.
 *
 * This function implements a modified Dijkstra's algorithm, which
 * can find the widest path widths from a source vertex to all
 * other vertices. This function allows specifying a set of source
 * and target vertices. The algorithm is run independently for each
 * source and the results are retained only for the specified targets.
 * This implementation uses a binary heap for efficiency.
 *
 * \param graph The input graph, can be directed.
 * \param res The result, a matrix. A pointer to an initialized matrix
 *    should be passed here. The matrix will be resized as needed.
 *    Each row contains the widths from a single source, to the
 *    vertices given in the \c to argument.
 *    Unreachable vertices have width \c IGRAPH_NEGINFINITY, and vertices
 *    have a width of \c IGRAPH_POSINFINITY to themselves.
 * \param from The source vertices.
 * \param to The target vertices. It is not allowed to include a
 *    vertex twice or more.
 * \param weights The edge weights. Edge weights can be negative. If this
 *        is a null pointer or if any edge weight is NaN, then an error
 *        is returned. Edges with positive infinite weight are ignored.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(s*(|E|log|E|+|V|)), where |V| is the number of
 * vertices in the graph, |E| the number of edges and s the number of sources.
 *
 * \sa \ref igraph_widest_path_widths_floyd_warshall() for a variant that runs faster
 * on dense graphs.
 */
igraph_error_t igraph_widest_path_widths_dijkstra(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    /* Implementation details: This is a Dijkstra algorithm with a
    binary heap, modified to support widest paths. The heap is indexed,
    so it stores both the widest path to a node, as well as it's index. We
    use a 2 way heap so that we can query indexes directly in the heap.

    To adapt a Dijkstra to handle widest path, instead of prioritising candidate
    nodes with the minimum distance, we prioritise those with the maximum
    width instead. When adding a node into our set of 'completed' nodes, we
    update all neighbouring nodes with a width that is equal to the min of the
    width to the current node and the width of the edge.

    We denote the widest path from a node to itself as infinity, and the widest
    path from a node to a node it cannot reach as negative infinity.
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
        IGRAPH_ERROR("Weight vector is required.", IGRAPH_EINVAL);
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
        for (i = 0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit)) {
            igraph_integer_t v = IGRAPH_VIT_GET(tovit);
            if (VECTOR(indexv)[v]) {
                IGRAPH_ERROR("Duplicate vertices in `to', this is not allowed.",
                             IGRAPH_EINVAL);
            }
            VECTOR(indexv)[v] = ++i;
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));
    igraph_matrix_fill(res, IGRAPH_NEGINFINITY);

    for (IGRAPH_VIT_RESET(fromvit), i = 0;
         !IGRAPH_VIT_END(fromvit);
         IGRAPH_VIT_NEXT(fromvit), i++) {

        igraph_integer_t reached = 0;
        igraph_integer_t source = IGRAPH_VIT_GET(fromvit);
        igraph_2wheap_clear(&Q);
        igraph_2wheap_push_with_index(&Q, source, IGRAPH_POSINFINITY);

        while (!igraph_2wheap_empty(&Q)) {
            igraph_integer_t maxnei = igraph_2wheap_max_index(&Q);
            igraph_real_t maxwidth = igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            igraph_integer_t nlen;

            IGRAPH_ALLOW_INTERRUPTION();

            if (all_to) {
                MATRIX(*res, i, maxnei) = maxwidth;
            } else {
                if (VECTOR(indexv)[maxnei]) {
                    MATRIX(*res, i, VECTOR(indexv)[maxnei] - 1) = maxwidth;
                    reached++;
                    if (reached == no_of_to) {
                        igraph_2wheap_clear(&Q);
                        break;
                    }
                }
            }

            /* Now check all neighbors of 'maxnei' for a wider path*/
            neis = igraph_lazy_inclist_get(&inclist, maxnei);
            IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_integer_t edge = VECTOR(*neis)[j];
                igraph_integer_t tto = IGRAPH_OTHER(graph, edge, maxnei);
                igraph_real_t edgewidth = VECTOR(*weights)[edge];
                igraph_real_t altwidth = maxwidth < edgewidth ? maxwidth : edgewidth;
                igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
                igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
                igraph_real_t curwidth = active ? igraph_2wheap_get(&Q, tto) : IGRAPH_POSINFINITY;
                if (edgewidth == IGRAPH_INFINITY) {
                    /* Ignore edges with infinite weight */
                } else if (!has) {
                    /* This is the first time assigning a width to this vertex */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, altwidth));
                } else if (altwidth > curwidth) {
                    /* This is a wider path */
                    igraph_2wheap_modify(&Q, tto, altwidth);
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
