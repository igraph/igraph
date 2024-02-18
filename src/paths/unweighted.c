/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team <igraph@igraph.org>

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
 * \ingroup structural
 * \function igraph_distances_cutoff
 * \brief Length of the shortest paths between vertices, with cutoff.
 *
 * \experimental
 *
 * This function is similar to \ref igraph_distances(), but
 * paths longer than \p cutoff will not be considered.
 *
 * \param graph The graph object.
 * \param res The result of the calculation, a matrix. A pointer to an
 *        initialized matrix, to be more precise. The matrix will be
 *        resized if needed. It will have the same
 *        number of rows as the length of the \p from
 *        argument, and its number of columns is the number of
 *        vertices in the \p to argument. One row of the matrix shows the
 *        distances from/to a given vertex to the ones in \p to.
 *        For the unreachable vertices \c IGRAPH_INFINITY is returned.
 * \param from The source vertices._d
 * \param to The target vertices. It is not allowed to include a
 *    vertex twice or more.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for
 *          the computation.
 *        \endclist
 * \param cutoff The maximal length of paths that will be considered.
 *    When the distance of two vertices is greater than this value,
 *    it will be returned as \c IGRAPH_INFINITY. Negative cutoffs are
 *    treated as infinity.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary
 *           data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(s |E| + |V|), where s is the number of source vertices to use,
 * and |V| and |E| are the number of vertices and edges in the graph.
 *
 * \sa  \ref igraph_distances_dijkstra_cutoff() for the weighted version with non-negative
 * weights.
 *
 * \example examples/simple/distances.c
 */
igraph_error_t igraph_distances_cutoff(const igraph_t *graph, igraph_matrix_t *res,
                          const igraph_vs_t from, const igraph_vs_t to,
                          igraph_neimode_t mode, igraph_real_t cutoff) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_from, no_of_to;
    igraph_integer_t *already_counted;
    igraph_adjlist_t adjlist;
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_bool_t all_to;

    igraph_integer_t i, j;
    igraph_vit_t fromvit, tovit;
    igraph_vector_int_t indexv;

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument.", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);
    no_of_from = IGRAPH_VIT_SIZE(fromvit);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    already_counted = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(already_counted, "Insufficient memory for graph distance calculation.");
    IGRAPH_FINALLY(igraph_free, already_counted);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

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
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, IGRAPH_VIT_GET(fromvit)));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        already_counted[ IGRAPH_VIT_GET(fromvit) ] = i + 1;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t act = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);

            if (cutoff >= 0 && actdist > cutoff) {
                continue;
            }

            if (all_to) {
                MATRIX(*res, i, act) = actdist;
            } else {
                if (VECTOR(indexv)[act]) {
                    MATRIX(*res, i, VECTOR(indexv)[act] - 1) = actdist;
                    reached++;
                    if (reached == no_of_to) {
                        igraph_dqueue_int_clear(&q);
                        break;
                    }
                }
            }

            neis = igraph_adjlist_get(&adjlist, act);
            igraph_integer_t nei_count = igraph_vector_int_size(neis);
            for (j = 0; j < nei_count; j++) {
                igraph_integer_t neighbor = VECTOR(*neis)[j];
                if (already_counted[neighbor] == i + 1) {
                    continue;
                }
                already_counted[neighbor] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
            }
        }
    }

    /* Clean */
    if (!all_to) {
        igraph_vit_destroy(&tovit);
        igraph_vector_int_destroy(&indexv);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_FREE(already_counted);
    igraph_dqueue_int_destroy(&q);
    igraph_vit_destroy(&fromvit);
    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_distances
 * \brief Length of the shortest paths between vertices.
 *
 * \param graph The graph object.
 * \param res The result of the calculation, a matrix. A pointer to an
 *        initialized matrix, to be more precise. The matrix will be
 *        resized if needed. It will have the same
 *        number of rows as the length of the \p from
 *        argument, and its number of columns is the number of
 *        vertices in the \p to argument. One row of the matrix shows the
 *        distances from/to a given vertex to the ones in \p to.
 *        For the unreachable vertices \c IGRAPH_INFINITY is returned.
 * \param from The source vertices.
 * \param to The target vertices. It is not allowed to include a
 *    vertex twice or more.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for
 *          the computation.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary
 *           data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n(|V|+|E|)),
 * n is the number of vertices to calculate,
 * |V| and |E| are the number of vertices and edges in the graph.
 *
 * \sa \ref igraph_get_shortest_paths() to get the paths themselves,
 * \ref igraph_distances_dijkstra() for the weighted version with non-negative
 * weights, \ref igraph_distances_bellman_ford() if you also have negative
 * weights.
 *
 * \example examples/simple/distances.c
 */
igraph_error_t igraph_distances(const igraph_t *graph, igraph_matrix_t *res,
                                const igraph_vs_t from, const igraph_vs_t to,
                                igraph_neimode_t mode) {
    return igraph_distances_cutoff(graph, res, from, to, mode, -1);
}

/**
 * \function igraph_shortest_paths
 * \brief Length of the shortest paths between vertices.
 *
 * \deprecated-by igraph_distances 0.10.0
 */
igraph_error_t igraph_shortest_paths(const igraph_t *graph,
                                     igraph_matrix_t *res,
                                     const igraph_vs_t from,
                                     const igraph_vs_t to,
                                     igraph_neimode_t mode) {
    return igraph_distances(graph, res, from, to, mode);
}

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths
 * \brief Shortest paths from a vertex.
 *
 * Finds unweighted shortest paths from a single source vertex to the specified
 * sets of target vertices. If there is more than one geodesic between two vertices,
 * this function gives only one of them. Use \ref igraph_get_all_shortest_paths()
 * to find \em all shortest paths.
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
 * \param from The ID of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the IDs of the vertices to/from which the
 *        shortest paths will be calculated. A vertex might be given multiple
 *        times.
 * \param mode The type of shortest paths to be used for the
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
 * \param parents A pointer to an initialized igraph vector or \c NULL.
 *        If not \c NULL, a vector containing the parent of each vertex in
 *        the single source shortest path tree is returned here. The
 *        parent of vertex \c i in the tree is the vertex from which vertex \c i
 *        was reached. The parent of the start vertex (in the \p from
 *        argument) is -1. If the parent is -2, it means
 *        that the given vertex was not reached from the source during the
 *        search. Note that the search terminates if all the vertices in
 *        \p to are reached.
 * \param inbound_edges A pointer to an initialized igraph vector or \c NULL.
 *        If not \c NULL, a vector containing the inbound edge of each vertex in
 *        the single source shortest path tree is returned here. The
 *        inbound edge of vertex \c i in the tree is the edge via which vertex \c i
 *        was reached. The start vertex and vertices that were not reached
 *        during the search will have -1 in the corresponding entry of the
 *        vector. Note that the search terminates if all the vertices in
 *        \p to are reached.
 *
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
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \sa \ref igraph_distances() if you only need the path lengths but
 * not the paths themselves; \ref igraph_get_shortest_paths_dijkstra()
 * for the weighted version; \ref igraph_get_all_shortest_paths() to
 * return all shortest paths between (source, target) pairs.
 *
 * \example examples/simple/igraph_get_shortest_paths.c
 */
igraph_error_t igraph_get_shortest_paths(const igraph_t *graph,
                              igraph_vector_int_list_t *vertices,
                              igraph_vector_int_list_t *edges,
                              igraph_integer_t from, const igraph_vs_t to,
                              igraph_neimode_t mode,
                              igraph_vector_int_t *parents,
                              igraph_vector_int_t *inbound_edges) {

    /* TODO: use inclist_t if to is long (longer than 1?) */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t *parent_eids;

    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;

    igraph_integer_t i, j, vsize;
    igraph_vector_int_t tmp = IGRAPH_VECTOR_NULL;

    igraph_vit_t vit;

    igraph_integer_t to_reach;
    igraph_integer_t reached = 0;

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("Index of source vertex is out of range.", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument.", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    if (vertices) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(vertices, IGRAPH_VIT_SIZE(vit)));
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(edges, IGRAPH_VIT_SIZE(vit)));
    }

    parent_eids = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(parent_eids, "Insufficient memory for shortest path calculation.");
    IGRAPH_FINALLY(igraph_free, parent_eids);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    /* Mark the vertices we need to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (parent_eids[ IGRAPH_VIT_GET(vit) ] == 0) {
            parent_eids[ IGRAPH_VIT_GET(vit) ] = -1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    /* Meaning of parent_eids[i]:
     *
     * - If parent_eids[i] < 0, it means that vertex i has to be reached and has not
     *   been reached yet.
     *
     * - If parent_eids[i] = 0, it means that vertex i does not have to be reached and
     *   it has not been reached yet.
     *
     * - If parent_eids[i] = 1, it means that vertex i is the start vertex.
     *
     * - Otherwise, parent_eids[i] is the ID of the edge from which vertex i was
     *   reached plus 2.
     */

    IGRAPH_CHECK(igraph_dqueue_int_push(&q, from + 1));
    if (parent_eids[ from ] < 0) {
        reached++;
    }
    parent_eids[ from ] = 1;

    while (!igraph_dqueue_int_empty(&q) && reached < to_reach) {
        igraph_integer_t act = igraph_dqueue_int_pop(&q) - 1;

        IGRAPH_CHECK(igraph_incident(graph, &tmp, act, mode));
        vsize = igraph_vector_int_size(&tmp);
        for (j = 0; j < vsize; j++) {
            igraph_integer_t edge = VECTOR(tmp)[j];
            igraph_integer_t neighbor = IGRAPH_OTHER(graph, edge, act);
            if (parent_eids[neighbor] > 0) {
                continue;
            } else if (parent_eids[neighbor] < 0) {
                reached++;
            }
            parent_eids[neighbor] = edge + 2;
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor + 1));
        }
    }

    if (reached < to_reach) {
        IGRAPH_WARNING("Couldn't reach some vertices");
    }

    /* Create `parents' if needed */
    if (parents) {
        IGRAPH_CHECK(igraph_vector_int_resize(parents, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (parent_eids[i] <= 0) {
                /* i was not reached */
                VECTOR(*parents)[i] = -2;
            } else if (parent_eids[i] == 1) {
                /* i is the start vertex */
                VECTOR(*parents)[i] = -1;
            } else {
                /* i was reached via the edge with ID = parent_eids[i] - 2 */
                VECTOR(*parents)[i] = IGRAPH_OTHER(graph, parent_eids[i] - 2, i);
            }
        }
    }

    /* Create `inbound_edges' if needed */
    if (inbound_edges) {
        IGRAPH_CHECK(igraph_vector_int_resize(inbound_edges, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (parent_eids[i] <= 1) {
                /* i was not reached or i is the start vertex */
                VECTOR(*inbound_edges)[i] = -1;
            } else {
                /* i was reached via the edge with ID = parent_eids[i] - 2 */
                VECTOR(*inbound_edges)[i] = parent_eids[i] - 2;
            }
        }
    }

    /* Create `vertices' and `edges' if needed */
    if (vertices || edges) {
        for (IGRAPH_VIT_RESET(vit), j = 0;
             !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), j++) {
            igraph_integer_t node = IGRAPH_VIT_GET(vit);
            igraph_vector_int_t *vvec = 0, *evec = 0;
            if (vertices) {
                vvec = igraph_vector_int_list_get_ptr(vertices, j);
                igraph_vector_int_clear(vvec);
            }
            if (edges) {
                evec = igraph_vector_int_list_get_ptr(edges, j);
                igraph_vector_int_clear(evec);
            }

            IGRAPH_ALLOW_INTERRUPTION();

            if (parent_eids[node] > 0) {
                igraph_integer_t act = node;
                igraph_integer_t size = 0;
                igraph_integer_t edge;
                while (parent_eids[act] > 1) {
                    size++;
                    edge = parent_eids[act] - 2;
                    act = IGRAPH_OTHER(graph, edge, act);
                }
                if (vvec) {
                    IGRAPH_CHECK(igraph_vector_int_resize(vvec, size + 1));
                    VECTOR(*vvec)[size] = node;
                }
                if (evec) {
                    IGRAPH_CHECK(igraph_vector_int_resize(evec, size));
                }
                act = node;
                while (parent_eids[act] > 1) {
                    size--;
                    edge = parent_eids[act] - 2;
                    act = IGRAPH_OTHER(graph, edge, act);
                    if (vvec) {
                        VECTOR(*vvec)[size] = act;
                    }
                    if (evec) {
                        VECTOR(*evec)[size] = edge;
                    }
                }
            }
        }
    }

    /* Clean */
    IGRAPH_FREE(parent_eids);
    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&tmp);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_shortest_path
 * \brief Shortest path from one vertex to another one.
 *
 * Calculates and returns a single unweighted shortest path from a
 * given vertex to another one. If there is more than one shortest
 * path between the two vertices, then an arbitrary one is returned.
 *
 * </para><para>
 * This function is a wrapper to \ref igraph_get_shortest_paths()
 * for the special case when only one target vertex is considered.
 *
 * \param graph The input graph, it can be directed or
 *        undirected. Directed paths are considered in directed
 *        graphs.
 * \param vertices Pointer to an initialized vector or a null
 *        pointer. If not a null pointer, then the vertex IDs along
 *        the path are stored here, including the source and target
 *        vertices.
 * \param edges Pointer to an initialized vector or a null
 *        pointer. If not a null pointer, then the edge IDs along the
 *        path are stored here.
 * \param from The ID of the source vertex.
 * \param to The ID of the target vertex.
 * \param mode A constant specifying how edge directions are
 *        considered in directed graphs. Valid modes are:
 *        \c IGRAPH_OUT, follows edge directions;
 *        \c IGRAPH_IN, follows the opposite directions; and
 *        \c IGRAPH_ALL, ignores edge directions. This argument is
 *        ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges in the graph.
 *
 * \sa \ref igraph_get_shortest_paths() for the version with more target
 * vertices.
 */

igraph_error_t igraph_get_shortest_path(const igraph_t *graph,
                             igraph_vector_int_t *vertices,
                             igraph_vector_int_t *edges,
                             igraph_integer_t from,
                             igraph_integer_t to,
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

    IGRAPH_CHECK(igraph_get_shortest_paths(graph, vp, ep, from,
                                           igraph_vss_1(to), mode, NULL, NULL));

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
