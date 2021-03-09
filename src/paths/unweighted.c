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
 * \function igraph_shortest_paths
 * \brief The length of the shortest paths between vertices.
 *
 * \param graph The graph object.
 * \param res The result of the calculation, a matrix. A pointer to an
 *        initialized matrix, to be more precise. The matrix will be
 *        resized if needed. It will have the same
 *        number of rows as the length of the \c from
 *        argument, and its number of columns is the number of
 *        vertices in the \c to argument. One row of the matrix shows the
 *        distances from/to a given vertex to the ones in \c to.
 *        For the unreachable vertices IGRAPH_INFINITY is returned.
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
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n(|V|+|E|)),
 * n is the
 * number of vertices to calculate, |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa \ref igraph_get_shortest_paths() to get the paths themselves,
 * \ref igraph_shortest_paths_dijkstra() for the weighted version.
 */
int igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res,
                          const igraph_vs_t from, const igraph_vs_t to,
                          igraph_neimode_t mode) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_from, no_of_to;
    long int *already_counted;
    igraph_adjlist_t adjlist;
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_bool_t all_to;

    long int i, j;
    igraph_vit_t fromvit, tovit;
    igraph_real_t my_infinity = IGRAPH_INFINITY;
    igraph_vector_t indexv;

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);
    no_of_from = IGRAPH_VIT_SIZE(fromvit);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    already_counted = IGRAPH_CALLOC(no_of_nodes, long int);
    if (already_counted == 0) {
        IGRAPH_ERROR("shortest paths failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_counted);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

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
        IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(fromvit)));
        IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
        already_counted[ (long int) IGRAPH_VIT_GET(fromvit) ] = i + 1;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int act = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);

            if (all_to) {
                MATRIX(*res, i, act) = actdist;
            } else {
                if (VECTOR(indexv)[act]) {
                    MATRIX(*res, i, (long int)(VECTOR(indexv)[act] - 1)) = actdist;
                    reached++;
                    if (reached == no_of_to) {
                        igraph_dqueue_clear(&q);
                        break;
                    }
                }
            }

            neis = igraph_adjlist_get(&adjlist, act);
            for (j = 0; j < igraph_vector_int_size(neis); j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
                if (already_counted[neighbor] == i + 1) {
                    continue;
                }
                already_counted[neighbor] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
            }
        }
    }

    /* Clean */
    if (!all_to) {
        igraph_vit_destroy(&tovit);
        igraph_vector_destroy(&indexv);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_FREE(already_counted);
    igraph_dqueue_destroy(&q);
    igraph_vit_destroy(&fromvit);
    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths
 * \brief Shortest paths from a vertex.
 *
 * </para><para>
 * If there is more than one geodesic between two vertices, this
 * function gives only one of them.
 * \param graph The graph object.
 * \param vertices The result, the ids of the vertices along the paths.
 *        This is a pointer vector, each element points to a vector
 *        object. These should be initialized before passing them to
 *        the function, which will properly clear and/or resize them
 *        and fill the ids of the vertices along the geodesics from/to
 *        the vertices. Supply a null pointer here if you don't need
 *        these vectors.
 * \param edges The result, the ids of the edges along the paths.
 *        This is a pointer vector, each element points to a vector
 *        object. These should be initialized before passing them to
 *        the function, which will properly clear and/or resize them
 *        and fill the ids of the vertices along the geodesics from/to
 *        the vertices. Supply a null pointer here if you don't need
 *        these vectors.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the ids of the vertices to/from which the
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
 *
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
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \sa \ref igraph_shortest_paths() if you only need the path length but
 * not the paths themselves.
 *
 * \example examples/simple/igraph_get_shortest_paths.c
 */
int igraph_get_shortest_paths(const igraph_t *graph,
                              igraph_vector_ptr_t *vertices,
                              igraph_vector_ptr_t *edges,
                              igraph_integer_t from, const igraph_vs_t to,
                              igraph_neimode_t mode,
                              igraph_vector_long_t *predecessors,
                              igraph_vector_long_t *inbound_edges) {

    /* TODO: use inclist_t if to is long (longer than 1?) */

    long int no_of_nodes = igraph_vcount(graph);
    long int *father;

    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;

    long int i, j, vsize;
    igraph_vector_t tmp = IGRAPH_VECTOR_NULL;

    igraph_vit_t vit;

    long int to_reach;
    long int reached = 0;

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("cannot get shortest paths", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    if (vertices && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(vertices)) {
        IGRAPH_ERROR("Size of the `vertices' and the `to' should match", IGRAPH_EINVAL);
    }
    if (edges && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(edges)) {
        IGRAPH_ERROR("Size of the `edges' and the `to' should match", IGRAPH_EINVAL);
    }

    father = IGRAPH_CALLOC(no_of_nodes, long int);
    if (father == 0) {
        IGRAPH_ERROR("cannot get shortest paths", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, father);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    /* Mark the vertices we need to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (father[ (long int) IGRAPH_VIT_GET(vit) ] == 0) {
            father[ (long int) IGRAPH_VIT_GET(vit) ] = -1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    /* Meaning of father[i]:
     *
     * - If father[i] < 0, it means that vertex i has to be reached and has not
     *   been reached yet.
     *
     * - If father[i] = 0, it means that vertex i does not have to be reached and
     *   it has not been reached yet.
     *
     * - If father[i] = 1, it means that vertex i is the start vertex.
     *
     * - Otherwise, father[i] is the ID of the edge from which vertex i was
     *   reached plus 2.
     */

    IGRAPH_CHECK(igraph_dqueue_push(&q, from + 1));
    if (father[ (long int) from ] < 0) {
        reached++;
    }
    father[ (long int)from ] = 1;

    while (!igraph_dqueue_empty(&q) && reached < to_reach) {
        long int act = (long int) igraph_dqueue_pop(&q) - 1;

        IGRAPH_CHECK(igraph_incident(graph, &tmp, (igraph_integer_t) act, mode));
        vsize = igraph_vector_size(&tmp);
        for (j = 0; j < vsize; j++) {
            long int edge = (long int) VECTOR(tmp)[j];
            long int neighbor = IGRAPH_OTHER(graph, edge, act);
            if (father[neighbor] > 0) {
                continue;
            } else if (father[neighbor] < 0) {
                reached++;
            }
            father[neighbor] = edge + 2;
            IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor + 1));
        }
    }

    if (reached < to_reach) {
        IGRAPH_WARNING("Couldn't reach some vertices");
    }

    /* Create `predecessors' if needed */
    if (predecessors) {
        IGRAPH_CHECK(igraph_vector_long_resize(predecessors, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (father[i] <= 0) {
                /* i was not reached */
                VECTOR(*predecessors)[i] = -1;
            } else if (father[i] == 1) {
                /* i is the start vertex */
                VECTOR(*predecessors)[i] = i;
            } else {
                /* i was reached via the edge with ID = father[i] - 2 */
                VECTOR(*predecessors)[i] = IGRAPH_OTHER(graph, father[i] - 2, i);
            }
        }
    }

    /* Create `inbound_edges' if needed */
    if (inbound_edges) {
        IGRAPH_CHECK(igraph_vector_long_resize(inbound_edges, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (father[i] <= 1) {
                /* i was not reached or i is the start vertex */
                VECTOR(*inbound_edges)[i] = -1;
            } else {
                /* i was reached via the edge with ID = father[i] - 2 */
                VECTOR(*inbound_edges)[i] = father[i] - 2;
            }
        }
    }

    /* Create `vertices' and `edges' if needed */
    if (vertices || edges) {
        for (IGRAPH_VIT_RESET(vit), j = 0;
             !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), j++) {
            long int node = IGRAPH_VIT_GET(vit);
            igraph_vector_t *vvec = 0, *evec = 0;
            if (vertices) {
                vvec = VECTOR(*vertices)[j];
                igraph_vector_clear(vvec);
            }
            if (edges) {
                evec = VECTOR(*edges)[j];
                igraph_vector_clear(evec);
            }

            IGRAPH_ALLOW_INTERRUPTION();

            if (father[node] > 0) {
                long int act = node;
                long int size = 0;
                long int edge;
                while (father[act] > 1) {
                    size++;
                    edge = father[act] - 2;
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
                while (father[act] > 1) {
                    size--;
                    edge = father[act] - 2;
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
    IGRAPH_FREE(father);
    igraph_dqueue_destroy(&q);
    igraph_vector_destroy(&tmp);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/**
 * \function igraph_get_shortest_path
 * \brief Shortest path from one vertex to another one.
 *
 * Calculates and returns a single unweighted shortest path from a
 * given vertex to another one. If there are more than one shortest
 * paths between the two vertices, then an arbitrary one is returned.
 *
 * </para><para>This function is a wrapper to \ref
 * igraph_get_shortest_paths(), for the special case when only one
 * target vertex is considered.
 * \param graph The input graph, it can be directed or
 *        undirected. Directed paths are considered in directed
 *        graphs.
 * \param vertices Pointer to an initialized vector or a null
 *        pointer. If not a null pointer, then the vertex ids along
 *        the path are stored here, including the source and target
 *        vertices.
 * \param edges Pointer to an uninitialized vector or a null
 *        pointer. If not a null pointer, then the edge ids along the
 *        path are stored here.
 * \param from The id of the source vertex.
 * \param to The id of the target vertex.
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

int igraph_get_shortest_path(const igraph_t *graph,
                             igraph_vector_t *vertices,
                             igraph_vector_t *edges,
                             igraph_integer_t from,
                             igraph_integer_t to,
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

    IGRAPH_CHECK(igraph_get_shortest_paths(graph, vp, ep, from,
                                           igraph_vss_1(to), mode, 0, 0));

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
