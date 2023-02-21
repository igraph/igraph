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

#include "igraph_neighborhood.h"

#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"

/**
 * \function igraph_neighborhood_size
 * \brief Calculates the size of the neighborhood of a given vertex.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. I.e., order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para>
 * This function calculates the size of the neighborhood
 * of the given order for the given vertices.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result will be
 *    stored here. It will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \c order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \c order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   Vertices reachable with a path shorter than this value are excluded.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood() for calculating the actual neighborhood,
 * \ref igraph_neighborhood_graphs() for creating separate graphs from
 * the neighborhoods.
 *
 * Time complexity: O(n*d*o), where n is the number vertices for which
 * the calculation is performed, d is the average degree, o is the order.
 */
igraph_error_t igraph_neighborhood_size(const igraph_t *graph, igraph_vector_int_t *res,
                             igraph_vs_t vids, igraph_integer_t order,
                             igraph_neimode_t mode,
                             igraph_integer_t mindist) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vit_t vit;
    igraph_integer_t i, j;
    igraph_integer_t *added;
    igraph_vector_int_t neis;

    if (order < 0) {
        IGRAPH_ERRORF("Negative order in neighborhood size: %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, order);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERRORF("Minimum distance should be between 0 and the neighborhood order (%" IGRAPH_PRId "), got %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, order, mindist);
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(added, "Cannot calculate neighborhood size.");
    IGRAPH_FINALLY(igraph_free, added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_vector_int_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        igraph_integer_t size = mindist == 0 ? 1 : 0;
        added[node] = i + 1;
        igraph_dqueue_int_clear(&q);
        if (order > 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, node));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            igraph_integer_t n;
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
            n = igraph_vector_int_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            size++;
                        }
                    }
                }
            } else {
                /* we just count them, but don't add them */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            size++;
                        }
                    }
                }
            }

        } /* while q not empty */

        VECTOR(*res)[i] = size;
    } /* for VIT, i */

    igraph_vector_int_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_neighborhood
 * \brief Calculate the neighborhood of vertices.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. I.e., order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para>
 * This function calculates the vertices within the
 * neighborhood of the specified vertices.
 *
 * \param graph The input graph.
 * \param res An initialized list of integer vectors. The result of the
 *    calculation will be stored here. The list will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \p order steps are included. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \p order steps are included. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   Vertices reachable with a path shorter than this value are excluded.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood_size() to calculate the size of the
 * neighborhood, \ref igraph_neighborhood_graphs() for creating
 * graphs from the neighborhoods.
 *
 * Time complexity: O(n*d*o), n is the number of vertices for which
 * the calculation is performed, d is the average degree, o is the
 * order.
 */
igraph_error_t igraph_neighborhood(const igraph_t *graph, igraph_vector_int_list_t *res,
                        igraph_vs_t vids, igraph_integer_t order,
                        igraph_neimode_t mode, igraph_integer_t mindist) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vit_t vit;
    igraph_integer_t i, j;
    igraph_integer_t *added;
    igraph_vector_int_t neis;
    igraph_vector_int_t tmp;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(added, "Cannot calculate neighborhood size.");
    IGRAPH_FINALLY(igraph_free, added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_int_list_reserve(res, IGRAPH_VIT_SIZE(vit)));
    igraph_vector_int_list_clear(res);

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_int_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, node));
        }
        if (order > 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, node));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            igraph_integer_t n;
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
            n = igraph_vector_int_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(res, &tmp));
    }

    igraph_vector_int_destroy(&tmp);
    igraph_vector_int_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_neighborhood_graphs
 * \brief Create graphs from the neighborhood(s) of some vertex/vertices.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. Ie. order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para>
 * This function finds every vertex in the neighborhood
 * of a given parameter vertex and creates the induced subgraph from these
 * vertices.
 *
 * </para><para>
 * The first version of this function was written by
 * Vincent Matossian, thanks Vincent.
 * \param graph The input graph.
 * \param res Pointer to a list of graphs, the result will be stored
 *   here. Each item in the list is an \c igraph_t object. The list will be
 *   resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \p order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \p order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   Vertices reachable with a path shorter than this value are excluded.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood_size() for calculating the neighborhood
 * sizes only, \ref igraph_neighborhood() for calculating the
 * neighborhoods (but not creating graphs).
 *
 * Time complexity: O(n*(|V|+|E|)), where n is the number vertices for
 * which the calculation is performed, |V| and |E| are the number of
 * vertices and edges in the original input graph.
 */
igraph_error_t igraph_neighborhood_graphs(const igraph_t *graph, igraph_graph_list_t *res,
                               igraph_vs_t vids, igraph_integer_t order,
                               igraph_neimode_t mode,
                               igraph_integer_t mindist) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vit_t vit;
    igraph_integer_t i, j;
    igraph_integer_t *added;
    igraph_vector_int_t neis;
    igraph_vector_int_t tmp;
    igraph_t newg;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(added, "Cannot calculate neighborhood size");
    IGRAPH_FINALLY(igraph_free, added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);

    igraph_graph_list_clear(res);

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_int_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, node));
        }
        if (order > 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, node));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            igraph_integer_t n;
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
            n = igraph_vector_int_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        if (igraph_vector_int_size(&tmp) < no_of_nodes) {
            IGRAPH_CHECK(igraph_induced_subgraph(graph, &newg,
                                                 igraph_vss_vector(&tmp),
                                                 IGRAPH_SUBGRAPH_AUTO));
        } else {
            IGRAPH_CHECK(igraph_copy(&newg, graph));
        }

        IGRAPH_FINALLY(igraph_destroy, &newg);
        IGRAPH_CHECK(igraph_graph_list_push_back(res, &newg));
        IGRAPH_FINALLY_CLEAN(1);  /* ownership of `newg' taken by `res' */
    }

    igraph_vector_int_destroy(&tmp);
    igraph_vector_int_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}
