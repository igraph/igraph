/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2020  The igraph development team <igraph@igraph.org>

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

#include "igraph_centrality.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_dqueue.h"

#include "core/indheap.h"
#include "core/interruption.h"
#include "core/math.h"

/***** Closeness centrality *****/

/**
 * \ingroup structural
 * \function igraph_closeness
 * \brief Closeness centrality calculations for some vertices.
 *
 * The closeness centrality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as
 * the inverse of the mean distance to (or from) all other vertices.
 *
 * </para><para>
 * Closeness centrality is meaningful only for connected graphs.
 * If the graph is not connected, igraph computes the inverse of the
 * mean distance to (or from) all \em reachable vertices. In undirected
 * graphs, this is equivalent to computing the closeness separately in
 * each connected component. The optional \p all_reachable output
 * parameter is provided to help detect when the graph is disconnected.
 *
 * </para><para>
 * While there is no universally adopted definition of closeness centrality
 * for disconnected graphs, there have been some attempts for generalizing
 * the concept to the disconnected case. One type of approach considers the mean distance
 * only to reachable vertices, then re-scales the obtained certrality score
 * by a factor that depends on the number of reachable vertices
 * (i.e. the size of the component in the undirected case).
 * To facilitate computing these generalizations of closeness centrality,
 * the number of reachable vertices (not including the starting vertex)
 * is returned in \p reachable_count.
 *
 * </para><para>
 * In disconnected graphs, consider using the harmonic centrality,
 * computable using \ref igraph_harmonic_centrality().
 *
 * </para><para>
 * For isolated vertices, i.e. those having no associated paths, NaN is returned.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        closeness centrality scores for the given vertices.
 * \param reachable_count If not \c NULL, this vector will contain the number of
 *        vertices reachable from each vertex for which the closeness is calculated
 *        (not including that vertex).
 * \param all_reachable Pointer to a Boolean. If not \c NULL, it indicates if all
 *        vertices of the graph were reachable from each vertex in \p vids.
 *        If false, the graph is non-connected If true, and the graph is undirected,
 *        or if the graph is directed and \p vids contains all vertices, then the
 *        graph is connected.
 * \param vids The vertices for which the closeness centrality will be computed.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \param weights An optional vector containing edge weights for
 *        weighted closeness. No edge weight may be NaN. Supply a null
 *        pointer here for traditional, unweighted closeness.
 * \param normalized If true, the inverse of the mean distance to reachable
 *        vetices is returned. If false, the inverse of the sum of distances
 *        is returned.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|),
 * n is the number
 * of vertices for which the calculation is done and
 * |E| is the number
 * of edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_betweenness(),
 *   \ref igraph_harmonic_centrality().
 *   See \ref igraph_closeness_cutoff() for the range-limited closeness centrality.
 */
int igraph_closeness(const igraph_t *graph, igraph_vector_t *res,
                     igraph_vector_t *reachable_count, igraph_bool_t *all_reachable,
                     const igraph_vs_t vids, igraph_neimode_t mode,
                     const igraph_vector_t *weights,
                     igraph_bool_t normalized) {
    return igraph_closeness_cutoff(graph, res, reachable_count, all_reachable, vids, mode, weights, normalized, -1);
}

static int igraph_i_closeness_cutoff_weighted(const igraph_t *graph,
                                                igraph_vector_t *res,
                                                igraph_vector_t *reachable_count,
                                                igraph_bool_t *all_reachable,
                                                const igraph_vs_t vids,
                                                igraph_neimode_t mode,
                                                igraph_real_t cutoff,
                                                const igraph_vector_t *weights,
                                                igraph_bool_t normalized) {

    /* See igraph_shortest_paths_dijkstra() for the implementation
       details and the dirty tricks. */

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);

    igraph_2wheap_t Q;
    igraph_vit_t vit;
    long int nodes_to_calc;

    igraph_lazy_inclist_t inclist;
    long int i, j;

    igraph_vector_t dist;
    igraph_vector_long_t which;
    long int nodes_reached;

    igraph_real_t mindist = 0;

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (no_of_edges > 0) {
        igraph_real_t minweight = igraph_vector_min(weights);
        if (minweight <= 0) {
            IGRAPH_ERROR("Weight vector must be positive.", IGRAPH_EINVAL);
        } else if (igraph_is_nan(minweight)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    if (reachable_count) {
        igraph_vector_resize(reachable_count, nodes_to_calc);
    }

    if (all_reachable) {
        *all_reachable = 1; /* be optimistic */
    }

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_long_init(&which, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &which);

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    igraph_vector_null(res);

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {

        long int source = IGRAPH_VIT_GET(vit);
        igraph_2wheap_clear(&Q);
        igraph_2wheap_push_with_index(&Q, source, -1.0);
        VECTOR(which)[source] = i + 1;
        VECTOR(dist)[source] = 1.0;     /* actual distance is zero but we need to store distance + 1 */
        nodes_reached = 0;

        while (!igraph_2wheap_empty(&Q)) {
            igraph_integer_t minnei = (igraph_integer_t) igraph_2wheap_max_index(&Q);
            /* Now check all neighbors of minnei for a shorter path */
            igraph_vector_int_t *neis = igraph_lazy_inclist_get(&inclist, minnei);
            long int nlen = igraph_vector_int_size(neis);

            mindist = -igraph_2wheap_delete_max(&Q);

            if (cutoff >= 0 && (mindist - 1.0) > cutoff) {
                continue;    /* NOT break!!! */
            }

            VECTOR(*res)[i] += (mindist - 1.0);
            nodes_reached++;

            for (j = 0; j < nlen; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                long int to = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_real_t curdist = VECTOR(dist)[to];

                if (VECTOR(which)[to] != i + 1) {
                    /* First non-infinite distance */
                    VECTOR(which)[to] = i + 1;
                    VECTOR(dist)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
                } else if (curdist == 0 /* this means curdist is infinity */ || altdist < curdist) {
                    /* This is a shorter path */
                    VECTOR(dist)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

        if (reachable_count) {
            VECTOR(*reachable_count)[i] = nodes_reached - 1;
        }

        if (normalized) {
            /* compute the inverse of the average distance, considering only reachable nodes */
            VECTOR(*res)[i] = VECTOR(*res)[i] == 0 ? IGRAPH_NAN : ((igraph_real_t) (nodes_reached-1)) / VECTOR(*res)[i];
        } else {
            /* compute the inverse of the sum of distances */
            VECTOR(*res)[i] = VECTOR(*res)[i] == 0 ? IGRAPH_NAN : 1.0 / VECTOR(*res)[i];
        }

        if (all_reachable) {
            if (nodes_reached < no_of_nodes) {
                *all_reachable = 0 /* false */;
            }
        }
    } /* !IGRAPH_VIT_END(vit) */

    igraph_vector_long_destroy(&which);
    igraph_vector_destroy(&dist);
    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_closeness_estimate
 * \brief Closeness centrality estimations for some vertices.
 *
 * \deprecated-by igraph_closeness_cutoff 0.9
 *
 * </para><para>
 * The closeness centrality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as
 * the number of vertices minus one divided by the sum of the
 * lengths of all geodesics from/to the given vertex. When estimating
 * closeness centrality, igraph considers paths having a length less than
 * or equal to a prescribed cutoff value.
 *
 * </para><para>
 * If the graph is not connected, and there is no such path between two
 * vertices, the number of vertices is used instead the length of the
 * geodesic. This is always longer than the longest possible geodesic.
 *
 * </para><para>
 * Since the estimation considers vertex pairs with a distance greater than
 * the given value as disconnected, the resulting estimation will always be
 * lower than the actual closeness centrality.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        closeness centrality scores for the given vertices.
 * \param vids The vertices for which the closeness centrality will be estimated.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \param cutoff The maximal length of paths that will be considered.
 *        If negative, the exact closeness will be calculated (no upper
 *        limit on path lengths).
 * \param weights An optional vector containing edge weights for
 *        weighted closeness. No edge weight may be NaN. Supply a
 *        null pointer here for traditional, unweighted closeness.
 * \param normalized Boolean, whether to normalize results by multiplying
 *        by the number of vertices minus one.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|),
 * n is the number
 * of vertices for which the calculation is done and
 * |E| is the number
 * of edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_betweenness().
 */

int igraph_closeness_estimate(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_vs_t vids, igraph_neimode_t mode,
                              igraph_real_t cutoff,
                              const igraph_vector_t *weights,
                              igraph_bool_t normalized) {
    IGRAPH_WARNING("igraph_closeness_estimate is deprecated, use igraph_closeness_cutoff.");
    return igraph_closeness_cutoff(graph, res, NULL, NULL, vids, mode, weights, normalized, cutoff);
}


/**
 * \ingroup structural
 * \function igraph_closeness_cutoff
 * \brief Range limited closeness centrality.
 *
 * This function computes a range-limited version of closeness centrality
 * by considering only those shortest paths whose length is no greater
 * then the given cutoff value.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        range-limited closeness centrality scores for the given vertices.
 * \param reachable_count If not \c NULL, this vector will contain the number of
 *        vertices reachable within the cutoff distance from each vertex for which
 *        the range-limited closeness is calculated (not including that vertex).
 * \param all_reachable Pointer to a Boolean. If not \c NULL, it indicates if all
 *        vertices of the graph were reachable from each vertex in \p vids within
 *        the given cutoff distance.
 * \param vids The vertices for which the range limited closeness centrality
 *             will be computed.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \param weights An optional vector containing edge weights for
 *        weighted closeness. No edge weight may be NaN. Supply a null
 *        pointer here for traditional, unweighted closeness.
 * \param normalized If true, the inverse of the mean distance to vertices
 *        reachable within the cutoff is returned. If false, the inverse
 *        of the sum of distances is returned.
 * \param cutoff The maximal length of paths that will be considered.
 *        If negative, the exact closeness will be calculated (no upper
 *        limit on path lengths).
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|),
 * n is the number
 * of vertices for which the calculation is done and
 * |E| is the number
 * of edges in the graph.
 *
 * \sa \ref igraph_closeness() to calculate the exact closeness centrality.
 */

int igraph_closeness_cutoff(const igraph_t *graph, igraph_vector_t *res,
                            igraph_vector_t *reachable_count, igraph_bool_t *all_reachable,
                            const igraph_vs_t vids, igraph_neimode_t mode,
                            const igraph_vector_t *weights,
                            igraph_bool_t normalized,
                            igraph_real_t cutoff) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t already_counted;
    igraph_vector_int_t *neis;
    long int i, j;
    long int nodes_reached;
    igraph_adjlist_t allneis;

    long int actdist = 0;

    igraph_dqueue_t q;

    long int nodes_to_calc;
    igraph_vit_t vit;

    if (weights) {
        return igraph_i_closeness_cutoff_weighted(graph, res, reachable_count, all_reachable, vids, mode, cutoff,
                weights, normalized);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    if (reachable_count) {
        igraph_vector_resize(reachable_count, nodes_to_calc);
    }

    if (all_reachable) {
        *all_reachable = 1; /* be optimistic */
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for closeness.", IGRAPH_EINVMODE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&already_counted, no_of_nodes);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    igraph_vector_null(res);

    for (IGRAPH_VIT_RESET(vit), i = 0;
         !IGRAPH_VIT_END(vit);
         IGRAPH_VIT_NEXT(vit), i++) {
        nodes_reached = 0;

        igraph_dqueue_clear(&q);
        IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(vit)));
        IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
        VECTOR(already_counted)[(long int)IGRAPH_VIT_GET(vit)] = i + 1;

        IGRAPH_PROGRESS("Closeness: ", 100.0 * i / nodes_to_calc, NULL);
        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int act = (long int) igraph_dqueue_pop(&q);
            actdist = (long int) igraph_dqueue_pop(&q);

            if (cutoff >= 0 && actdist > cutoff) {
                continue;    /* NOT break!!! */
            }

            VECTOR(*res)[i] += actdist;
            nodes_reached++;

            /* check the neighbors */
            neis = igraph_adjlist_get(&allneis, act);
            for (j = 0; j < igraph_vector_int_size(neis); j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
                if (VECTOR(already_counted)[neighbor] == i + 1) {
                    continue;
                }
                VECTOR(already_counted)[neighbor] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
            }
        }

        if (reachable_count) {
            VECTOR(*reachable_count)[i] = nodes_reached - 1;
        }

        if (normalized) {
            /* compute the inverse of the average distance, considering only reachable nodes */
            VECTOR(*res)[i] = VECTOR(*res)[i] == 0 ? IGRAPH_NAN : ((igraph_real_t) (nodes_reached-1)) / VECTOR(*res)[i];
        } else {
            /* compute the inverse of the sum of distances */
            VECTOR(*res)[i] = VECTOR(*res)[i] == 0 ? IGRAPH_NAN : 1.0 / VECTOR(*res)[i];
        }

        if (all_reachable) {
            if (nodes_reached < no_of_nodes) {
                *all_reachable = 0 /* false */;
            }
        }
    }

    IGRAPH_PROGRESS("Closeness: ", 100.0, NULL);

    /* Clean */
    igraph_dqueue_destroy(&q);
    igraph_vector_destroy(&already_counted);
    igraph_vit_destroy(&vit);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


/***** Harmonic centrality *****/

static int igraph_i_harmonic_centrality_unweighted(const igraph_t *graph, igraph_vector_t *res,
                                                   const igraph_vs_t vids, igraph_neimode_t mode,
                                                   igraph_bool_t normalized,
                                                   igraph_real_t cutoff) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t already_counted;
    igraph_vector_int_t *neis;
    long int i, j;
    igraph_adjlist_t allneis;

    long int actdist = 0;

    igraph_dqueue_t q;

    long int nodes_to_calc;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for harmonic centrality.", IGRAPH_EINVMODE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&already_counted, no_of_nodes);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    igraph_vector_null(res);

    for (IGRAPH_VIT_RESET(vit), i = 0;
         !IGRAPH_VIT_END(vit);
         IGRAPH_VIT_NEXT(vit), i++)
    {
        long int source = IGRAPH_VIT_GET(vit);

        igraph_dqueue_clear(&q);
        IGRAPH_CHECK(igraph_dqueue_push(&q, source));
        IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
        VECTOR(already_counted)[source] = i + 1;

        IGRAPH_PROGRESS("Harmonic centrality: ", 100.0 * i / nodes_to_calc, NULL);
        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int act = (long int) igraph_dqueue_pop(&q);
            actdist = (long int) igraph_dqueue_pop(&q);

            if (cutoff >= 0 && actdist > cutoff) {
                continue;    /* NOT break!!! */
            }

            /* Exclude self-distance, which is zero. */
            if (source != act) {
                VECTOR(*res)[i] += 1.0/actdist;
            }

            /* check the neighbors */
            neis = igraph_adjlist_get(&allneis, act);
            for (j = 0; j < igraph_vector_int_size(neis); j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
                if (VECTOR(already_counted)[neighbor] == i + 1) {
                    continue;
                }
                VECTOR(already_counted)[neighbor] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
            }
        }
    }

    if (normalized && no_of_nodes > 1 /* not a null graph or singleton graph */) {
        igraph_vector_scale(res, 1.0 / (no_of_nodes - 1));
    }

    IGRAPH_PROGRESS("Harmonic centrality: ", 100.0, NULL);

    /* Clean */
    igraph_dqueue_destroy(&q);
    igraph_vector_destroy(&already_counted);
    igraph_vit_destroy(&vit);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


static int igraph_i_harmonic_centrality_weighted(const igraph_t *graph,
                                                 igraph_vector_t *res,
                                                 const igraph_vs_t vids,
                                                 igraph_neimode_t mode,
                                                 const igraph_vector_t *weights,
                                                 igraph_bool_t normalized,
                                                 igraph_real_t cutoff) {

    /* See igraph_shortest_paths_dijkstra() for the implementation
       details and the dirty tricks. */

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);

    igraph_2wheap_t Q;
    igraph_vit_t vit;
    long int nodes_to_calc;

    igraph_lazy_inclist_t inclist;
    long int i, j;

    igraph_vector_t dist;
    igraph_vector_long_t which;

    igraph_real_t mindist = 0;

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (no_of_edges > 0) {
        igraph_real_t minweight = igraph_vector_min(weights);
        if (minweight <= 0) {
            IGRAPH_ERROR("Weight vector must be positive.", IGRAPH_EINVAL);
        } else if (igraph_is_nan(minweight)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_long_init(&which, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &which);

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    igraph_vector_null(res);

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {

        long int source = IGRAPH_VIT_GET(vit);
        igraph_2wheap_clear(&Q);
        igraph_2wheap_push_with_index(&Q, source, -1.0);
        VECTOR(which)[source] = i + 1;
        VECTOR(dist)[source] = 1.0;     /* actual distance is zero but we need to store distance + 1 */

        while (!igraph_2wheap_empty(&Q)) {
            igraph_integer_t minnei = (igraph_integer_t) igraph_2wheap_max_index(&Q);
            /* Now check all neighbors of minnei for a shorter path */
            igraph_vector_int_t *neis = igraph_lazy_inclist_get(&inclist, minnei);
            long int nlen = igraph_vector_int_size(neis);

            mindist = -igraph_2wheap_delete_max(&Q);

            if (cutoff >= 0 && (mindist - 1.0) > cutoff) {
                continue;    /* NOT break!!! */
            }

            /* Exclude self-distance, which is zero. */
            if (source != minnei) {
                VECTOR(*res)[i] += 1.0 / (mindist - 1.0);
            }

            for (j = 0; j < nlen; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                long int to = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_real_t curdist = VECTOR(dist)[to];

                if (VECTOR(which)[to] != i + 1) {
                    /* First non-infinite distance */
                    VECTOR(which)[to] = i + 1;
                    VECTOR(dist)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
                } else if (curdist == 0 /* this means curdist is infinity */ || altdist < curdist) {
                    /* This is a shorter path */
                    VECTOR(dist)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

    } /* !IGRAPH_VIT_END(vit) */

    if (normalized && no_of_nodes > 1 /* not a null graph or singleton graph */) {
        igraph_vector_scale(res, 1.0 / (no_of_nodes - 1));
    }

    igraph_vector_long_destroy(&which);
    igraph_vector_destroy(&dist);
    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_harmonic_centrality_cutoff
 * \brief Range limited harmonic centrality.
 *
 * This function computes the range limited version of harmonic centrality:
 * only those shortest paths are considered whose length is not above the given cutoff.
 * The inverse distance to vertices not reachable within the cutoff is considered
 * to be zero.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        range limited harmonic centrality scores for the given vertices.
 * \param vids The vertices for which the harmonic centrality will be computed.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \param weights An optional vector containing edge weights for
 *        weighted harmonic centrality. No edge weight may be NaN.
 *        If \c NULL, all weights are considered to be one.
 * \param normalized Boolean, whether to normalize the result. If true,
 *        the result is the mean inverse path length to other vertices.
 *        i.e. it is normalized by the number of vertices minus one.
 *        If false, the result is the sum of inverse path lengths to other
 *        vertices.
 * \param cutoff The maximal length of paths that will be considered.
 *        The inverse distance to vertices that are not reachable within
 *        the cutoff path length is considered to be zero.
 *        Supply a negative value to compute the exact harmonic centrality,
 *        without any upper limit on the length of paths.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|), where
 * n is the number of vertices for which the calculation is done and
 * |E| is the number of edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_closeness(), \ref igraph_betweenness().
 */

int igraph_harmonic_centrality_cutoff(const igraph_t *graph, igraph_vector_t *res,
                                      const igraph_vs_t vids, igraph_neimode_t mode,
                                      const igraph_vector_t *weights,
                                      igraph_bool_t normalized,
                                      igraph_real_t cutoff) {
    if (weights) {
        return igraph_i_harmonic_centrality_weighted(graph, res, vids, mode, weights, normalized, cutoff);
    } else {
        return igraph_i_harmonic_centrality_unweighted(graph, res, vids, mode, normalized, cutoff);
    }
}


/**
 * \ingroup structural
 * \function igraph_harmonic_centrality
 * \brief Harmonic centrality for some vertices.
 *
 * The harmonic centrality of a vertex is the mean inverse distance to
 * all other vertices. The inverse distance to an unreachable vertex
 * is considered to be zero.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * M. Marchiori and V. Latora, Harmony in the small-world, Physica A 285, pp. 539-546 (2000).
 * https://doi.org/10.1016/S0378-4371%2800%2900311-3
 *
 * </para><para>
 * Y. Rochat, Closeness Centrality Extended to Unconnected Graphs: the Harmonic Centrality Index, ASNA 2009.
 * https://infoscience.epfl.ch/record/200525
 *
 * </para><para>
 * S. Vigna and P. Boldi, Axioms for Centrality, Internet Mathematics 10, (2014).
 * https://doi.org/10.1080/15427951.2013.865686
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        harmonic centrality scores for the given vertices.
 * \param vids The vertices for which the harmonic centrality will be computed.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \param weights An optional vector containing edge weights for
 *        weighted harmonic centrality. No edge weight may be NaN.
 *        If \c NULL, all weights are considered to be one.
 * \param normalized Boolean, whether to normalize the result. If true,
 *        the result is the mean inverse path length to other vertices,
 *        i.e. it is normalized by the number of vertices minus one.
 *        If false, the result is the sum of inverse path lengths to other
 *        vertices.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|), where
 * n is the numberof vertices for which the calculation is done and
 * |E| is the number of edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_closeness(), \ref igraph_degree(), \ref igraph_betweenness().
 */

int igraph_harmonic_centrality(const igraph_t *graph, igraph_vector_t *res,
                               const igraph_vs_t vids, igraph_neimode_t mode,
                               const igraph_vector_t *weights,
                               igraph_bool_t normalized) {
    return igraph_harmonic_centrality_cutoff(graph, res, vids, mode, weights, normalized, /* cutoff= */ -1);
}
