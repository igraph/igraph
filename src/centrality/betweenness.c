/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2021  The igraph development team <igraph@igraph.org>

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

#include "igraph_memory.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_stack.h"
#include "igraph_dqueue.h"

#include "core/indheap.h"
#include "core/interruption.h"
#include "core/math.h"

static int igraph_i_sspf( const igraph_t *graph, long int source, igraph_vector_t *dist, 
    double *nrgeo,
    igraph_stack_t *S,
    igraph_inclist_t *fathers,
    igraph_inclist_t *inclist,
    igraph_real_t cutoff) {
    
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_dqueue_t Q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    long int nlen;

    IGRAPH_DQUEUE_INIT_FINALLY(&Q, 100);
    
    IGRAPH_ALLOW_INTERRUPTION();

    IGRAPH_CHECK(igraph_dqueue_push(&Q, source));
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_dqueue_empty(&Q)) {
        long int actnode = (long int) igraph_dqueue_pop(&Q);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && VECTOR(*dist)[actnode] > cutoff + 1) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            igraph_vector_int_clear(igraph_inclist_get(fathers, actnode));
            continue;
        }

        IGRAPH_CHECK(igraph_stack_push(S, actnode));
        neis = igraph_inclist_get(inclist, actnode);
        nlen = igraph_vector_int_size(neis);
        for (int j = 0; j < nlen; j++) {
            long int edge = (long int) VECTOR(*neis)[j];
            long int neighbor = IGRAPH_OTHER(graph, edge, actnode);
            if (VECTOR(*dist)[neighbor] == 0) {
                VECTOR(*dist)[neighbor] = VECTOR(*dist)[actnode] + 1;
                IGRAPH_CHECK(igraph_dqueue_push(&Q, neighbor));
            }
            if (VECTOR(*dist)[neighbor] == VECTOR(*dist)[actnode] + 1 &&
                (VECTOR(*dist)[neighbor] <= cutoff + 1 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                igraph_vector_int_t *v = igraph_inclist_get(fathers,
                                            neighbor);
                igraph_vector_int_push_back(v, edge);
                nrgeo[neighbor] += nrgeo[actnode];
            }
        }
    }
    igraph_dqueue_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_SUCCESS; 
}

static int igraph_i_sspf_weighted(
    const igraph_t *graph, long int source, igraph_vector_t *dist, 
    double *nrgeo, 
    const igraph_vector_t *weights,
    igraph_stack_t *S,
    igraph_inclist_t *fathers,
    igraph_inclist_t *inclist,
    igraph_real_t cutoff) {
    
    int cmp_result;
    igraph_2wheap_t Q;
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;
    
    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_ALLOW_INTERRUPTION();

    igraph_2wheap_push_with_index(&Q, source, -1.0);
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_2wheap_empty(&Q)) {
        long int minnei = igraph_2wheap_max_index(&Q);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;
        long int nlen;

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && mindist > cutoff + 1.0) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[minnei] = 0;
            nrgeo[minnei] = 0;
            igraph_vector_int_clear(igraph_inclist_get(fathers, minnei));
            continue;
        }

        igraph_stack_push(S, minnei);

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_inclist_get(inclist, minnei);
        nlen = igraph_vector_int_size(neis);
        for (int j = 0; j < nlen; j++) {
            long int edge = (long int) VECTOR(*neis)[j];
            long int to = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
            igraph_real_t curdist = VECTOR(*dist)[to];

            if (curdist == 0) {
                /* this means curdist is infinity */
                cmp_result = -1;
            } else {
                cmp_result = igraph_cmp_epsilon(altdist, curdist, eps);
            }

            if (curdist == 0) {
                /* This is the first non-infinite distance */
                igraph_vector_int_t *v = igraph_inclist_get(fathers, to);
                igraph_vector_int_resize(v, 1);
                VECTOR(*v)[0] = edge;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
            } else if (cmp_result < 0) {
                /* This is a shorter path */
                igraph_vector_int_t *v = igraph_inclist_get(fathers, to);
                igraph_vector_int_resize(v, 1);
                VECTOR(*v)[0] = edge;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
            } else if (cmp_result == 0 &&
                (altdist <= cutoff + 1.0 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                igraph_vector_int_t *v = igraph_inclist_get(fathers, to);
                igraph_vector_int_push_back(v, edge);
                nrgeo[to] += nrgeo[minnei];
            }
        }
    }
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_SUCCESS;
}

/***** Vertex betweenness *****/

/**
 * \ingroup structural
 * \function igraph_betweenness
 * \brief Betweenness centrality of some vertices.
 *
 * </para><para>
 * The betweenness centrality of a vertex is the number of geodesics
 * going through it. If there are more than one geodesic between two
 * vertices, the value of these geodesics are weighted by one over the
 * number of geodesics.
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        betweenness scores for the specified vertices.
 * \param vids The vertices of which the betweenness centrality scores
 *        will be calculated.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional vector containing edge weights for
 *        calculating weighted betweenness. No edge weight may be NaN.
 *        Supply a null pointer here for unweighted betweenness.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *        \c IGRAPH_EINVVID, invalid vertex id passed in
 *        \p vids.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 * Note that the time complexity is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph. See \ref igraph_betweenness_cutoff() to
 *     calculate the range-limited betweenness of the vertices in a graph.
 */
int igraph_betweenness(const igraph_t *graph, igraph_vector_t *res,
                       const igraph_vs_t vids, igraph_bool_t directed,
                       const igraph_vector_t* weights) {
    return igraph_betweenness_cutoff(graph, res, vids, directed, weights, -1);
}

/**
 * \ingroup structural
 * \function igraph_betweenness_cutoff
 * \brief Range-limited betweenness centrality.
 *
 * </para><para>
 * This function computes a range-limited version of betweenness centrality
 * by considering only those shortest paths whose length is no greater
 * then the given cutoff value.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        range-limited betweenness scores for the specified vertices.
 * \param vids The vertices for which the range-limited betweenness centrality
 *        scores will be computed.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional vector containing edge weights for
 *        calculating weighted betweenness. No edge weight may be NaN.
 *        Supply a null pointer here for unweighted betweenness.
 * \param cutoff The maximal length of paths that will be considered.
 *        If negative, the exact betweenness will be calculated, and
 *        there will be no upper limit on path lengths.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *        \c IGRAPH_EINVVID, invalid vertex id passed in
 *        \p vids.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 * Note that the time complexity is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa \ref igraph_betweenness() to calculate the exact betweenness and
 * \ref igraph_edge_betweenness_cutoff() to calculate the range-limited
 * edge betweenness.
 */
int igraph_betweenness_cutoff(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_vs_t vids, igraph_bool_t directed,
                              const igraph_vector_t *weights, igraph_real_t cutoff) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_inclist_t inclist;
    igraph_inclist_t fathers;
    long int source, j, neighbor;
    igraph_stack_t S;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    /* Note: nrgeo holds the number of shortest paths, which may be very large in some cases,
     * e.g. in a grid graph. If using an integer type, this results in overflow.
     * With a 'long long int', overflow already affects the result for a grid as small as 36*36.
     * Therefore, we use a 'double' instead. While a 'double' holds fewer digits than a 'long long int',
     * i.e. its precision is lower, it is effectively immune to overflow. The impact on the precision
     * of the final result is negligible. The max betweenness is correct to 14 decimal digits,
     * i.e. the precision limit of 'double', even for a 101*101 grid graph. */
    double *nrgeo = 0;
    double *tmpscore;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_vit_t vit;
    int cmp_result;
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length must agree with number of edges.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t minweight = igraph_vector_min(weights);
            if (minweight <= 0) {
                IGRAPH_ERROR("Weight vector must be positive.", IGRAPH_EINVAL);
            } else if (igraph_is_nan(minweight)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            } else if (minweight <= eps) {
                IGRAPH_WARNING("Some weights are smaller than epsilon, calculations may suffer from numerical precision.");
            }
        }
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_inclist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_inclist_destroy, &fathers);

    IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &S);
    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    nrgeo = igraph_Calloc(no_of_nodes, double);
    if (nrgeo == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
    tmpscore = igraph_Calloc(no_of_nodes, double);
    if (tmpscore == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    if (igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        igraph_vector_null(res);
        tmpres = res;
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
    }

    for (source = 0; source < no_of_nodes; source++) {
        if (weights) {
            igraph_i_sspf_weighted(graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, cutoff);
        } else {
            igraph_i_sspf (graph, source, &dist, nrgeo, &S, &fathers, &inclist, cutoff);
        }

        while (!igraph_stack_empty(&S)) {
            long int actnode = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *neis = igraph_inclist_get(&fathers, actnode);
            long nneis = igraph_vector_int_size(neis);
            for (j = 0; j < nneis; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                neighbor = IGRAPH_OTHER(graph, edge, actnode);
                tmpscore[neighbor] +=  (tmpscore[actnode] + 1) * nrgeo[neighbor] / nrgeo[actnode];
            }

            if (actnode != source) {
                VECTOR(*tmpres)[actnode] += tmpscore[actnode];
            }

            /* Reset variables */
            VECTOR(dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            tmpscore[actnode] = 0;
            igraph_vector_int_clear(igraph_inclist_get(&fathers, actnode));
        }

    } /* for source < no_of_nodes */

    if (!igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

        for (j = 0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), j++) {
            long int node = IGRAPH_VIT_GET(vit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[node];
        }

        no_of_nodes = (igraph_integer_t) j;

        igraph_vit_destroy(&vit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (!directed || !igraph_is_directed(graph)) {
        for (j = 0; j < no_of_nodes; j++) {
            VECTOR(*res)[j] /= 2.0;
        }
    }

    igraph_Free(nrgeo);
    igraph_Free(tmpscore);
    igraph_vector_destroy(&dist);
    igraph_stack_destroy(&S);
    igraph_inclist_destroy(&fathers);
    igraph_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(6);

    IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_betweenness_estimate
 * \brief Estimated betweenness centrality of some vertices.
 *
 * \deprecated-by igraph_betweenness_cutoff 0.9
 *
 * </para><para>
 * The betweenness centrality of a vertex is the number of geodesics
 * going through it. If there are more than one geodesic between two
 * vertices, the value of these geodesics are weighted by one over the
 * number of geodesics. When estimating betweenness centrality, igraph
 * takes into consideration only those paths that are shorter than or
 * equal to a prescribed length. Note that the estimated centrality
 * will always be less than the real one.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        estimated betweenness scores for the specified vertices.
 * \param vids The vertices of which the betweenness centrality scores
 *        will be estimated.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param cutoff The maximal length of paths that will be considered.
 *        If negative, the exact betweenness will be calculated, and
 *        there will be no upper limit on path lengths.
 * \param weights An optional vector containing edge weights for
 *        calculating weighted betweenness. No edge weight may be NaN.
 *        Supply a null pointer here for unweighted betweenness.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *        \c IGRAPH_EINVVID, invalid vertex id passed in
 *        \p vids.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 * Note that the time complexity is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph.
 */

int igraph_betweenness_estimate(const igraph_t *graph, igraph_vector_t *res,
                                const igraph_vs_t vids, igraph_bool_t directed,
                                igraph_real_t cutoff, const igraph_vector_t *weights) {
    IGRAPH_WARNING("igraph_betweenness_estimate is deprecated, use igraph_betweenness_cutoff.");
    return igraph_betweenness_cutoff(graph, res, vids, directed, weights, cutoff);
}

/***** Edge betweenness *****/


/**
 * \ingroup structural
 * \function igraph_edge_betweenness
 * \brief Betweenness centrality of the edges.
 *
 * </para><para>
 * The betweenness centrality of an edge is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the
 * number of geodesics.
 * \param graph The graph object.
 * \param result The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional weight vector for weighted edge
 *        betweenness. No edge weight may be NaN. Supply a null
 *        pointer here for the unweighted version.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph. See \ref igraph_edge_betweenness_cutoff() to
 *     compute the range-limited betweenness score of the edges in a graph.
 */
int igraph_edge_betweenness(const igraph_t *graph, igraph_vector_t *result,
                            igraph_bool_t directed,
                            const igraph_vector_t *weights) {
    return igraph_edge_betweenness_cutoff(graph, result, directed,
                                          weights, -1);
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness_cutoff
 * \brief Range-limited betweenness centrality of the edges.
 *
 * </para><para>
 * </para><para>
 * This function computes a range-limited version of edge betweenness centrality
 * by considering only those shortest paths whose length is no greater
 * then the given cutoff value.
 *
 * \param graph The graph object.
 * \param result The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional weight vector for weighted
 *        betweenness. No edge weight may be NaN. Supply a null
 *        pointer here for unweighted betweenness.
 * \param cutoff The maximal length of paths that will be considered.
 *        If negative, the exact betweenness will be calculated (no
 *        upper limit on path lengths).
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa \ref igraph_edge_betweenness() to compute the exact edge betweenness and
 * \ref igraph_betweenness_cutoff() to compute the range-limited vertex betweenness.
 */
int igraph_edge_betweenness_cutoff(const igraph_t *graph, igraph_vector_t *result,
                                   igraph_bool_t directed,
                                   const igraph_vector_t *weights, igraph_real_t cutoff) {
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_inclist_t inclist;
    igraph_inclist_t fathers;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    double *nrgeo;
    double *tmpscore;
    long int source, j;
    int cmp_result;
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;
    igraph_stack_t S;

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length must match number of edges.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t minweight = igraph_vector_min(weights);
            if (minweight <= 0) {
                IGRAPH_ERROR("Weight vector must be positive.", IGRAPH_EINVAL);
            } else if (igraph_is_nan(minweight)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            } else if (minweight <= eps) {
                IGRAPH_WARNING("Some weights are smaller than epsilon, calculations may suffer from numerical precision.");
            }
        }
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_inclist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_inclist_destroy, &fathers);
    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    nrgeo = igraph_Calloc(no_of_nodes, double);
    if (nrgeo == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
    tmpscore = igraph_Calloc(no_of_nodes, double);
    if (tmpscore == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &S);

    IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
    igraph_vector_null(result);

    for (source = 0; source < no_of_nodes; source++) {
        if (weights) {
            igraph_i_sspf_weighted (graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, cutoff);
        } else {
            igraph_i_sspf (graph, source, &dist, nrgeo, &S, &fathers, &inclist, cutoff);
        }

        while (!igraph_stack_empty(&S)) {
            long int w = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&fathers, w);
            long int fatv_len = igraph_vector_int_size(fatv);
            for (j = 0; j < fatv_len; j++) {
                long int fedge = (long int) VECTOR(*fatv)[j];
                long int neighbor = IGRAPH_OTHER(graph, fedge, w);
                tmpscore[neighbor] += (nrgeo[neighbor]) /
                                              nrgeo[w] * (1.0 + tmpscore[w]);
                VECTOR(*result)[fedge] +=
                    (tmpscore[w] + 1) * nrgeo[neighbor] /
                    nrgeo[w];
            }

            /* Reset variables */
            tmpscore[w] = 0;
            VECTOR(dist)[w] = 0;
            nrgeo[w] = 0;
            igraph_vector_int_clear(fatv);
        }

    } /* source < no_of_nodes */

    if (!directed || !igraph_is_directed(graph)) {
        for (j = 0; j < no_of_edges; j++) {
            VECTOR(*result)[j] /= 2.0;
        }
    }

    igraph_stack_destroy(&S);
    IGRAPH_FINALLY_CLEAN(2);

    igraph_inclist_destroy(&inclist);
    igraph_inclist_destroy(&fathers);
    igraph_vector_destroy(&dist);
    igraph_Free(tmpscore);
    igraph_Free(nrgeo);
    IGRAPH_FINALLY_CLEAN(5);

    IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness_estimate
 * \brief Estimated betweenness centrality of the edges.
 *
 * \deprecated-by igraph_edge_betweenness_cutoff 0.9
 *
 * </para><para>
 * The betweenness centrality of an edge is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the
 * number of geodesics. When estimating betweenness centrality, igraph
 * takes into consideration only those paths that are shorter than or
 * equal to a prescribed length. Note that the estimated centrality
 * will always be less than the real one.
 *
 * \param graph The graph object.
 * \param result The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param cutoff The maximal length of paths that will be considered.
 *        If negative, the exact betweenness will be calculated (no
 *        upper limit on path lengths).
 * \param weights An optional weight vector for weighted betweenness.
 *        No edge weight may be NaN. Supply a null pointer here for
 *        unweighted betweenness.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_betweenness() for calculating the betweenness score
 *     of the vertices in a graph.
 */
int igraph_edge_betweenness_estimate(const igraph_t *graph, igraph_vector_t *result,
                                   igraph_bool_t directed, igraph_real_t cutoff,
                                   const igraph_vector_t *weights) {
    IGRAPH_WARNING("igraph_edge_betweenness_estimate is deprecated, use igraph_edge_betweenness_cutoff.");
    return igraph_edge_betweenness_cutoff(graph, result, directed, weights, cutoff);
}

/* vertex subset*/

/**
 * \ingroup structural
 * \function igraph_betweenness_subset
 * \brief betweenness centrality for subset of vertices.
 *
 * </para><para>
 * This function computes the subset version of betweenness centrality
 * by considering only those shortest paths between vertices in a given
 * subset.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *         betweenness score for the subset of vertices.
 * \param vids The vertices for which the range-limited betweenness centrality
 *        scores will be computed.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional vector containing edge weights for
 *        calculating weighted betweenness. No edge weight may be NaN.
 *        Supply a null pointer here for unweighted betweenness.
 * \param sources A vector selector of the vertices which will be the sources
 *        of the shortest paths taken into considuration in the betweenness
 *        calculation.
 * \param targetes A vector selector of the vertices which will be the targets
 *        of the shortest paths taken into considuration in the betweenness
 *        calculation.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *        \c IGRAPH_EINVVID, invalid vertex id passed in
 *        \p vids.
 *
 * Time complexity: O(|V||E|),
 * |S| The number of vertices in the group
 * |E| The number of edges in the graph
 *
 * \sa \ref igraph_betweenness() to calculate the exact betweenness and
 * \ref igraph_edge_betweenness_cutoff() to calculate the range-limited
 * edge betweenness.
 */
int igraph_betweenness_subset(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_vs_t vids, igraph_bool_t directed,
                              const igraph_vs_t sources, const igraph_vs_t targets,
                              const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_inclist_t inclist, fathers;
    long int source, j;
    igraph_stack_t S;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    long int f;
    igraph_vector_t dist;
    double *nrgeo;
    double *tmpscore;
    igraph_vit_t vit, vit_source, vit_target;
    unsigned char *is_target;
    int cmp_result;
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t minweight = igraph_vector_min(weights);
            if (minweight <= 0) {
                IGRAPH_ERROR("Weight vector must be positive", IGRAPH_EINVAL);
            } else if (igraph_is_nan(minweight)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values", IGRAPH_EINVAL);
            } else if (minweight <= eps) {
                IGRAPH_WARNING("Some weights are smaller than epsilon, calculations may suffer from numerical precision.");
            }
        }
    }

    IGRAPH_CHECK(igraph_vit_create(graph, sources, &vit_source));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit_source);
    IGRAPH_CHECK(igraph_vit_create(graph, targets, &vit_target));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit_target); 

    if (!igraph_vs_is_all(&vids)) {
        /* subset */
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
    } else {
        /* only  */
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        igraph_vector_null(res);
        tmpres = res;
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_inclist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_inclist_destroy, &fathers);

    IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &S);
    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    nrgeo = igraph_Calloc(no_of_nodes, double);
    if (nrgeo == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
    tmpscore = igraph_Calloc(no_of_nodes, double);
    if (tmpscore == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    is_target = igraph_Calloc(no_of_nodes, unsigned char);
    if (is_target == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, is_target); 

    for (IGRAPH_VIT_RESET(vit_target); !IGRAPH_VIT_END(vit_target); IGRAPH_VIT_NEXT(vit_target)) {
        is_target[(long int) IGRAPH_VIT_GET(vit_target)] = 1;
    }    

    for (IGRAPH_VIT_RESET(vit_source); !IGRAPH_VIT_END(vit_source); IGRAPH_VIT_NEXT(vit_source)) {
        source = IGRAPH_VIT_GET(vit_source);

        if (weights) {
            igraph_i_sspf_weighted (graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, -1);
        } else {
            igraph_i_sspf (graph, source, &dist, nrgeo, &S, &fathers, &inclist, -1);
        }

        while (!igraph_stack_empty(&S)) {
            long int w = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&fathers, w);
            long int fatv_len = igraph_vector_int_size(fatv);
            for (j = 0; j < fatv_len; j++) {
                long int edge = (long int) VECTOR(*fatv)[j];
                f = IGRAPH_OTHER(graph, edge, w);
                if (is_target[w]){
                    tmpscore[f] += nrgeo[f] / nrgeo[w] * (1 + tmpscore[w]);
                }
                else{
                   tmpscore[f] += nrgeo[f] / nrgeo[w] * (tmpscore[w]);
                }
            }
            if (w != source) {
                VECTOR(*tmpres)[w] += tmpscore[w];
            }

            /* Reset variables */
            tmpscore[w] = 0;
            VECTOR(dist)[w] = 0;
            nrgeo[w] = 0;
            igraph_vector_int_clear(igraph_inclist_get(&fathers, w));
        }

    } 

    /* Keep only the requested vertices */
    if (!igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

        for (j = 0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), j++) {
            long int node = IGRAPH_VIT_GET(vit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[node];
        }

        igraph_vit_destroy(&vit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }

   if (!directed || !igraph_is_directed(graph)) {
        for (j = 0; j < no_of_nodes; j++) {
            VECTOR(*res)[j] /= 2.0;
        }
    }

    igraph_Free(nrgeo);
    igraph_Free(tmpscore);
    igraph_vector_destroy(&dist);
    igraph_stack_destroy(&S);
    igraph_inclist_destroy(&fathers);
    igraph_inclist_destroy(&inclist);
    igraph_Free(is_target);
    igraph_vit_destroy(&vit_source);
    igraph_vit_destroy(&vit_target);
    IGRAPH_FINALLY_CLEAN(9);

    IGRAPH_SUCCESS;
}

/***** Edge betweenness subset*****/

/**
 * \ingroup structural
 * \function igraph_edge_betweenness_subset
 * \brief Find the shortest path edge betweenness centrality for subset of vertices
 *
 * </para><para>
 * </para><para>
 * This function computes a subset version of edge betweenness centrality
 * by considering only those shortest paths betweenn two vertices in a subset.
 *
 * \param graph The graph object.
 * \param res The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param eids The edges for which the range-limited betweenness centrality
 *        scores will be computed.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional weight vector for weighted
 *        betweenness. No edge weight may be NaN. Supply a null
 *        pointer here for unweighted betweenness.
 * \param sources A vector selector of the vertices which will be the sources
 *        of the shortest paths taken into considuration in the betweenness
 *        calculation.
 * \param targetes A vector selector of the vertices which will be the targets
 *        of the shortest paths taken into considuration in the betweenness
 *        calculation.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa \ref igraph_edge_betweenness() to compute the exact edge betweenness and
 * \ref igraph_betweenness_cutoff() to compute the range-limited vertex betweenness.
 */
int igraph_edge_betweenness_subset(const igraph_t *graph, igraph_vector_t *res,
                                   const igraph_es_t eids, igraph_bool_t directed,
                                   const igraph_vs_t sources, const igraph_vs_t targets,
                                   const igraph_vector_t *weights) {
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_inclist_t inclist, fathers;
    igraph_vit_t vit_source, vit_target;
    igraph_eit_t eit;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    double *nrgeo;
    double *tmpscore;
    long int source, j;
    unsigned char *is_target;
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;
    igraph_stack_t S;

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length must match number of edges.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t minweight = igraph_vector_min(weights);
            if (minweight <= 0) {
                IGRAPH_ERROR("Weight vector must be positive.", IGRAPH_EINVAL);
            } else if (igraph_is_nan(minweight)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            } else if (minweight <= eps) {
                IGRAPH_WARNING("Some weights are smaller than epsilon, calculations may suffer from numerical precision.");
            }
        }
    }

    IGRAPH_CHECK(igraph_vit_create(graph, sources, &vit_source));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit_source);
    IGRAPH_CHECK(igraph_vit_create(graph, targets, &vit_target));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit_target); 

    if (!igraph_es_is_all(&eids)) {
        /* subset */
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_edges);
    } else {
        /* only  */
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_edges));
        igraph_vector_null(res);
        tmpres = res;
    }

    is_target = igraph_Calloc(no_of_nodes, unsigned char);
    if (is_target == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, is_target);

    for (IGRAPH_VIT_RESET(vit_target); !IGRAPH_VIT_END(vit_target); IGRAPH_VIT_NEXT(vit_target)) {
        is_target[(long int) IGRAPH_VIT_GET(vit_target)] = 1;
    }   

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_inclist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_inclist_destroy, &fathers);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    nrgeo = igraph_Calloc(no_of_nodes, double);
    if (nrgeo == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
    tmpscore = igraph_Calloc(no_of_nodes, double);
    if (tmpscore == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &S);

    for (IGRAPH_VIT_RESET(vit_source); !IGRAPH_VIT_END(vit_source); IGRAPH_VIT_NEXT(vit_source)) { 
        source = IGRAPH_VIT_GET(vit_source);

        if (weights) {
            igraph_i_sspf_weighted (graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, -1);
        } else {
            igraph_i_sspf (graph, source, &dist, nrgeo, &S, &fathers, &inclist, -1);
        }

        while (!igraph_stack_empty(&S)) {
            long int w = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&fathers, w);
            long int fatv_len = igraph_vector_int_size(fatv);
            for (j = 0; j < fatv_len; j++) {
                long int fedge = (long int) VECTOR(*fatv)[j];
                long int neighbor = IGRAPH_OTHER(graph, fedge, w);
                if (is_target[w]) {
                    tmpscore[neighbor] += (nrgeo[neighbor]) /
                                                nrgeo[w] * (1.0 + tmpscore[w]);
                    VECTOR(*tmpres)[fedge] +=
                        (tmpscore[w] + 1) * nrgeo[neighbor] /
                        nrgeo[w];
                }
                else{
                    tmpscore[neighbor] += (nrgeo[neighbor]) /
                                                nrgeo[w] * (tmpscore[w]);
                    VECTOR(*tmpres)[fedge] +=
                        (tmpscore[w]) * nrgeo[neighbor] /
                        nrgeo[w];
                }
            }

            /* Reset variables */
            tmpscore[w] = 0;
            VECTOR(dist)[w] = 0;
            nrgeo[w] = 0;
            igraph_vector_int_clear(fatv);
        }

    } /* source < no_of_nodes */

        /* Keep only the requested edges */
    if (!igraph_es_is_all(&eids)) {
        IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);
        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

        for (j = 0, IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit);
             IGRAPH_EIT_NEXT(eit), j++) {
            long int node = IGRAPH_EIT_GET(eit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[node];
        }

        igraph_eit_destroy(&eit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }


    if (!directed || !igraph_is_directed(graph)) {
        for (j = 0; j < no_of_edges; j++) {
            VECTOR(*res)[j] /= 2.0;
        }
    }

    igraph_stack_destroy(&S);
    igraph_inclist_destroy(&inclist);
    igraph_inclist_destroy(&fathers);
    igraph_vector_destroy(&dist);
    igraph_Free(tmpscore);
    igraph_Free(nrgeo);
    igraph_Free(is_target);
    igraph_vit_destroy(&vit_source);
    igraph_vit_destroy(&vit_target);
    IGRAPH_FINALLY_CLEAN(9);

    IGRAPH_SUCCESS;
}
