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

#include "igraph_memory.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_stack.h"
#include "igraph_dqueue.h"

#include "core/indheap.h"
#include "core/interruption.h"
#include "core/math.h"

/***** Vertex betweenness *****/

/**
 * \ingroup structural
 * \function igraph_betweenness
 * \brief Betweenness centrality of some vertices.
 *
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

static int igraph_i_betweenness_cutoff_weighted(
        const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vs_t vids,
        igraph_bool_t directed,
        igraph_real_t cutoff,
        const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_inclist_t inclist;
    igraph_adjlist_t fathers;
    long int source, j;
    igraph_stack_t S;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist, nrgeo, tmpscore;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_vit_t vit;
    int cmp_result;
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;

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

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_adjlist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &fathers);

    IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &S);
    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&tmpscore, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&nrgeo, no_of_nodes);

    if (igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        igraph_vector_null(res);
        tmpres = res;
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
    }

    for (source = 0; source < no_of_nodes; source++) {
        IGRAPH_PROGRESS("Betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        igraph_2wheap_push_with_index(&Q, source, -1.0);
        VECTOR(dist)[source] = 1.0;
        VECTOR(nrgeo)[source] = 1;

        while (!igraph_2wheap_empty(&Q)) {
            long int minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_delete_max(&Q);
            igraph_vector_int_t *neis;
            long int nlen;

            /* Ignore vertices that are more distant than the cutoff */
            if (cutoff >= 0 && mindist > cutoff + 1.0) {
                /* Reset variables if node is too distant */
                VECTOR(tmpscore)[minnei] = 0;
                VECTOR(dist)[minnei] = 0;
                VECTOR(nrgeo)[minnei] = 0;
                igraph_vector_int_clear(igraph_adjlist_get(&fathers, minnei));
                continue;
            }

            igraph_stack_push(&S, minnei);

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_inclist_get(&inclist, minnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                long int to = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_real_t curdist = VECTOR(dist)[to];

                if (curdist == 0) {
                    /* this means curdist is infinity */
                    cmp_result = -1;
                } else {
                    cmp_result = igraph_cmp_epsilon(altdist, curdist, eps);
                }

                if (curdist == 0) {
                    /* This is the first non-infinite distance */
                    igraph_vector_int_t *v = igraph_adjlist_get(&fathers, to);
                    igraph_vector_int_resize(v, 1);
                    VECTOR(*v)[0] = minnei;
                    VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];

                    VECTOR(dist)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
                } else if (cmp_result < 0) {
                    /* This is a shorter path */
                    igraph_vector_int_t *v = igraph_adjlist_get(&fathers, to);
                    igraph_vector_int_resize(v, 1);
                    VECTOR(*v)[0] = minnei;
                    VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];

                    VECTOR(dist)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
                } else if (cmp_result == 0 &&
                    (altdist <= cutoff + 1.0 || cutoff < 0)) {
                    /* Only add if the node is not more distant than the cutoff */
                    igraph_vector_int_t *v = igraph_adjlist_get(&fathers, to);
                    igraph_vector_int_push_back(v, minnei);
                    VECTOR(nrgeo)[to] += VECTOR(nrgeo)[minnei];
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

        while (!igraph_stack_empty(&S)) {
            long int w = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_adjlist_get(&fathers, w);
            long int fatv_len = igraph_vector_int_size(fatv);
            for (j = 0; j < fatv_len; j++) {
                long int f = (long int) VECTOR(*fatv)[j];
                VECTOR(tmpscore)[f] += VECTOR(nrgeo)[f] / VECTOR(nrgeo)[w] * (1 + VECTOR(tmpscore)[w]);
            }
            if (w != source) {
                VECTOR(*tmpres)[w] += VECTOR(tmpscore)[w];
            }

            /* Reset variables */
            VECTOR(tmpscore)[w] = 0;
            VECTOR(dist)[w] = 0;
            VECTOR(nrgeo)[w] = 0;
            igraph_vector_int_clear(igraph_adjlist_get(&fathers, w));
        }

    } /* source < no_of_nodes */

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

    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0, 0);

    igraph_vector_destroy(&nrgeo);
    igraph_vector_destroy(&tmpscore);
    igraph_vector_destroy(&dist);
    igraph_stack_destroy(&S);
    igraph_adjlist_destroy(&fathers);
    igraph_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(7);

    return 0;
}


/**
 * \ingroup structural
 * \function igraph_betweenness_cutoff
 * \brief Range-limited betweenness centrality.
 *
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

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    long int *distance;
    /* Note: nrgeo holds the number of shortest paths, which may be very large in some cases,
     * e.g. in a grid graph. If using an integer type, this results in overflow.
     * With a 'long long int', overflow already affects the result for a grid as small as 36*36.
     * Therefore, we use a 'double' instead. While a 'double' holds fewer digits than a 'long long int',
     * i.e. its precision is lower, it is effectively immune to overflow. The impact on the precision
     * of the final result is negligible. The max betweenness is correct to 14 decimal digits,
     * i.e. the precision limit of 'double', even for a 101*101 grid graph. */
    double *nrgeo = 0;
    double *tmpscore;
    igraph_stack_t stack = IGRAPH_STACK_NULL;
    long int source;
    long int j, k, nneis;
    igraph_vector_int_t *neis;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_vit_t vit;

    igraph_adjlist_t adjlist_out, adjlist_in;
    igraph_adjlist_t *adjlist_out_p, *adjlist_in_p;

    if (weights) {
        return igraph_i_betweenness_cutoff_weighted(graph, res, vids, directed,
                                                    cutoff, weights);
    }

    if (!igraph_vs_is_all(&vids)) {
        /* subset */
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
    } else {
        /* only  */
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        igraph_vector_null(res);
        tmpres = res;
    }

    directed = directed && igraph_is_directed(graph);
    if (directed) {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_in, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_in);
        adjlist_out_p = &adjlist_out;
        adjlist_in_p = &adjlist_in;
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_in, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_in);
        adjlist_out_p = &adjlist_out;
        adjlist_in_p = &adjlist_in;
    }
    for (j = 0; j < no_of_nodes; j++) {
        igraph_vector_int_clear(igraph_adjlist_get(adjlist_in_p, j));
    }

    distance = IGRAPH_CALLOC(no_of_nodes, long int);
    if (distance == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, distance);
    nrgeo = IGRAPH_CALLOC(no_of_nodes, double);
    if (nrgeo == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
    tmpscore = IGRAPH_CALLOC(no_of_nodes, double);
    if (tmpscore == 0) {
        IGRAPH_ERROR("Insufficient memory for betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    igraph_stack_init(&stack, no_of_nodes);
    IGRAPH_FINALLY(igraph_stack_destroy, &stack);

    /* here we go */

    for (source = 0; source < no_of_nodes; source++) {
        IGRAPH_PROGRESS("Betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_dqueue_push(&q, source));
        nrgeo[source] = 1;
        distance[source] = 1;

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);

            /* Ignore vertices that are more distant than the cutoff */
            if (cutoff >= 0 && distance[actnode] > cutoff + 1) {
                /* Reset variables if node is too distant */
                distance[actnode] = 0;
                nrgeo[actnode] = 0;
                tmpscore[actnode] = 0;
                igraph_vector_int_clear(igraph_adjlist_get(adjlist_in_p, actnode));
                continue;
            }

            IGRAPH_CHECK(igraph_stack_push(&stack, actnode));
            neis = igraph_adjlist_get(adjlist_out_p, actnode);
            nneis = igraph_vector_int_size(neis);
            for (j = 0; j < nneis; j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
                if (distance[neighbor] == 0) {
                    distance[neighbor] = distance[actnode] + 1;
                    IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                }
                if (distance[neighbor] == distance[actnode] + 1 &&
                    (distance[neighbor] <= cutoff + 1 || cutoff < 0)) {
                    /* Only add if the node is not more distant than the cutoff */
                    igraph_vector_int_t *v = igraph_adjlist_get(adjlist_in_p,
                                             neighbor);
                    igraph_vector_int_push_back(v, actnode);
                    nrgeo[neighbor] += nrgeo[actnode];
                }
            }
        } /* while !igraph_dqueue_empty */

        /* Ok, we've the distance of each node and also the number of
           shortest paths to them. Now we do an inverse search, starting
           with the farthest nodes. */
        while (!igraph_stack_empty(&stack)) {
            long int actnode = (long int) igraph_stack_pop(&stack);
            neis = igraph_adjlist_get(adjlist_in_p, actnode);
            nneis = igraph_vector_int_size(neis);
            for (j = 0; j < nneis; j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
                tmpscore[neighbor] +=  (tmpscore[actnode] + 1) * nrgeo[neighbor] / nrgeo[actnode];
            }

            if (actnode != source) {
                VECTOR(*tmpres)[actnode] += tmpscore[actnode];
            }

            /* Reset variables */
            distance[actnode] = 0;
            nrgeo[actnode] = 0;
            tmpscore[actnode] = 0;
            igraph_vector_int_clear(igraph_adjlist_get(adjlist_in_p, actnode));
        }

    } /* for source < no_of_nodes */

    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0, 0);

    /* clean  */
    IGRAPH_FREE(distance);
    IGRAPH_FREE(nrgeo);
    IGRAPH_FREE(tmpscore);

    igraph_dqueue_destroy(&q);
    igraph_stack_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(5);

    /* Keep only the requested vertices */
    if (!igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

        for (k = 0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), k++) {
            long int node = IGRAPH_VIT_GET(vit);
            VECTOR(*res)[k] = VECTOR(*tmpres)[node];
        }

        igraph_vit_destroy(&vit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }

    /* divide by 2 for undirected graph */
    if (!directed) {
        nneis = igraph_vector_size(res);
        for (j = 0; j < nneis; j++) {
            VECTOR(*res)[j] /= 2.0;
        }
    }

    igraph_adjlist_destroy(&adjlist_out);
    igraph_adjlist_destroy(&adjlist_in);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
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

static int igraph_i_edge_betweenness_cutoff_weighted(
        const igraph_t *graph,
        igraph_vector_t *result,
        igraph_bool_t directed,
        igraph_real_t cutoff,
        const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_inclist_t inclist;
    igraph_inclist_t fathers;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t distance, tmpscore;
    igraph_vector_long_t nrgeo;
    long int source, j;
    int cmp_result;
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;
    igraph_stack_t S;

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

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_inclist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_inclist_destroy, &fathers);

    IGRAPH_VECTOR_INIT_FINALLY(&distance, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&tmpscore, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_long_init(&nrgeo, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &nrgeo);

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &S);

    IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
    igraph_vector_null(result);

    for (source = 0; source < no_of_nodes; source++) {
        IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        /*     printf("source: %li\n", source); */

        igraph_2wheap_push_with_index(&Q, source, -1.0);
        VECTOR(distance)[source] = 1.0;
        VECTOR(nrgeo)[source] = 1;

        while (!igraph_2wheap_empty(&Q)) {
            long int minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_delete_max(&Q);
            igraph_vector_int_t *neis;
            long int nlen;

            /* printf("SP to %li is final, dist: %g, nrgeo: %li\n", minnei, */
            /* VECTOR(distance)[minnei]-1.0, VECTOR(nrgeo)[minnei]); */

            /* Ignore vertices that are more distant than the cutoff */
            if (cutoff >= 0 && VECTOR(distance)[minnei] > cutoff + 1.0) {
                /* Reset variables if node is too distant */
                VECTOR(tmpscore)[minnei] = 0;
                VECTOR(distance)[minnei] = 0;
                VECTOR(nrgeo)[minnei] = 0;
                igraph_vector_int_clear(igraph_inclist_get(&fathers, minnei));
                continue;
            }

            igraph_stack_push(&S, minnei);

            neis = igraph_inclist_get(&inclist, minnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                long int to = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_real_t curdist = VECTOR(distance)[to];

                if (curdist == 0) {
                    /* this means curdist is infinity */
                    cmp_result = -1;
                } else {
                    cmp_result = igraph_cmp_epsilon(altdist, curdist, eps);
                }

                /* printf("to=%ld, altdist = %lg, curdist = %lg, cmp = %d\n",
                  to, altdist, curdist-1, cmp_result); */
                if (curdist == 0) {
                    /* This is the first finite distance to 'to' */
                    igraph_vector_int_t *v = igraph_inclist_get(&fathers, to);
                    /* printf("Found first path to %li (from %li)\n", to, minnei); */
                    igraph_vector_int_resize(v, 1);
                    VECTOR(*v)[0] = edge;
                    VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];
                    VECTOR(distance)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
                } else if (cmp_result < 0) {
                    /* This is a shorter path */
                    igraph_vector_int_t *v = igraph_inclist_get(&fathers, to);
                    /* printf("Found a shorter path to %li (from %li)\n", to, minnei); */
                    igraph_vector_int_resize(v, 1);
                    VECTOR(*v)[0] = edge;
                    VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];
                    VECTOR(distance)[to] = altdist;
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
                } else if (cmp_result == 0 &&
                    (altdist <= cutoff + 1.0 || cutoff < 0)) {
                    /* Only add if the edge is not more distant than the cutoff */
                    igraph_vector_int_t *v = igraph_inclist_get(&fathers, to);
                    /* printf("Found a second SP to %li (from %li)\n", to, minnei); */
                    IGRAPH_CHECK(igraph_vector_int_push_back(v, edge));
                    VECTOR(nrgeo)[to] += VECTOR(nrgeo)[minnei];
                }
            }

        } /* igraph_2wheap_empty(&Q) */

        while (!igraph_stack_empty(&S)) {
            long int w = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&fathers, w);
            long int fatv_len = igraph_vector_int_size(fatv);
            /* printf("Popping %li.\n", w); */
            for (j = 0; j < fatv_len; j++) {
                long int fedge = (long int) VECTOR(*fatv)[j];
                long int neighbor = IGRAPH_OTHER(graph, fedge, w);
                VECTOR(tmpscore)[neighbor] += ((double)VECTOR(nrgeo)[neighbor]) /
                                              VECTOR(nrgeo)[w] * (1.0 + VECTOR(tmpscore)[w]);
                /* printf("Scoring %li (edge %li)\n", neighbor, fedge); */
                VECTOR(*result)[fedge] +=
                    ((VECTOR(tmpscore)[w] + 1) * VECTOR(nrgeo)[neighbor]) /
                    VECTOR(nrgeo)[w];
            }

            /* Reset variables */
            VECTOR(tmpscore)[w] = 0;
            VECTOR(distance)[w] = 0;
            VECTOR(nrgeo)[w] = 0;
            igraph_vector_int_clear(fatv);
        }

    } /* source < no_of_nodes */

    if (!directed || !igraph_is_directed(graph)) {
        for (j = 0; j < no_of_edges; j++) {
            VECTOR(*result)[j] /= 2.0;
        }
    }

    IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0, 0);

    igraph_stack_destroy(&S);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    igraph_inclist_destroy(&inclist);
    igraph_inclist_destroy(&fathers);
    igraph_vector_destroy(&distance);
    igraph_vector_destroy(&tmpscore);
    igraph_vector_long_destroy(&nrgeo);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness
 * \brief Betweenness centrality of the edges.
 *
 * The betweenness centrality of an edge is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the
 * number of geodesics.
 *
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
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    long int *distance;
    double *nrgeo;
    double *tmpscore;
    igraph_stack_t stack = IGRAPH_STACK_NULL;
    long int source;
    long int j;

    igraph_inclist_t elist_out, elist_in;
    igraph_inclist_t *elist_out_p, *elist_in_p;
    igraph_vector_int_t *neip;
    long int neino;
    long int i;

    if (weights) {
        return igraph_i_edge_betweenness_cutoff_weighted(graph, result,
                                                         directed, cutoff, weights);
    }

    directed = directed && igraph_is_directed(graph);
    if (directed) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_out, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_in, IGRAPH_IN, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_in);
        elist_out_p = &elist_out;
        elist_in_p = &elist_in;
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_out, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
        elist_out_p = elist_in_p = &elist_out;
    }

    distance = IGRAPH_CALLOC(no_of_nodes, long int);
    if (distance == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, distance);
    nrgeo = IGRAPH_CALLOC(no_of_nodes, double);
    if (nrgeo == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
    tmpscore = IGRAPH_CALLOC(no_of_nodes, double);
    if (tmpscore == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_stack_init(&stack, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &stack);

    IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));

    igraph_vector_null(result);

    /* here we go */

    for (source = 0; source < no_of_nodes; source++) {
        IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_dqueue_push(&q, source));

        nrgeo[source] = 1;
        distance[source] = 0;

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);

            if (cutoff >= 0 && distance[actnode] > cutoff ) {
                /* Reset variables if node is too distant */
                distance[actnode] = 0;
                tmpscore[actnode] = 0;
                nrgeo[actnode] = 0;
                continue;
            }

            IGRAPH_CHECK(igraph_stack_push(&stack, actnode));

            /* check the neighbors and add to them to the queue if unseen before */
            neip = igraph_inclist_get(elist_out_p, actnode);
            neino = igraph_vector_int_size(neip);
            for (i = 0; i < neino; i++) {
                igraph_integer_t edge = (igraph_integer_t) VECTOR(*neip)[i];
                long int neighbor = (long int) IGRAPH_OTHER(graph, edge, actnode);
                if (nrgeo[neighbor] != 0) {
                    /* we've already seen this node, another shortest path? */
                    if (distance[neighbor] == distance[actnode] + 1) {
                        nrgeo[neighbor] += nrgeo[actnode];
                    }
                } else if (distance[actnode] + 1 <= cutoff || cutoff < 0) {
                    /* we haven't seen this node yet, but we only consider
                     * it if it is not more distant than the cutoff. */
                    nrgeo[neighbor] += nrgeo[actnode];
                    distance[neighbor] = distance[actnode] + 1;
                    IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                }
            }
        } /* while !igraph_dqueue_empty */

        /* Ok, we've the distance of each node and also the number of
           shortest paths to them. Now we do an inverse search, starting
           with the farthest nodes. */
        while (!igraph_stack_empty(&stack)) {
            long int actnode = (long int) igraph_stack_pop(&stack);
            if (distance[actnode] < 1) {
                distance[actnode] = 0;
                tmpscore[actnode] = 0;
                nrgeo[actnode] = 0;
                continue;    /* skip source node */
            }
            /* set the temporary score of the friends */
            neip = igraph_inclist_get(elist_in_p, actnode);
            neino = igraph_vector_int_size(neip);
            for (i = 0; i < neino; i++) {
                igraph_integer_t edgeno = (igraph_integer_t) VECTOR(*neip)[i];
                long int neighbor = (long int) IGRAPH_OTHER(graph, edgeno, actnode);
                if (distance[neighbor] == distance[actnode] - 1 &&
                    nrgeo[neighbor] != 0) {
                    tmpscore[neighbor] +=
                        (tmpscore[actnode] + 1) * nrgeo[neighbor] / nrgeo[actnode];
                    VECTOR(*result)[edgeno] +=
                        (tmpscore[actnode] + 1) * nrgeo[neighbor] / nrgeo[actnode];
                }
            }
            /* Reset variables */
            distance[actnode] = 0;
            tmpscore[actnode] = 0;
            nrgeo[actnode] = 0;
        }
        /* Ok, we've the scores for this source */
    } /* for source <= no_of_nodes */
    IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0, 0);

    /* clean and return */
    IGRAPH_FREE(distance);
    IGRAPH_FREE(nrgeo);
    IGRAPH_FREE(tmpscore);
    igraph_dqueue_destroy(&q);
    igraph_stack_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(5);

    if (directed) {
        igraph_inclist_destroy(&elist_out);
        igraph_inclist_destroy(&elist_in);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        igraph_inclist_destroy(&elist_out);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* divide by 2 for undirected graph */
    if (!directed || !igraph_is_directed(graph)) {
        for (j = 0; j < igraph_vector_size(result); j++) {
            VECTOR(*result)[j] /= 2.0;
        }
    }

    return 0;
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
