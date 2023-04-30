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

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_nongraph.h"
#include "igraph_progress.h"
#include "igraph_stack.h"

#include "core/indheap.h"
#include "core/interruption.h"

/*
 * We provide separate implementations of single-source shortest path searches,
 * one with incidence lists and one with adjacency lists. We use the implementation
 * based on adjacency lists when possible (i.e. when weights are not needed) to
 * avoid an expensive IGRAPH_OTHER() lookup on edge IDs. The cost of this macro
 * comes from the inability of the branch predictor to predict accurately whether
 * the condition in the macro will be true or not.
 *
 * The following four functions are very similar in their structure. If you make
 * a modification to one of them, consider whether the same modification makes
 * sense in the context of the remaining three functions as well.
 */

/**
 * Internal function to calculate the single source shortest paths for the
 * vertex unweighted case.
 *
 * \param  graph   the graph to calculate the single source shortest paths on
 * \param  source  the source node
 * \param  dist    distance of each node from the source node \em plus one;
 *                 must be filled with zeros initially
 * \param  nrgeo   vector storing the number of geodesics from the source node
 *                 to each node; must be filled with zeros initially
 * \param  stack   stack in which the nodes are pushed in the order they are
 *                 discovered during the traversal
 * \param  parents adjacent list that starts empty and that stores the IDs
 *                 of the vertices that lead to a given node during the traversal
 * \param  adjlist the adjacency list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static igraph_error_t igraph_i_sspf(
        igraph_integer_t source,
        igraph_vector_t *dist,
        igraph_real_t *nrgeo,
        igraph_stack_int_t *stack,
        igraph_adjlist_t *parents,
        const igraph_adjlist_t *adjlist,
        igraph_real_t cutoff) {

    igraph_dqueue_int_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    igraph_integer_t nlen;

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&queue, 100);

    IGRAPH_CHECK(igraph_dqueue_int_push(&queue, source));
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_dqueue_int_empty(&queue)) {
        igraph_integer_t actnode = igraph_dqueue_int_pop(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && VECTOR(*dist)[actnode] > cutoff + 1) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            igraph_vector_int_clear(igraph_adjlist_get(parents, actnode));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_int_push(stack, actnode));

        /* Examine the neighbors of this node */
        neis = igraph_adjlist_get(adjlist, actnode);
        nlen = igraph_vector_int_size(neis);
        for (igraph_integer_t j = 0; j < nlen; j++) {
            igraph_integer_t neighbor = VECTOR(*neis)[j];

            if (VECTOR(*dist)[neighbor] == 0) {
                /* We have found 'neighbor' for the first time */
                VECTOR(*dist)[neighbor] = VECTOR(*dist)[actnode] + 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(&queue, neighbor));
            }

            if (VECTOR(*dist)[neighbor] == VECTOR(*dist)[actnode] + 1 &&
                (VECTOR(*dist)[neighbor] <= cutoff + 1 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_adjlist_get(parents, neighbor);
                IGRAPH_CHECK(igraph_vector_int_push_back(v, actnode));
                nrgeo[neighbor] += nrgeo[actnode];
            }
        }
    }

    igraph_dqueue_int_destroy(&queue);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * Internal function to calculate the single source shortest paths for the
 * edge unweighted case.
 *
 * \param  graph   the graph to calculate the single source shortest paths on
 * \param  source  the source node
 * \param  dist    distance of each node from the source node \em plus one;
 *                 must be filled with zeros initially
 * \param  nrgeo   vector storing the number of geodesics from the source node
 *                 to each node; must be filled with zeros initially
 * \param  stack   stack in which the nodes are pushed in the order they are
 *                 discovered during the traversal
 * \param  parents incidence list that starts empty and that stores the IDs
 *                 of the edges that lead to a given node during the traversal
 * \param  inclist the incidence list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static igraph_error_t igraph_i_sspf_edge(
        const igraph_t *graph,
        igraph_integer_t source,
        igraph_vector_t *dist,
        igraph_real_t *nrgeo,
        igraph_stack_int_t *stack,
        igraph_inclist_t *parents,
        const igraph_inclist_t *inclist,
        igraph_real_t cutoff) {

    igraph_dqueue_int_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    igraph_integer_t nlen;

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&queue, 100);

    IGRAPH_CHECK(igraph_dqueue_int_push(&queue, source));
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_dqueue_int_empty(&queue)) {
        igraph_integer_t actnode = igraph_dqueue_int_pop(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && VECTOR(*dist)[actnode] > cutoff + 1) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            igraph_vector_int_clear(igraph_inclist_get(parents, actnode));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_int_push(stack, actnode));

        /* Examine the neighbors of this node */
        neis = igraph_inclist_get(inclist, actnode);
        nlen = igraph_vector_int_size(neis);
        for (igraph_integer_t j = 0; j < nlen; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t neighbor = IGRAPH_OTHER(graph, edge, actnode);

            if (VECTOR(*dist)[neighbor] == 0) {
                /* We have found 'neighbor' for the first time */
                VECTOR(*dist)[neighbor] = VECTOR(*dist)[actnode] + 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(&queue, neighbor));
            }

            if (VECTOR(*dist)[neighbor] == VECTOR(*dist)[actnode] + 1 &&
                (VECTOR(*dist)[neighbor] <= cutoff + 1 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_inclist_get(parents, neighbor);
                IGRAPH_CHECK(igraph_vector_int_push_back(v, edge));
                nrgeo[neighbor] += nrgeo[actnode];
            }
        }
    }

    igraph_dqueue_int_destroy(&queue);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * Internal function to calculate the single source shortest paths for the vertex
 * weighted case.
 *
 * \param  graph   the graph to calculate the single source shortest paths on
 * \param  weights the weights of the edges
 * \param  source  the source node
 * \param  dist    distance of each node from the source node \em plus one;
 *                 must be filled with zeros initially
 * \param  nrgeo   vector storing the number of geodesics from the source node
 *                 to each node; must be filled with zeros initially
 * \param  stack   stack in which the nodes are pushed in the order they are
 *                 discovered during the traversal
 * \param  parents adjacency list that starts empty and that stores the IDs
 *                 of the vertices that lead to a given node during the traversal
 * \param  inclist the incidence list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static igraph_error_t igraph_i_sspf_weighted(
        const igraph_t *graph,
        igraph_integer_t source,
        igraph_vector_t *dist,
        igraph_real_t *nrgeo,
        const igraph_vector_t *weights,
        igraph_stack_int_t *stack,
        igraph_adjlist_t *parents,
        const igraph_inclist_t *inclist,
        igraph_real_t cutoff) {

    const igraph_real_t eps = IGRAPH_SHORTEST_PATH_EPSILON;

    int cmp_result;
    igraph_2wheap_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    igraph_integer_t nlen;

    /* TODO: this is an O|V| step here. We could save some time by pre-allocating
     * the two-way heap in the caller and re-using it here */
    IGRAPH_CHECK(igraph_2wheap_init(&queue, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &queue);

    igraph_2wheap_push_with_index(&queue, source, -1.0);
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_2wheap_empty(&queue)) {
        igraph_integer_t minnei = igraph_2wheap_max_index(&queue);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && mindist > cutoff + 1.0) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[minnei] = 0;
            nrgeo[minnei] = 0;
            igraph_vector_int_clear(igraph_adjlist_get(parents, minnei));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_int_push(stack, minnei));

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_inclist_get(inclist, minnei);
        nlen = igraph_vector_int_size(neis);
        for (igraph_integer_t j = 0; j < nlen; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t to = IGRAPH_OTHER(graph, edge, minnei);
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
                v = igraph_adjlist_get(parents, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = minnei;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&queue, to, -altdist));
            } else if (cmp_result < 0) {
                /* This is a shorter path */
                v = igraph_adjlist_get(parents, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = minnei;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                igraph_2wheap_modify(&queue, to, -altdist);
            } else if (cmp_result == 0 && (altdist <= cutoff + 1.0 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_adjlist_get(parents, to);
                IGRAPH_CHECK(igraph_vector_int_push_back(v, minnei));
                nrgeo[to] += nrgeo[minnei];
            }
        }
    }

    igraph_2wheap_destroy(&queue);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * Internal function to calculate the single source shortest paths for the edge
 * weighted case.
 *
 * \param  graph   the graph to calculate the single source shortest paths on
 * \param  weights the weights of the edges
 * \param  source  the source node
 * \param  dist    distance of each node from the source node \em plus one;
 *                 must be filled with zeros initially
 * \param  nrgeo   vector storing the number of geodesics from the source node
 *                 to each node; must be filled with zeros initially
 * \param  stack   stack in which the nodes are pushed in the order they are
 *                 discovered during the traversal
 * \param  parents incidence list that starts empty and that stores the IDs
 *                 of the edges that lead to a given node during the traversal
 * \param  inclist the incidence list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static igraph_error_t igraph_i_sspf_weighted_edge(
        const igraph_t *graph,
        igraph_integer_t source,
        igraph_vector_t *dist,
        igraph_real_t *nrgeo,
        const igraph_vector_t *weights,
        igraph_stack_int_t *stack,
        igraph_inclist_t *parents,
        const igraph_inclist_t *inclist,
        igraph_real_t cutoff) {

    const igraph_real_t eps = IGRAPH_SHORTEST_PATH_EPSILON;

    int cmp_result;
    igraph_2wheap_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    igraph_integer_t nlen;

    /* TODO: this is an O|V| step here. We could save some time by pre-allocating
     * the two-way heap in the caller and re-using it here */
    IGRAPH_CHECK(igraph_2wheap_init(&queue, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &queue);

    igraph_2wheap_push_with_index(&queue, source, -1.0);
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_2wheap_empty(&queue)) {
        igraph_integer_t minnei = igraph_2wheap_max_index(&queue);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && mindist > cutoff + 1.0) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[minnei] = 0;
            nrgeo[minnei] = 0;
            igraph_vector_int_clear(igraph_inclist_get(parents, minnei));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_int_push(stack, minnei));

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_inclist_get(inclist, minnei);
        nlen = igraph_vector_int_size(neis);
        for (igraph_integer_t j = 0; j < nlen; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t to = IGRAPH_OTHER(graph, edge, minnei);
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
                v = igraph_inclist_get(parents, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = edge;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&queue, to, -altdist));
            } else if (cmp_result < 0) {
                /* This is a shorter path */
                v = igraph_inclist_get(parents, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = edge;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                igraph_2wheap_modify(&queue, to, -altdist);
            } else if (cmp_result == 0 && (altdist <= cutoff + 1.0 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_inclist_get(parents, to);
                IGRAPH_CHECK(igraph_vector_int_push_back(v, edge));
                nrgeo[to] += nrgeo[minnei];
            }
        }
    }

    igraph_2wheap_destroy(&queue);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_betweenness_check_weights(
    const igraph_vector_t* weights, igraph_integer_t no_of_edges
) {
    igraph_real_t minweight;

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length must match the number of edges.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            minweight = igraph_vector_min(weights);
            if (minweight <= 0) {
                IGRAPH_ERROR("Weight vector must be positive.", IGRAPH_EINVAL);
            } else if (isnan(minweight)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            } else if (minweight <= IGRAPH_SHORTEST_PATH_EPSILON) {
                IGRAPH_WARNING(
                    "Some weights are smaller than epsilon, calculations may "
                    "suffer from numerical precision issues."
                );
            }
        }
    }

    return IGRAPH_SUCCESS;
}

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
 *
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
 *        \c IGRAPH_EINVVID, invalid vertex ID passed in
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
igraph_error_t igraph_betweenness(const igraph_t *graph, igraph_vector_t *res,
                       const igraph_vs_t vids, igraph_bool_t directed,
                       const igraph_vector_t* weights) {
    return igraph_betweenness_cutoff(graph, res, vids, directed, weights, -1);
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
 *        \c IGRAPH_EINVVID, invalid vertex ID passed in
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
igraph_error_t igraph_betweenness_cutoff(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_vs_t vids, igraph_bool_t directed,
                              const igraph_vector_t *weights, igraph_real_t cutoff) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_adjlist_t adjlist, parents;
    igraph_inclist_t inclist;
    igraph_integer_t source, j, neighbor;
    igraph_stack_int_t S;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    /* Note: nrgeo holds the number of shortest paths, which may be very large in some cases,
     * e.g. in a grid graph. If using an integer type, this results in overflow.
     * With a 'long long int', overflow already affects the result for a grid as small as 36*36.
     * Therefore, we use a 'igraph_real_t' instead. While a 'igraph_real_t' holds fewer digits than a
     * 'long long int', i.e. its precision is lower, it is effectively immune to overflow.
     * The impact on the precision of the final result is negligible. The max betweenness
     * is correct to 14 decimal digits, i.e. the precision limit of 'igraph_real_t', even
     * for a 101*101 grid graph. */
    igraph_real_t *nrgeo = 0;
    igraph_real_t *tmpscore;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    }

    IGRAPH_CHECK(igraph_adjlist_init_empty(&parents, no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &parents);

    IGRAPH_CHECK(igraph_stack_int_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &S);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    nrgeo = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(nrgeo, "Insufficient memory for betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, nrgeo);

    tmpscore = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(tmpscore, "Insufficient memory for betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, tmpscore);

    if (igraph_vs_is_all(&vids)) {
        /* result covers all vertices */
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        igraph_vector_null(res);
        tmpres = res;
    } else {
        /* result needed only for a subset of the vertices */
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
    }

    for (source = 0; source < no_of_nodes; source++) {

        /* Loop invariant that is valid at this point:
         *
         * - the stack S is empty
         * - the 'dist' vector contains zeros only
         * - the 'nrgeo' array contains zeros only
         * - the 'tmpscore' array contains zeros only
         * - the 'parents' adjacency list contains empty vectors only
         */

        IGRAPH_PROGRESS("Betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted(graph, source, &dist, nrgeo, weights, &S, &parents, &inclist, cutoff));
        } else {
            IGRAPH_CHECK(igraph_i_sspf(source, &dist, nrgeo, &S, &parents, &adjlist, cutoff));
        }

        /* Aggregate betweenness scores for the nodes we have reached in this
         * traversal */
        while (!igraph_stack_int_empty(&S)) {
            igraph_integer_t actnode = igraph_stack_int_pop(&S);
            igraph_vector_int_t *neis = igraph_adjlist_get(&parents, actnode);
            igraph_integer_t nneis = igraph_vector_int_size(neis);
            igraph_real_t coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];

            for (j = 0; j < nneis; j++) {
                neighbor = VECTOR(*neis)[j];
                tmpscore[neighbor] += nrgeo[neighbor] * coeff;
            }

            if (actnode != source) {
                VECTOR(*tmpres)[actnode] += tmpscore[actnode];
            }

            /* Reset variables to ensure that the 'for' loop invariant will
             * still be valid in the next iteration */

            VECTOR(dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            tmpscore[actnode] = 0;
            igraph_vector_int_clear(neis);
        }

    } /* for source < no_of_nodes */

    /* Keep only the requested vertices */
    if (!igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

        for (j = 0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), j++) {
            igraph_integer_t node = IGRAPH_VIT_GET(vit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[node];
        }

        igraph_vit_destroy(&vit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (!directed || !igraph_is_directed(graph)) {
        igraph_vector_scale(res, 0.5);
    }

    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0, 0);

    IGRAPH_FREE(nrgeo);
    IGRAPH_FREE(tmpscore);
    igraph_vector_destroy(&dist);
    igraph_stack_int_destroy(&S);
    igraph_adjlist_destroy(&parents);
    if (weights) {
        igraph_inclist_destroy(&inclist);
    } else {
        igraph_adjlist_destroy(&adjlist);
    }
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/***** Edge betweenness *****/


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
igraph_error_t igraph_edge_betweenness(const igraph_t *graph, igraph_vector_t *result,
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
igraph_error_t igraph_edge_betweenness_cutoff(const igraph_t *graph, igraph_vector_t *result,
                                   igraph_bool_t directed,
                                   const igraph_vector_t *weights, igraph_real_t cutoff) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_inclist_t inclist, parents;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    igraph_real_t *nrgeo;
    igraph_real_t *tmpscore;
    igraph_integer_t source, j;
    igraph_stack_int_t S;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_inclist_init_empty(&parents, no_of_nodes));
    IGRAPH_FINALLY(igraph_inclist_destroy, &parents);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    nrgeo = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(nrgeo, "Insufficient memory for edge betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, nrgeo);

    tmpscore = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    if (tmpscore == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    IGRAPH_CHECK(igraph_stack_int_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &S);

    IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
    igraph_vector_null(result);

    for (source = 0; source < no_of_nodes; source++) {

        /* Loop invariant that is valid at this point:
         *
         * - the stack S is empty
         * - the 'dist' vector contains zeros only
         * - the 'nrgeo' array contains zeros only
         * - the 'tmpscore' array contains zeros only
         * - the 'parents' incidence list contains empty vectors only
         */

        IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted_edge(graph, source, &dist, nrgeo, weights, &S, &parents, &inclist, cutoff));
        } else {
            IGRAPH_CHECK(igraph_i_sspf_edge(graph, source, &dist, nrgeo, &S, &parents, &inclist, cutoff));
        }

        /* Aggregate betweenness scores for the edges we have reached in this
         * traversal */
        while (!igraph_stack_int_empty(&S)) {
            igraph_integer_t actnode = igraph_stack_int_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&parents, actnode);
            igraph_integer_t fatv_len = igraph_vector_int_size(fatv);
            igraph_real_t coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];

            for (j = 0; j < fatv_len; j++) {
                igraph_integer_t fedge = VECTOR(*fatv)[j];
                igraph_integer_t neighbor = IGRAPH_OTHER(graph, fedge, actnode);
                tmpscore[neighbor] += nrgeo[neighbor] * coeff;
                VECTOR(*result)[fedge] += nrgeo[neighbor] * coeff;
            }

            /* Reset variables to ensure that the 'for' loop invariant will
             * still be valid in the next iteration */

            VECTOR(dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            tmpscore[actnode] = 0;
            igraph_vector_int_clear(fatv);
        }
    } /* source < no_of_nodes */

    if (!directed || !igraph_is_directed(graph)) {
        igraph_vector_scale(result, 0.5);
    }

    IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0, 0);

    igraph_stack_int_destroy(&S);
    igraph_inclist_destroy(&inclist);
    igraph_inclist_destroy(&parents);
    igraph_vector_destroy(&dist);
    IGRAPH_FREE(tmpscore);
    IGRAPH_FREE(nrgeo);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_betweenness_subset
 * \brief Betweenness centrality for a subset of source and target vertices.
 *
 * This function computes the subset-limited version of betweenness centrality
 * by considering only those shortest paths that lie between vertices in a given
 * source and target subset.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *         betweenness score for the subset of vertices.
 * \param vids The vertices for which the subset-limited betweenness centrality
 *        scores will be computed.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional vector containing edge weights for
 *        calculating weighted betweenness. No edge weight may be NaN.
 *        Supply a null pointer here for unweighted betweenness.
 * \param sources A vertex selector for the sources of the shortest paths taken
 *        into considuration in the betweenness calculation.
 * \param targets A vertex selector for the targets of the shortest paths taken
 *        into considuration in the betweenness calculation.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for temporary data.
 *        \c IGRAPH_EINVVID, invalid vertex ID passed in \p vids,
 *        \p sources or \p targets
 *
 * Time complexity: O(|S||E|), where
 * |S| is the number of vertices in the subset and
 * |E| is the number of edges in the graph.
 *
 * \sa \ref igraph_betweenness() to calculate the exact vertex betweenness and
 * \ref igraph_betweenness_cutoff() to calculate the range-limited vertex
 * betweenness.
 */
igraph_error_t igraph_betweenness_subset(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_vs_t vids, igraph_bool_t directed,
                              const igraph_vs_t sources, const igraph_vs_t targets,
                              const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_sources;
    igraph_integer_t no_of_processed_sources;
    igraph_adjlist_t adjlist, parents;
    igraph_inclist_t inclist;
    igraph_integer_t source, j;
    igraph_stack_int_t S;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_integer_t father;
    igraph_vector_t dist;
    igraph_real_t *nrgeo;
    igraph_real_t *tmpscore;
    igraph_vit_t vit;
    bool *is_target;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    IGRAPH_CHECK(igraph_vs_size(graph, &sources, &no_of_sources));

    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    }

    IGRAPH_CHECK(igraph_adjlist_init_empty(&parents, no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &parents);

    IGRAPH_CHECK(igraph_stack_int_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &S);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    nrgeo = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(nrgeo, "Insufficient memory for subset betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, nrgeo);

    tmpscore = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(tmpscore, "Insufficient memory for subset betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, tmpscore);

    is_target = IGRAPH_CALLOC(no_of_nodes, bool);
    IGRAPH_CHECK_OOM(is_target, "Insufficient memory for subset betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, is_target);

    IGRAPH_CHECK(igraph_vit_create(graph, targets, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        is_target[IGRAPH_VIT_GET(vit)] = true;
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    if (!igraph_vs_is_all(&vids)) {
        /* result needed only for a subset of the vertices */
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
    } else {
        /* result covers all vertices */
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        igraph_vector_null(res);
        tmpres = res;
    }

    IGRAPH_CHECK(igraph_vit_create(graph, sources, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    for (
        no_of_processed_sources = 0, IGRAPH_VIT_RESET(vit);
        !IGRAPH_VIT_END(vit);
        IGRAPH_VIT_NEXT(vit), no_of_processed_sources++
    ) {
        source = IGRAPH_VIT_GET(vit);

        IGRAPH_PROGRESS(
            "Betweenness centrality (subset): ",
            100.0 * no_of_processed_sources / no_of_sources, 0
        );
        IGRAPH_ALLOW_INTERRUPTION();

        /* Loop invariant that is valid at this point:
         *
         * - the stack S is empty
         * - the 'dist' vector contains zeros only
         * - the 'nrgeo' array contains zeros only
         * - the 'tmpscore' array contains zeros only
         * - the 'parents' adjacency list contains empty vectors only
         */

        /* TODO: there is more room for optimization here; the single-source
         * shortest path search runs until it reaches all the nodes in the
         * component of the source node even if we are only interested in a
         * smaller target subset. We could stop the search when all target
         * nodes were reached.
         */

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted(graph, source, &dist, nrgeo, weights, &S, &parents, &inclist, -1));
        } else {
            IGRAPH_CHECK(igraph_i_sspf(source, &dist, nrgeo, &S, &parents, &adjlist, -1));
        }

        /* Aggregate betweenness scores for the nodes we have reached in this
         * traversal */
        while (!igraph_stack_int_empty(&S)) {
            igraph_integer_t actnode = igraph_stack_int_pop(&S);
            igraph_vector_int_t *fatv = igraph_adjlist_get(&parents, actnode);
            igraph_integer_t fatv_len = igraph_vector_int_size(fatv);
            igraph_real_t coeff;

            if (is_target[actnode]) {
                coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];
            } else {
                coeff = tmpscore[actnode] / nrgeo[actnode];
            }

            for (j = 0; j < fatv_len; j++) {
                father = VECTOR(*fatv)[j];
                tmpscore[father] += nrgeo[father] * coeff;
            }

            if (actnode != source) {
                VECTOR(*tmpres)[actnode] += tmpscore[actnode];
            }

            /* Reset variables to ensure that the 'for' loop invariant will
             * still be valid in the next iteration */

            VECTOR(dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            tmpscore[actnode] = 0;
            igraph_vector_int_clear(fatv);
        }
    }

    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    /* Keep only the requested vertices */
    if (!igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);

        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));
        for (j = 0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), j++) {
            igraph_integer_t node = IGRAPH_VIT_GET(vit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[node];
        }

        igraph_vit_destroy(&vit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }

   if (!directed || !igraph_is_directed(graph)) {
        igraph_vector_scale(res, 0.5);
    }

    IGRAPH_FREE(is_target);
    IGRAPH_FREE(tmpscore);
    IGRAPH_FREE(nrgeo);
    igraph_vector_destroy(&dist);
    igraph_stack_int_destroy(&S);
    igraph_adjlist_destroy(&parents);
    if (weights) {
        igraph_inclist_destroy(&inclist);
    } else {
        igraph_adjlist_destroy(&adjlist);
    }
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness_subset
 * \brief Edge betweenness centrality for a subset of source and target vertices.
 *
 * This function computes the subset-limited version of edge betweenness centrality
 * by considering only those shortest paths that lie between vertices in a given
 * source and target subset.
 *
 * \param graph The graph object.
 * \param res The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param eids The edges for which the subset-limited betweenness centrality
 *        scores will be computed.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional weight vector for weighted
 *        betweenness. No edge weight may be NaN. Supply a null
 *        pointer here for unweighted betweenness.
 * \param sources A vertex selector for the sources of the shortest paths taken
 *        into considuration in the betweenness calculation.
 * \param targets A vertex selector for the targets of the shortest paths taken
 *        into considuration in the betweenness calculation.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for temporary data.
 *        \c IGRAPH_EINVVID, invalid vertex ID passed in \p sources or \p targets
 *
 * Time complexity: O(|S||E|), where
 * |S| is the number of vertices in the subset and
 * |E| is the number of edges in the graph.
 *
 * \sa \ref igraph_edge_betweenness() to compute the exact edge betweenness and
 * \ref igraph_edge_betweenness_cutoff() to compute the range-limited edge betweenness.
 */
igraph_error_t igraph_edge_betweenness_subset(const igraph_t *graph, igraph_vector_t *res,
                                   const igraph_es_t eids, igraph_bool_t directed,
                                   const igraph_vs_t sources, const igraph_vs_t targets,
                                   const igraph_vector_t *weights) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_sources;
    igraph_integer_t no_of_processed_sources;
    igraph_inclist_t inclist, parents;
    igraph_vit_t vit;
    igraph_eit_t eit;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_real_t *nrgeo;
    igraph_real_t *tmpscore;
    igraph_integer_t source, j;
    bool *is_target;
    igraph_stack_int_t S;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    IGRAPH_CHECK(igraph_vs_size(graph, &sources, &no_of_sources));

    is_target = IGRAPH_CALLOC(no_of_nodes, bool);
    IGRAPH_CHECK_OOM(is_target, "Insufficient memory for subset edge betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, is_target);

    IGRAPH_CHECK(igraph_vit_create(graph, targets, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        is_target[IGRAPH_VIT_GET(vit)] = true;
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_inclist_init_empty(&parents, no_of_nodes));
    IGRAPH_FINALLY(igraph_inclist_destroy, &parents);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    nrgeo = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(nrgeo, "Insufficient memory for subset edge betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, nrgeo);

    tmpscore = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(tmpscore, "Insufficient memory for subset edge betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, tmpscore);

    IGRAPH_CHECK(igraph_stack_int_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &S);

    if (!igraph_es_is_all(&eids)) {
        /* result needed only for a subset of the vertices */
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_edges);
    } else {
        /* result covers all vertices */
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_edges));
        igraph_vector_null(res);
        tmpres = res;
    }

    IGRAPH_CHECK(igraph_vit_create(graph, sources, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    for (
        no_of_processed_sources = 0, IGRAPH_VIT_RESET(vit);
        !IGRAPH_VIT_END(vit);
        IGRAPH_VIT_NEXT(vit), no_of_processed_sources++
    ) {
        source = IGRAPH_VIT_GET(vit);

        IGRAPH_PROGRESS(
            "Edge betweenness centrality (subset): ",
            100.0 * no_of_processed_sources / no_of_sources, 0
        );
        IGRAPH_ALLOW_INTERRUPTION();

        /* Loop invariant that is valid at this point:
         *
         * - the stack S is empty
         * - the 'dist' vector contains zeros only
         * - the 'nrgeo' array contains zeros only
         * - the 'tmpscore' array contains zeros only
         * - the 'parents' incidence list contains empty vectors only
         */

        /* TODO: there is more room for optimization here; the single-source
         * shortest path search runs until it reaches all the nodes in the
         * component of the source node even if we are only interested in a
         * smaller target subset. We could stop the search when all target
         * nodes were reached.
         */

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted_edge(graph, source, &dist, nrgeo, weights, &S, &parents, &inclist, -1));
        } else {
            IGRAPH_CHECK(igraph_i_sspf_edge(graph, source, &dist, nrgeo, &S, &parents, &inclist, -1));
        }

        /* Aggregate betweenness scores for the nodes we have reached in this
         * traversal */
        while (!igraph_stack_int_empty(&S)) {
            igraph_integer_t actnode = igraph_stack_int_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&parents, actnode);
            igraph_integer_t fatv_len = igraph_vector_int_size(fatv);
            igraph_real_t coeff;

            if (is_target[actnode]) {
                coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];
            } else {
                coeff = tmpscore[actnode] / nrgeo[actnode];
            }

            for (j = 0; j < fatv_len; j++) {
                igraph_integer_t father_edge = VECTOR(*fatv)[j];
                igraph_integer_t neighbor = IGRAPH_OTHER(graph, father_edge, actnode);
                tmpscore[neighbor] += nrgeo[neighbor] * coeff;
                VECTOR(*tmpres)[father_edge] += nrgeo[neighbor] * coeff;
            }

            /* Reset variables to ensure that the 'for' loop invariant will
             * still be valid in the next iteration */

            VECTOR(dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            tmpscore[actnode] = 0;
            igraph_vector_int_clear(fatv);
        }
    }

    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    /* Keep only the requested edges */
    if (!igraph_es_is_all(&eids)) {
        IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);

        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

        for (j = 0, IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit);
             IGRAPH_EIT_NEXT(eit), j++) {
            igraph_integer_t edge = IGRAPH_EIT_GET(eit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[edge];
        }

        igraph_eit_destroy(&eit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }


    if (!directed || !igraph_is_directed(graph)) {
        igraph_vector_scale(res, 0.5);
    }

    igraph_stack_int_destroy(&S);
    IGRAPH_FREE(tmpscore);
    IGRAPH_FREE(nrgeo);
    igraph_vector_destroy(&dist);
    igraph_inclist_destroy(&parents);
    igraph_inclist_destroy(&inclist);
    IGRAPH_FREE(is_target);
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}
