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
#include "igraph_progress.h"
#include "igraph_stack.h"

#include "core/indheap.h"
#include "core/interruption.h"
#include "core/math.h"

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

/*
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
 * \param  fathers adjacent list that starts empty and that stores the IDs
 *                 of the vertices that lead to a given node during the traversal
 * \param  adjlist the adjacency list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static int igraph_i_sspf(
    const igraph_t *graph,
    long int source,
    igraph_vector_t *dist, 
    double *nrgeo,
    igraph_stack_t *stack,
    igraph_adjlist_t *fathers,
    igraph_adjlist_t *adjlist,
    igraph_real_t cutoff
) {
    igraph_dqueue_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    long int nlen;

    IGRAPH_DQUEUE_INIT_FINALLY(&queue, 100);
    
    IGRAPH_CHECK(igraph_dqueue_push(&queue, source));
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_dqueue_empty(&queue)) {
        long int actnode = (long int) igraph_dqueue_pop(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && VECTOR(*dist)[actnode] > cutoff + 1) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            igraph_vector_int_clear(igraph_adjlist_get(fathers, actnode));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_push(stack, actnode));

        /* Examine the neighbors of this node */
        neis = igraph_adjlist_get(adjlist, actnode);
        nlen = igraph_vector_int_size(neis);
        for (int j = 0; j < nlen; j++) {
            long int neighbor = (long int) VECTOR(*neis)[j];

            if (VECTOR(*dist)[neighbor] == 0) {
                /* We have found 'neighbor' for the first time */
                VECTOR(*dist)[neighbor] = VECTOR(*dist)[actnode] + 1;
                IGRAPH_CHECK(igraph_dqueue_push(&queue, neighbor));
            }

            if (VECTOR(*dist)[neighbor] == VECTOR(*dist)[actnode] + 1 &&
                (VECTOR(*dist)[neighbor] <= cutoff + 1 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_adjlist_get(fathers, neighbor);
                IGRAPH_CHECK(igraph_vector_int_push_back(v, actnode));
                nrgeo[neighbor] += nrgeo[actnode];
            }
        }
    }

    igraph_dqueue_destroy(&queue);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS; 
}

/*
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
 * \param  fathers incidence list that starts empty and that stores the IDs
 *                 of the edges that lead to a given node during the traversal
 * \param  inclist the incidence list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static int igraph_i_sspf_edge( const igraph_t *graph, long int source, igraph_vector_t *dist, 
    double *nrgeo,
    igraph_stack_t *stack,
    igraph_inclist_t *fathers,
    const igraph_inclist_t *inclist,
    igraph_real_t cutoff
) {
    igraph_dqueue_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    long int nlen;

    IGRAPH_DQUEUE_INIT_FINALLY(&queue, 100);
    
    IGRAPH_CHECK(igraph_dqueue_push(&queue, source));
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_dqueue_empty(&queue)) {
        long int actnode = (long int) igraph_dqueue_pop(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && VECTOR(*dist)[actnode] > cutoff + 1) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            igraph_vector_int_clear(igraph_inclist_get(fathers, actnode));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_push(stack, actnode));

        /* Examine the neighbors of this node */
        neis = igraph_inclist_get(inclist, actnode);
        nlen = igraph_vector_int_size(neis);
        for (int j = 0; j < nlen; j++) {
            long int edge = (long int) VECTOR(*neis)[j];
            long int neighbor = IGRAPH_OTHER(graph, edge, actnode);

            if (VECTOR(*dist)[neighbor] == 0) {
                /* We have found 'neighbor' for the first time */
                VECTOR(*dist)[neighbor] = VECTOR(*dist)[actnode] + 1;
                IGRAPH_CHECK(igraph_dqueue_push(&queue, neighbor));
            }

            if (VECTOR(*dist)[neighbor] == VECTOR(*dist)[actnode] + 1 &&
                (VECTOR(*dist)[neighbor] <= cutoff + 1 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_inclist_get(fathers, neighbor);
                IGRAPH_CHECK(igraph_vector_int_push_back(v, edge));
                nrgeo[neighbor] += nrgeo[actnode];
            }
        }
    }

    igraph_dqueue_destroy(&queue);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS; 
}

/*
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
 * \param  fathers adjacency list that starts empty and that stores the IDs
 *                 of the vertices that lead to a given node during the traversal
 * \param  inclist the incidence list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static int igraph_i_sspf_weighted(
    const igraph_t *graph,
    long int source, 
    igraph_vector_t *dist, 
    double *nrgeo, 
    const igraph_vector_t *weights,
    igraph_stack_t *stack,
    igraph_adjlist_t *fathers,
    igraph_inclist_t *inclist,
    igraph_real_t cutoff
) {
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;
    
    int cmp_result;
    igraph_2wheap_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    long int nlen;
    
    /* TODO: this is an O|V| step here. We could save some time by pre-allocating
     * the two-way heap in the caller and re-using it here */
    IGRAPH_CHECK(igraph_2wheap_init(&queue, (long int) igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &queue);

    igraph_2wheap_push_with_index(&queue, source, -1.0);
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_2wheap_empty(&queue)) {
        long int minnei = igraph_2wheap_max_index(&queue);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && mindist > cutoff + 1.0) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[minnei] = 0;
            nrgeo[minnei] = 0;
            igraph_vector_int_clear(igraph_adjlist_get(fathers, minnei));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_push(stack, minnei));

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
                v = igraph_adjlist_get(fathers, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = minnei;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&queue, to, -altdist));
            } else if (cmp_result < 0) {
                /* This is a shorter path */
                v = igraph_adjlist_get(fathers, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = minnei;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_modify(&queue, to, -altdist));
            } else if (cmp_result == 0 && (altdist <= cutoff + 1.0 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_adjlist_get(fathers, to);
                IGRAPH_CHECK(igraph_vector_int_push_back(v, minnei));
                nrgeo[to] += nrgeo[minnei];
            }
        }
    }

    igraph_2wheap_destroy(&queue);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/*
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
 * \param  fathers incidence list that starts empty and that stores the IDs
 *                 of the edges that lead to a given node during the traversal
 * \param  inclist the incidence list of the graph
 * \param  cutoff  cutoff length of shortest paths
 */
static int igraph_i_sspf_weighted_edge(
    const igraph_t *graph,
    long int source,
    igraph_vector_t *dist, 
    double *nrgeo, 
    const igraph_vector_t *weights,
    igraph_stack_t *stack,
    igraph_inclist_t *fathers,
    const igraph_inclist_t *inclist,
    igraph_real_t cutoff
) {
    const double eps = IGRAPH_SHORTEST_PATH_EPSILON;

    int cmp_result;
    igraph_2wheap_t queue;
    const igraph_vector_int_t *neis;
    igraph_vector_int_t *v;
    long int nlen;
    
    /* TODO: this is an O|V| step here. We could save some time by pre-allocating
     * the two-way heap in the caller and re-using it here */
    IGRAPH_CHECK(igraph_2wheap_init(&queue, (long int) igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &queue);

    igraph_2wheap_push_with_index(&queue, source, -1.0);
    VECTOR(*dist)[source] = 1.0;
    nrgeo[source] = 1;

    while (!igraph_2wheap_empty(&queue)) {
        long int minnei = igraph_2wheap_max_index(&queue);
        igraph_real_t mindist = -igraph_2wheap_delete_max(&queue);

        /* Ignore vertices that are more distant than the cutoff */
        if (cutoff >= 0 && mindist > cutoff + 1.0) {
            /* Reset variables if node is too distant */
            VECTOR(*dist)[minnei] = 0;
            nrgeo[minnei] = 0;
            igraph_vector_int_clear(igraph_inclist_get(fathers, minnei));
            continue;
        }

        /* Record that we have visited this node */
        IGRAPH_CHECK(igraph_stack_push(stack, minnei));

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
                v = igraph_inclist_get(fathers, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = edge;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&queue, to, -altdist));
            } else if (cmp_result < 0) {
                /* This is a shorter path */
                v = igraph_inclist_get(fathers, to);
                IGRAPH_CHECK(igraph_vector_int_resize(v, 1));
                VECTOR(*v)[0] = edge;
                nrgeo[to] = nrgeo[minnei];
                VECTOR(*dist)[to] = altdist;
                IGRAPH_CHECK(igraph_2wheap_modify(&queue, to, -altdist));
            } else if (cmp_result == 0 && (altdist <= cutoff + 1.0 || cutoff < 0)) {
                /* Only add if the node is not more distant than the cutoff */
                v = igraph_inclist_get(fathers, to);
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
            } else if (igraph_is_nan(minweight)) {
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

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_adjlist_t adjlist, fathers;
    igraph_inclist_t inclist;
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

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    }

    IGRAPH_CHECK(igraph_adjlist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &fathers);

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
         * - the 'fathers' adjacency list contains empty vectors only
         */

        IGRAPH_PROGRESS("Betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted(graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, cutoff));
        } else {
            IGRAPH_CHECK(igraph_i_sspf(graph, source, &dist, nrgeo, &S, &fathers, &adjlist, cutoff));
        }

        /* Aggregate betweenness scores for the nodes we have reached in this
         * traversal */
        while (!igraph_stack_empty(&S)) {
            long int actnode = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *neis = igraph_adjlist_get(&fathers, actnode);
            long int nneis = igraph_vector_int_size(neis);
            double coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];

            for (j = 0; j < nneis; j++) {
                neighbor = (long int) VECTOR(*neis)[j];
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
            long int node = IGRAPH_VIT_GET(vit);
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

    igraph_Free(nrgeo);
    igraph_Free(tmpscore);
    igraph_vector_destroy(&dist);
    igraph_stack_destroy(&S);
    igraph_adjlist_destroy(&fathers);
    if (weights) {
        igraph_inclist_destroy(&inclist);
    } else {
        igraph_adjlist_destroy(&adjlist);
    }
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
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
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_inclist_t inclist, fathers;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    double *nrgeo;
    double *tmpscore;
    long int source, j;
    igraph_stack_t S;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

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

        /* Loop invariant that is valid at this point:
         *
         * - the stack S is empty
         * - the 'dist' vector contains zeros only
         * - the 'nrgeo' array contains zeros only
         * - the 'tmpscore' array contains zeros only
         * - the 'fathers' incidence list contains empty vectors only
         */

        IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted_edge(graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, cutoff));
        } else {
            IGRAPH_CHECK(igraph_i_sspf_edge(graph, source, &dist, nrgeo, &S, &fathers, &inclist, cutoff));
        }

        /* Aggregate betweenness scores for the edges we have reached in this
         * traversal */
        while (!igraph_stack_empty(&S)) {
            long int actnode = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&fathers, actnode);
            long int fatv_len = igraph_vector_int_size(fatv);
            double coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];

            for (j = 0; j < fatv_len; j++) {
                long int fedge = (long int) VECTOR(*fatv)[j];
                long int neighbor = IGRAPH_OTHER(graph, fedge, actnode);
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

    igraph_stack_destroy(&S);
    igraph_inclist_destroy(&inclist);
    igraph_inclist_destroy(&fathers);
    igraph_vector_destroy(&dist);
    igraph_Free(tmpscore);
    igraph_Free(nrgeo);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness_estimate
 * \brief Estimated betweenness centrality of the edges.
 *
 * \deprecated-by igraph_edge_betweenness_cutoff 0.9
 *
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

/**
 * \ingroup structural
 * \function igraph_betweenness_subset
 * \brief Betweenness centrality for a subset of source and target vertices.
 *
 * </para><para>
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
 * Time complexity: O(|S||E|),
 * |S| The number of vertices in the subset
 * |E| The number of edges in the graph
 *
 * \sa \ref igraph_betweenness() to calculate the exact vertex betweenness and
 * \ref igraph_betweenness_cutoff() to calculate the range-limited vertex
 * betweenness.
 */
int igraph_betweenness_subset(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_vs_t vids, igraph_bool_t directed,
                              const igraph_vs_t sources, const igraph_vs_t targets,
                              const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_sources;
    igraph_integer_t no_of_processed_sources;
    igraph_adjlist_t adjlist, fathers;
    igraph_inclist_t inclist;
    long int source, j;
    igraph_stack_t S;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    long int father;
    igraph_vector_t dist;
    double *nrgeo;
    double *tmpscore;
    igraph_vit_t vit;
    unsigned char *is_target;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    IGRAPH_CHECK(igraph_vs_size(graph, &sources, &no_of_sources));

    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    }

    IGRAPH_CHECK(igraph_adjlist_init_empty(&fathers, no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &fathers);

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

    IGRAPH_CHECK(igraph_vit_create(graph, targets, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        is_target[(long int) IGRAPH_VIT_GET(vit)] = 1;
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
         * - the 'fathers' adjacency list contains empty vectors only
         */

        /* TODO: there is more room for optimization here; the single-source
         * shortest path search runs until it reaches all the nodes in the
         * component of the source node even if we are only interested in a
         * smaller target subset. We could stop the search when all target
         * nodes were reached.
         */

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted(graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, -1));
        } else {
            IGRAPH_CHECK(igraph_i_sspf(graph, source, &dist, nrgeo, &S, &fathers, &adjlist, -1));
        }

        /* Aggregate betweenness scores for the nodes we have reached in this
         * traversal */
        while (!igraph_stack_empty(&S)) {
            long int actnode = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_adjlist_get(&fathers, actnode);
            long int fatv_len = igraph_vector_int_size(fatv);
            double coeff;

            if (is_target[actnode]) {
                coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];
            } else {
                coeff = tmpscore[actnode] / nrgeo[actnode];
            }

            for (j = 0; j < fatv_len; j++) {
                father = (long int) VECTOR(*fatv)[j];
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
            long int node = IGRAPH_VIT_GET(vit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[node];
        }

        igraph_vit_destroy(&vit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }

   if (!directed || !igraph_is_directed(graph)) {
        igraph_vector_scale(res, 0.5);
    }

    igraph_Free(is_target);
    igraph_Free(tmpscore);
    igraph_Free(nrgeo);
    igraph_vector_destroy(&dist);
    igraph_stack_destroy(&S);
    igraph_adjlist_destroy(&fathers);
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
 * </para><para>
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
 * Time complexity: O(|S||E|),
 * |S| The number of vertices in the subset
 * |E| The number of edges in the graph
 *
 * \sa \ref igraph_edge_betweenness() to compute the exact edge betweenness and
 * \ref igraph_edge_betweenness_cutoff() to compute the range-limited edge betweenness.
 */
int igraph_edge_betweenness_subset(const igraph_t *graph, igraph_vector_t *res,
                                   const igraph_es_t eids, igraph_bool_t directed,
                                   const igraph_vs_t sources, const igraph_vs_t targets,
                                   const igraph_vector_t *weights) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_sources;
    igraph_integer_t no_of_processed_sources;
    igraph_inclist_t inclist, fathers;
    igraph_vit_t vit;
    igraph_eit_t eit;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    double *nrgeo;
    double *tmpscore;
    long int source, j;
    unsigned char *is_target;
    igraph_stack_t S;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    IGRAPH_CHECK(igraph_vs_size(graph, &sources, &no_of_sources));

    is_target = igraph_Calloc(no_of_nodes, unsigned char);
    if (is_target == 0) {
        IGRAPH_ERROR("Insufficient memory for edge betweenness calculation.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, is_target);

    IGRAPH_CHECK(igraph_vit_create(graph, targets, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        is_target[(long int) IGRAPH_VIT_GET(vit)] = 1;
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

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
         * - the 'fathers' incidence list contains empty vectors only
         */

        /* TODO: there is more room for optimization here; the single-source
         * shortest path search runs until it reaches all the nodes in the
         * component of the source node even if we are only interested in a
         * smaller target subset. We could stop the search when all target
         * nodes were reached.
         */

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted_edge(graph, source, &dist, nrgeo, weights, &S, &fathers, &inclist, -1));
        } else {
            IGRAPH_CHECK(igraph_i_sspf_edge(graph, source, &dist, nrgeo, &S, &fathers, &inclist, -1));
        }

        /* Aggregate betweenness scores for the nodes we have reached in this
         * traversal */
        while (!igraph_stack_empty(&S)) {
            long int actnode = (long int) igraph_stack_pop(&S);
            igraph_vector_int_t *fatv = igraph_inclist_get(&fathers, actnode);
            long int fatv_len = igraph_vector_int_size(fatv);
            double coeff;

            if (is_target[actnode]) {
                coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];
            } else {
                coeff = tmpscore[actnode] / nrgeo[actnode];
                // coeff = tmpscore[actnode] / fatv_len;
            }

            for (j = 0; j < fatv_len; j++) {
                long int father_edge = (long int) VECTOR(*fatv)[j];
                long int neighbor = IGRAPH_OTHER(graph, father_edge, actnode);
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
            long int edge = IGRAPH_EIT_GET(eit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[edge];
        }

        igraph_eit_destroy(&eit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }


    if (!directed || !igraph_is_directed(graph)) {
        igraph_vector_scale(res, 0.5);
    }

    igraph_stack_destroy(&S);
    igraph_Free(tmpscore);
    igraph_Free(nrgeo);
    igraph_vector_destroy(&dist);
    igraph_inclist_destroy(&fathers);
    igraph_inclist_destroy(&inclist);
    igraph_Free(is_target);
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}
