/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_scan.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_stack.h"
#include "igraph_structural.h"

#include "core/interruption.h"

/**
 * \section about_local_scan
 *
 * <para>
 * The scan statistic is a summary of the locality statistics that is computed
 * from the local neighborhood of each vertex. For details, see
 * Priebe, C. E., Conroy, J. M., Marchette, D. J., Park, Y. (2005).
 * Scan Statistics on Enron Graphs. Computational and Mathematical Organization Theory.
 * </para>
 */

/**
 * \function igraph_local_scan_0
 * Local scan-statistics, k=0
 *
 * K=0 scan-statistics is arbitrarily defined as the vertex degree for
 * unweighted, and the vertex strength for weighted graphs. See \ref
 * igraph_degree() and \ref igraph_strength().
 *
 * \param graph The input graph
 * \param res An initialized vector, the results are stored here.
 * \param weights Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param mode Type of the neighborhood, \c IGRAPH_OUT means outgoing,
 *        \c IGRAPH_IN means incoming and \c IGRAPH_ALL means all edges.
 * \return Error code.
 *
 */

igraph_error_t igraph_local_scan_0(const igraph_t *graph, igraph_vector_t *res,
                        const igraph_vector_t *weights,
                        igraph_neimode_t mode) {
    return igraph_strength(graph, res, igraph_vss_all(), mode, /*loops=*/ 1,
                    weights);
}

/* This removes loop, multiple edges and edges that point
   "backwards" according to the rank vector. It works on
   edge lists */

static igraph_error_t igraph_i_trans4_il_simplify(const igraph_t *graph, igraph_inclist_t *il,
                                       const igraph_vector_int_t *rank) {

    igraph_integer_t i;
    igraph_integer_t n = il->length;
    igraph_vector_int_t mark;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&mark, n);

    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &il->incs[i];
        igraph_integer_t j, l = igraph_vector_int_size(v);
        igraph_integer_t irank = VECTOR(*rank)[i];
        VECTOR(mark)[i] = i + 1;
        for (j = 0; j < l; /* nothing */) {
            igraph_integer_t edge = VECTOR(*v)[j];
            igraph_integer_t e = IGRAPH_OTHER(graph, edge, i);
            if (VECTOR(*rank)[e] > irank && VECTOR(mark)[e] != i + 1) {
                VECTOR(mark)[e] = i + 1;
                j++;
            } else {
                VECTOR(*v)[j] = igraph_vector_int_tail(v);
                igraph_vector_int_pop_back(v);
                l--;
            }
        }
    }

    igraph_vector_int_destroy(&mark);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;

}

/* This one handles both weighted and unweighted cases */

static igraph_error_t igraph_i_local_scan_1_directed(const igraph_t *graph,
                                          igraph_vector_t *res,
                                          const igraph_vector_t *weights,
                                          igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_inclist_t incs;
    igraph_integer_t i, node;

    igraph_vector_int_t neis;

    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {
        igraph_vector_int_t *edges1 = igraph_inclist_get(&incs, node);
        igraph_integer_t edgeslen1 = igraph_vector_int_size(edges1);

        IGRAPH_ALLOW_INTERRUPTION();

        /* Mark neighbors and self */
        VECTOR(neis)[node] = node + 1;
        for (i = 0; i < edgeslen1; i++) {
            igraph_integer_t e = VECTOR(*edges1)[i];
            igraph_integer_t nei = IGRAPH_OTHER(graph, e, node);
            igraph_real_t w = weights ? VECTOR(*weights)[e] : 1;
            VECTOR(neis)[nei] = node + 1;
            VECTOR(*res)[node] += w;
        }

        /* Crawl neighbors */
        for (i = 0; i < edgeslen1; i++) {
            igraph_integer_t e2 = VECTOR(*edges1)[i];
            igraph_integer_t nei = IGRAPH_OTHER(graph, e2, node);
            if (nei == node) {
                break;
            }
            igraph_vector_int_t *edges2 = igraph_inclist_get(&incs, nei);
            igraph_integer_t j, edgeslen2 = igraph_vector_int_size(edges2);
            for (j = 0; j < edgeslen2; j++) {
                igraph_integer_t e2 = VECTOR(*edges2)[j];
                igraph_integer_t nei2 = IGRAPH_OTHER(graph, e2, nei);
                igraph_real_t w2 = weights ? VECTOR(*weights)[e2] : 1;
                if (VECTOR(neis)[nei2] == node + 1) {
                    VECTOR(*res)[node] += w2;
                }
            }
        }

    } /* node < no_of_nodes */

    igraph_vector_int_destroy(&neis);
    igraph_inclist_destroy(&incs);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_local_scan_1_directed_all(const igraph_t *graph,
                                              igraph_vector_t *res,
                                              const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_inclist_t incs;
    igraph_integer_t i, node;

    igraph_vector_int_t neis;

    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, IGRAPH_ALL, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {
        igraph_vector_int_t *edges1 = igraph_inclist_get(&incs, node);
        igraph_integer_t edgeslen1 = igraph_vector_int_size(edges1);

        IGRAPH_ALLOW_INTERRUPTION();

        /* Mark neighbors. We also count the edges that are incident on ego.
           Note that this time we do not mark ego, because we don't want to
           double count its incident edges later, when we are going over the
           incident edges of ego's neighbors. */
        for (i = 0; i < edgeslen1; i++) {
            igraph_integer_t e = VECTOR(*edges1)[i];
            igraph_integer_t nei = IGRAPH_OTHER(graph, e, node);
            igraph_real_t w = weights ? VECTOR(*weights)[e] : 1;
            VECTOR(neis)[nei] = node + 1;
            VECTOR(*res)[node] += w;
        }

        /* Explicitly unmark ego in case it had a loop edge */
        VECTOR(neis)[node] = 0;

        /* Crawl neighbors. We make sure that each neighbor of 'node' is
           only crawled once. We count all qualifying edges of ego, and
           then unmark ego to avoid double counting. */
        for (i = 0; i < edgeslen1; i++) {
            igraph_integer_t e2 = VECTOR(*edges1)[i];
            igraph_integer_t nei = IGRAPH_OTHER(graph, e2, node);
            igraph_vector_int_t *edges2;
            igraph_integer_t j, edgeslen2;
            if (VECTOR(neis)[nei] != node + 1) {
                continue;
            }
            edges2 = igraph_inclist_get(&incs, nei);
            edgeslen2 = igraph_vector_int_size(edges2);
            for (j = 0; j < edgeslen2; j++) {
                igraph_integer_t e2 = VECTOR(*edges2)[j];
                igraph_integer_t nei2 = IGRAPH_OTHER(graph, e2, nei);
                igraph_real_t w2 = weights ? VECTOR(*weights)[e2] : 1;
                if (VECTOR(neis)[nei2] == node + 1) {
                    VECTOR(*res)[node] += w2;
                }
            }
            VECTOR(neis)[nei] = 0;
        }

    } /* node < no_of_nodes */

    igraph_vector_int_destroy(&neis);
    igraph_inclist_destroy(&incs);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_local_scan_1_ecount
 * Local scan-statistics, k=1, edge count and sum of weights
 *
 * Count the number of edges or the sum the edge weights in the
 * 1-neighborhood of vertices.
 *
 * \param graph The input graph
 * \param res An initialized vector, the results are stored here.
 * \param weights Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param mode Type of the neighborhood, \c IGRAPH_OUT means outgoing,
 *        \c IGRAPH_IN means incoming and \c IGRAPH_ALL means all edges.
 * \return Error code.
 *
 */

igraph_error_t igraph_local_scan_1_ecount(const igraph_t *graph, igraph_vector_t *res,
                               const igraph_vector_t *weights,
                               igraph_neimode_t mode) {

    if (igraph_is_directed(graph)) {
        if (mode != IGRAPH_ALL) {
            return igraph_i_local_scan_1_directed(graph, res, weights, mode);
        } else {
            return igraph_i_local_scan_1_directed_all(graph, res, weights);
        }
    } else {
        return igraph_local_scan_k_ecount(graph, 1, res, weights, mode);
    }
}

static igraph_error_t igraph_i_local_scan_0_them_w(const igraph_t *us, const igraph_t *them,
                                        igraph_vector_t *res,
                                        const igraph_vector_t *weights_them,
                                        igraph_neimode_t mode) {

    igraph_t is;
    igraph_vector_int_t map2;
    igraph_vector_t weights;
    igraph_integer_t i, m;

    if (!weights_them) {
        IGRAPH_ERROR("Edge weights not given for weighted scan-0",
                     IGRAPH_EINVAL);
    }
    if (igraph_vector_size(weights_them) != igraph_ecount(them)) {
        IGRAPH_ERROR("Invalid weights length for scan-0", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&map2, 0);
    IGRAPH_CHECK(igraph_intersection(&is, us, them, /* edge_map1= */ 0, &map2));
    IGRAPH_FINALLY(igraph_destroy, &is);

    /* Rewrite the map as edge weights */
    m = igraph_vector_int_size(&map2);
    IGRAPH_VECTOR_INIT_FINALLY(&weights, m);
    for (i = 0; i < m; i++) {
        VECTOR(weights)[i] = VECTOR(*weights_them)[ VECTOR(map2)[i] ];
    }

    IGRAPH_CHECK(igraph_strength(&is, res, igraph_vss_all(), mode, IGRAPH_LOOPS,
                    /*weights=*/ &weights));

    igraph_destroy(&is);
    igraph_vector_int_destroy(&map2);
    igraph_vector_destroy(&weights);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_local_scan_0_them
 * Local THEM scan-statistics, k=0
 *
 * K=0 scan-statistics is arbitrarily defined as the vertex degree for
 * unweighted, and the vertex strength for weighted graphs. See \ref
 * igraph_degree() and \ref igraph_strength().
 *
 * \param us The input graph, to use to extract the neighborhoods.
 * \param them The input graph to use for the actually counting.
 * \param res An initialized vector, the results are stored here.
 * \param weights_them Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param mode Type of the neighborhood, \c IGRAPH_OUT means outgoing,
 *        \c IGRAPH_IN means incoming and \c IGRAPH_ALL means all edges.
 * \return Error code.
 *
 */

igraph_error_t igraph_local_scan_0_them(const igraph_t *us, const igraph_t *them,
                             igraph_vector_t *res,
                             const igraph_vector_t *weights_them,
                             igraph_neimode_t mode) {

    igraph_t is;

    if (igraph_vcount(us) != igraph_vcount(them)) {
        IGRAPH_ERROR("Number of vertices don't match in scan-0", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(us) != igraph_is_directed(them)) {
        IGRAPH_ERROR("Directedness don't match in scan-0", IGRAPH_EINVAL);
    }

    if (weights_them) {
        return igraph_i_local_scan_0_them_w(us, them, res, weights_them, mode);
    }

    IGRAPH_CHECK(igraph_intersection(&is, us, them, /*edge_map1=*/ 0, /*edge_map2=*/ 0));
    IGRAPH_FINALLY(igraph_destroy, &is);

    IGRAPH_CHECK(igraph_strength(&is, res, igraph_vss_all(), mode, IGRAPH_LOOPS, /* weights = */ 0));

    igraph_destroy(&is);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_local_scan_1_ecount_them
 * Local THEM scan-statistics, k=1, edge count and sum of weights
 *
 * Count the number of edges or the sum the edge weights in the
 * 1-neighborhood of vertices.
 *
 * \param us The input graph to extract the neighborhoods.
 * \param them The input graph to perform the counting.
 * \param weights_them Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param mode Type of the neighborhood, \c IGRAPH_OUT means outgoing,
 *        \c IGRAPH_IN means incoming and \c IGRAPH_ALL means all edges.
 * \return Error code.
 *
 * \sa \ref igraph_local_scan_1_ecount() for the US statistics.
 */

igraph_error_t igraph_local_scan_1_ecount_them(const igraph_t *us, const igraph_t *them,
                                    igraph_vector_t *res,
                                    const igraph_vector_t *weights_them,
                                    igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(us);
    igraph_adjlist_t adj_us;
    igraph_inclist_t incs_them;
    igraph_vector_int_t neis;
    igraph_integer_t node;

    if (igraph_vcount(them) != no_of_nodes) {
        IGRAPH_ERROR("Number of vertices must match in scan-1", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(us) != igraph_is_directed(them)) {
        IGRAPH_ERROR("Directedness must match in scan-1", IGRAPH_EINVAL);
    }
    if (weights_them &&
        igraph_vector_size(weights_them) != igraph_ecount(them)) {
        IGRAPH_ERROR("Invalid weight vector length in scan-1 (them)",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_adjlist_init(
        us, &adj_us, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adj_us);
    IGRAPH_CHECK(igraph_inclist_init(them, &incs_them, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs_them);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {
        igraph_vector_int_t *neis_us = igraph_adjlist_get(&adj_us, node);
        igraph_vector_int_t *edges1_them = igraph_inclist_get(&incs_them, node);
        igraph_integer_t len1_us = igraph_vector_int_size(neis_us);
        igraph_integer_t len1_them = igraph_vector_int_size(edges1_them);
        igraph_integer_t i;

        IGRAPH_ALLOW_INTERRUPTION();

        /* Mark neighbors and self in us */
        VECTOR(neis)[node] = node + 1;
        for (i = 0; i < len1_us; i++) {
            igraph_integer_t nei = VECTOR(*neis_us)[i];
            VECTOR(neis)[nei] = node + 1;
        }

        /* Crawl neighbors in them, first ego */
        for (i = 0; i < len1_them; i++) {
            igraph_integer_t e = VECTOR(*edges1_them)[i];
            igraph_integer_t nei = IGRAPH_OTHER(them, e, node);
            if (VECTOR(neis)[nei] == node + 1) {
                igraph_real_t w = weights_them ? VECTOR(*weights_them)[e] : 1;
                VECTOR(*res)[node] += w;
            }
        }
        /* Then the rest */
        for (i = 0; i < len1_us; i++) {
            igraph_integer_t nei = VECTOR(*neis_us)[i];
            igraph_vector_int_t *edges2_them = igraph_inclist_get(&incs_them, nei);
            igraph_integer_t j, len2_them = igraph_vector_int_size(edges2_them);
            for (j = 0; j < len2_them; j++) {
                igraph_integer_t e2 = VECTOR(*edges2_them)[j];
                igraph_integer_t nei2 = IGRAPH_OTHER(them, e2, nei);
                if (VECTOR(neis)[nei2] == node + 1) {
                    igraph_real_t w = weights_them ? VECTOR(*weights_them)[e2] : 1;
                    VECTOR(*res)[node] += w;
                }
            }
        }

        /* For undirected, it was double counted */
        if (mode == IGRAPH_ALL || ! igraph_is_directed(us)) {
            VECTOR(*res)[node] /= 2.0;
        }

    } /* node < no_of_nodes */

    igraph_vector_int_destroy(&neis);
    igraph_inclist_destroy(&incs_them);
    igraph_adjlist_destroy(&adj_us);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_local_scan_k_ecount
 * \brief Sum the number of edges or the weights in k-neighborhood of every vertex.
 *
 * \param graph The input graph.
 * \param k The size of the neighborhood, non-negative integer.
 *        The k=0 case is special, see \ref igraph_local_scan_0().
 * \param res An initialized vector, the results are stored here.
 * \param weights Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param mode Type of the neighborhood, \c IGRAPH_OUT means outgoing,
 *        \c IGRAPH_IN means incoming and \c IGRAPH_ALL means all edges.
 * \return Error code.
 *
 */

igraph_error_t igraph_local_scan_k_ecount(const igraph_t *graph, igraph_integer_t k,
                               igraph_vector_t *res,
                               const igraph_vector_t *weights,
                               igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t node;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t marked;
    igraph_inclist_t incs;

    if (k < 0) {
        IGRAPH_ERROR("k must be non-negative in k-scan.", IGRAPH_EINVAL);
    }
    if (weights && igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERRORF("The weight vector length (%" IGRAPH_PRId ") in k-scan should equal "
                      "the number of edges of the graph (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(weights),
                      igraph_ecount(graph));
    }

    if (k == 0) {
        return igraph_local_scan_0(graph, res, weights, mode);
    }
    if (k == 1 && igraph_is_directed(graph)) {
        return igraph_local_scan_1_ecount(graph, res, weights, mode);
    }

    /* We do a BFS form each node, and simply count the number
       of edges on the way */

    IGRAPH_CHECK(igraph_dqueue_int_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &Q);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&marked, no_of_nodes);
    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0 ; node < no_of_nodes ; node++) {
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, node));
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, 0));
        VECTOR(marked)[node] = node + 1;
        while (!igraph_dqueue_int_empty(&Q)) {
            igraph_integer_t act = igraph_dqueue_int_pop(&Q);
            igraph_integer_t dist = igraph_dqueue_int_pop(&Q) + 1;
            igraph_vector_int_t *edges = igraph_inclist_get(&incs, act);
            igraph_integer_t i, edgeslen = igraph_vector_int_size(edges);
            for (i = 0; i < edgeslen; i++) {
                igraph_integer_t edge = VECTOR(*edges)[i];
                igraph_integer_t nei = IGRAPH_OTHER(graph, edge, act);
                if (dist <= k || VECTOR(marked)[nei] == node + 1) {
                    igraph_real_t w = weights ? VECTOR(*weights)[edge] : 1;
                    VECTOR(*res)[node] += w;
                }
                if (dist <= k && VECTOR(marked)[nei] != node + 1) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, nei));
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, dist));
                    VECTOR(marked)[nei] = node + 1;
                }
            }
        }

        if (mode == IGRAPH_ALL || ! igraph_is_directed(graph)) {
            VECTOR(*res)[node] /= 2.0;
        }

    } /* node < no_of_nodes */

    igraph_inclist_destroy(&incs);
    igraph_vector_int_destroy(&marked);
    igraph_dqueue_int_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_local_scan_k_ecount_them
 * \brief Local THEM scan-statistics, edge count or sum of weights.
 *
 * Count the number of edges or the sum the edge weights in the
 * k-neighborhood of vertices.
 *
 * \param us The input graph to extract the neighborhoods.
 * \param them The input graph to perform the counting.
 * \param k The size of the neighborhood, non-negative integer.
 *        The k=0 case is special, see \ref igraph_local_scan_0_them().
 * \param weights_them Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param mode Type of the neighborhood, \c IGRAPH_OUT means outgoing,
 *        \c IGRAPH_IN means incoming and \c IGRAPH_ALL means all edges.
 * \return Error code.
 *
 * \sa \ref igraph_local_scan_1_ecount() for the US statistics.
 */

igraph_error_t igraph_local_scan_k_ecount_them(const igraph_t *us, const igraph_t *them,
                                    igraph_integer_t k, igraph_vector_t *res,
                                    const igraph_vector_t *weights_them,
                                    igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(us);
    igraph_integer_t node;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t marked;
    igraph_stack_int_t ST;
    igraph_inclist_t incs_us, incs_them;

    if (igraph_vcount(them) != no_of_nodes) {
        IGRAPH_ERROR("The number of vertices in the two graphs must "
                "match in scan-k.",
                IGRAPH_EINVAL);
    }
    if (igraph_is_directed(us) != igraph_is_directed(them)) {
        IGRAPH_ERROR("Directedness in the two graphs must match "
                "in scan-k", IGRAPH_EINVAL);
    }
    if (k < 0) {
        IGRAPH_ERRORF("k must be non-negative in k-scan, got %" IGRAPH_PRId
                ".", IGRAPH_EINVAL, k);
    }
    if (weights_them &&
        igraph_vector_size(weights_them) != igraph_ecount(them)) {
        IGRAPH_ERRORF("The weight vector length (%" IGRAPH_PRId
            ") must be equal to the number of edges (%" IGRAPH_PRId
            ").", IGRAPH_EINVAL, igraph_vector_size(weights_them),
            igraph_ecount(them));
    }

    if (k == 0) {
        return igraph_local_scan_0_them(us, them, res, weights_them, mode);
    }
    if (k == 1) {
        return igraph_local_scan_1_ecount_them(us, them, res, weights_them, mode);
    }

    /* We mark the nodes in US in a BFS. Then we check the outgoing edges
       of all marked nodes in THEM. */

    IGRAPH_CHECK(igraph_dqueue_int_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &Q);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&marked, no_of_nodes);
    IGRAPH_CHECK(igraph_inclist_init(us, &incs_us, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs_us);
    IGRAPH_CHECK(igraph_inclist_init(them, &incs_them, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs_them);
    IGRAPH_CHECK(igraph_stack_int_init(&ST, 100));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &ST);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {

        /* BFS to mark the nodes in US */
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, node));
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, 0));
        IGRAPH_CHECK(igraph_stack_int_push(&ST, node));
        VECTOR(marked)[node] = node + 1;
        while (!igraph_dqueue_int_empty(&Q)) {
            igraph_integer_t act = igraph_dqueue_int_pop(&Q);
            igraph_integer_t dist = igraph_dqueue_int_pop(&Q) + 1;
            igraph_vector_int_t *edges = igraph_inclist_get(&incs_us, act);
            igraph_integer_t i, edgeslen = igraph_vector_int_size(edges);
            for (i = 0; i < edgeslen; i++) {
                igraph_integer_t edge = VECTOR(*edges)[i];
                igraph_integer_t nei = IGRAPH_OTHER(us, edge, act);
                if (dist <= k && VECTOR(marked)[nei] != node + 1) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, nei));
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, dist));
                    VECTOR(marked)[nei] = node + 1;
                    IGRAPH_CHECK(igraph_stack_int_push(&ST, nei));
                }
            }
        }

        /* Now check the edges of all nodes in THEM */
        while (!igraph_stack_int_empty(&ST)) {
            igraph_integer_t act = igraph_stack_int_pop(&ST);
            igraph_vector_int_t *edges = igraph_inclist_get(&incs_them, act);
            igraph_integer_t i, edgeslen = igraph_vector_int_size(edges);
            for (i = 0; i < edgeslen; i++) {
                igraph_integer_t edge = VECTOR(*edges)[i];
                igraph_integer_t nei = IGRAPH_OTHER(them, edge, act);
                if (VECTOR(marked)[nei] == node + 1) {
                    igraph_real_t w = weights_them ? VECTOR(*weights_them)[edge] : 1;
                    VECTOR(*res)[node] += w;
                }
            }
        }

        if (mode == IGRAPH_ALL || ! igraph_is_directed(us)) {
            VECTOR(*res)[node] /= 2;
        }

    } /* node < no_of_nodes */

    igraph_stack_int_destroy(&ST);
    igraph_inclist_destroy(&incs_them);
    igraph_inclist_destroy(&incs_us);
    igraph_vector_int_destroy(&marked);
    igraph_dqueue_int_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}
/**
 * \function igraph_local_scan_subset_ecount
 * \brief Local scan-statistics of subgraphs induced by subsets of vertices.
 *
 * Count the number of edges, or sum the edge weights in
 * induced subgraphs from vertices given as a parameter.
 *
 * \param graph The graph to perform the counting/summing in.
 * \param res Initialized vector, the result is stored here.
 * \param weights Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param subsets List of \type igraph_vector_int_t
 *        objects, the vertex subsets.
 * \return Error code.
 */

igraph_error_t igraph_local_scan_subset_ecount(const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vector_t *weights,
        const igraph_vector_int_list_t *subsets) {

    igraph_integer_t subset, no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_subsets = igraph_vector_int_list_size(subsets);
    igraph_inclist_t incs;
    igraph_vector_int_t marked;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (weights && igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length in local scan.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&marked, no_of_nodes);
    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, IGRAPH_OUT, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_subsets));
    igraph_vector_null(res);

    for (subset = 0; subset < no_of_subsets; subset++) {
        igraph_vector_int_t *nei = igraph_vector_int_list_get_ptr(subsets, subset);
        igraph_integer_t i, neilen = igraph_vector_int_size(nei);
        for (i = 0; i < neilen; i++) {
            igraph_integer_t vertex = VECTOR(*nei)[i];
            if (vertex < 0 || vertex >= no_of_nodes) {
                IGRAPH_ERROR("Invalid vertex ID in neighborhood list in local scan.",
                             IGRAPH_EINVAL);
            }
            VECTOR(marked)[vertex] = subset + 1;
        }

        for (i = 0; i < neilen; i++) {
            igraph_integer_t vertex = VECTOR(*nei)[i];
            igraph_vector_int_t *edges = igraph_inclist_get(&incs, vertex);
            igraph_integer_t j, edgeslen = igraph_vector_int_size(edges);
            for (j = 0; j < edgeslen; j++) {
                igraph_integer_t edge = VECTOR(*edges)[j];
                igraph_integer_t nei2 = IGRAPH_OTHER(graph, edge, vertex);
                if (VECTOR(marked)[nei2] == subset + 1) {
                    igraph_real_t w = weights ? VECTOR(*weights)[edge] : 1;
                    VECTOR(*res)[subset] += w;
                }
            }
        }
        if (!directed) {
            VECTOR(*res)[subset] /= 2.0;
        }
    }

    igraph_inclist_destroy(&incs);
    igraph_vector_int_destroy(&marked);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_local_scan_neighborhood_ecount
 * Local scan-statistics with pre-calculated neighborhoods
 *
 * Count the number of edges, or sum the edge weights in
 * neighborhoods given as a parameter.
 *
 * \deprecated-by igraph_local_scan_subset_ecount 0.10.0
 *
 * \param graph The graph to perform the counting/summing in.
 * \param res Initialized vector, the result is stored here.
 * \param weights Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param neighborhoods List of \type igraph_vector_int_t
 *        objects, the neighborhoods, one for each vertex in the
 *        graph.
 * \return Error code.
 */

igraph_error_t igraph_local_scan_neighborhood_ecount(const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vector_t *weights,
        const igraph_vector_int_list_t *neighborhoods) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (igraph_vector_int_list_size(neighborhoods) != no_of_nodes) {
        IGRAPH_ERROR("Invalid neighborhood list length in local scan.",
                     IGRAPH_EINVAL);
    }

    return igraph_local_scan_subset_ecount(graph, res, weights, neighborhoods);
}
