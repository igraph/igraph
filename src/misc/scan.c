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
#include "igraph_arpack.h"
#include "igraph_centrality.h"
#include "igraph_dqueue.h"
#include "igraph_eigen.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_stack.h"
#include "igraph_structural.h"

#include "core/interruption.h"
#include "properties/properties_internal.h"

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

int igraph_local_scan_0(const igraph_t *graph, igraph_vector_t *res,
                        const igraph_vector_t *weights,
                        igraph_neimode_t mode) {
    if (weights) {
        igraph_strength(graph, res, igraph_vss_all(), mode, /*loops=*/ 1,
                        weights);
    } else {
        igraph_degree(graph, res, igraph_vss_all(), mode, /*loops=*/ 1);
    }
    return 0;
}

/* This removes loop, multiple edges and edges that point
   "backwards" according to the rank vector. It works on
   edge lists */

static int igraph_i_trans4_il_simplify(const igraph_t *graph, igraph_inclist_t *il,
                                       const igraph_vector_int_t *rank) {

    long int i;
    long int n = il->length;
    igraph_vector_int_t mark;
    igraph_vector_int_init(&mark, n);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &mark);

    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &il->incs[i];
        int j, l = igraph_vector_int_size(v);
        int irank = VECTOR(*rank)[i];
        VECTOR(mark)[i] = i + 1;
        for (j = 0; j < l; /* nothing */) {
            long int edge = (long int) VECTOR(*v)[j];
            long int e = IGRAPH_OTHER(graph, edge, i);
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
    return 0;

}

/* This one handles both weighted and unweighted cases */

static int igraph_i_local_scan_1_directed(const igraph_t *graph,
                                          igraph_vector_t *res,
                                          const igraph_vector_t *weights,
                                          igraph_neimode_t mode) {

    int no_of_nodes = igraph_vcount(graph);
    igraph_inclist_t incs;
    int i, node;

    igraph_vector_int_t neis;

    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    igraph_vector_int_init(&neis, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

    igraph_vector_resize(res, no_of_nodes);
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {
        igraph_vector_int_t *edges1 = igraph_inclist_get(&incs, node);
        int edgeslen1 = igraph_vector_int_size(edges1);

        IGRAPH_ALLOW_INTERRUPTION();

        /* Mark neighbors and self*/
        VECTOR(neis)[node] = node + 1;
        for (i = 0; i < edgeslen1; i++) {
            int e = VECTOR(*edges1)[i];
            int nei = IGRAPH_OTHER(graph, e, node);
            igraph_real_t w = weights ? VECTOR(*weights)[e] : 1;
            VECTOR(neis)[nei] = node + 1;
            VECTOR(*res)[node] += w;
        }

        /* Crawl neighbors */
        for (i = 0; i < edgeslen1; i++) {
            int e2 = VECTOR(*edges1)[i];
            int nei = IGRAPH_OTHER(graph, e2, node);
            igraph_vector_int_t *edges2 = igraph_inclist_get(&incs, nei);
            int j, edgeslen2 = igraph_vector_int_size(edges2);
            for (j = 0; j < edgeslen2; j++) {
                int e2 = VECTOR(*edges2)[j];
                int nei2 = IGRAPH_OTHER(graph, e2, nei);
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

    return 0;
}

static int igraph_i_local_scan_1_directed_all(const igraph_t *graph,
                                              igraph_vector_t *res,
                                              const igraph_vector_t *weights) {

    int no_of_nodes = igraph_vcount(graph);
    igraph_inclist_t incs;
    int i, node;

    igraph_vector_int_t neis;

    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    igraph_vector_int_init(&neis, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

    igraph_vector_resize(res, no_of_nodes);
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {
        igraph_vector_int_t *edges1 = igraph_inclist_get(&incs, node);
        int edgeslen1 = igraph_vector_int_size(edges1);

        IGRAPH_ALLOW_INTERRUPTION();

        /* Mark neighbors. We also count the edges that are incident to ego.
           Note that this time we do not mark ego, because we don't want to
           double count its incident edges later, when we are going over the
           incident edges of ego's neighbors. */
        for (i = 0; i < edgeslen1; i++) {
            int e = VECTOR(*edges1)[i];
            int nei = IGRAPH_OTHER(graph, e, node);
            igraph_real_t w = weights ? VECTOR(*weights)[e] : 1;
            VECTOR(neis)[nei] = node + 1;
            VECTOR(*res)[node] += w;
        }

        /* Crawl neighbors. We make sure that each neighbor of 'node' is
           only crawed once. We count all qualifying edges of ego, and
           then unmark ego to avoid double counting. */
        for (i = 0; i < edgeslen1; i++) {
            int e2 = VECTOR(*edges1)[i];
            int nei = IGRAPH_OTHER(graph, e2, node);
            igraph_vector_int_t *edges2;
            int j, edgeslen2;
            if (VECTOR(neis)[nei] != node + 1) {
                continue;
            }
            edges2 = igraph_inclist_get(&incs, nei);
            edgeslen2 = igraph_vector_int_size(edges2);
            for (j = 0; j < edgeslen2; j++) {
                int e2 = VECTOR(*edges2)[j];
                int nei2 = IGRAPH_OTHER(graph, e2, nei);
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

    return 0;
}

static int igraph_i_local_scan_1_sumweights(const igraph_t *graph,
                                            igraph_vector_t *res,
                                            const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    long int node, i, j, nn;
    igraph_inclist_t allinc;
    igraph_vector_int_t *neis1, *neis2;
    long int neilen1, neilen2;
    long int *neis;
    long int maxdegree;

    igraph_vector_int_t order;
    igraph_vector_int_t rank;
    igraph_vector_t degree, *edge1 = &degree; /* reuse degree as edge1 */

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    igraph_vector_int_init(&order, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &order);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL,
                               IGRAPH_LOOPS));
    maxdegree = (long int) igraph_vector_max(&degree) + 1;
    igraph_vector_order1_int(&degree, &order, maxdegree);
    igraph_vector_int_init(&rank, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &rank);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(rank)[ VECTOR(order)[i] ] = no_of_nodes - i - 1;
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &allinc, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &allinc);
    IGRAPH_CHECK(igraph_i_trans4_il_simplify(graph, &allinc, &rank));

    neis = IGRAPH_CALLOC(no_of_nodes, long int);
    if (neis == 0) {
        IGRAPH_ERROR("undirected local transitivity failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, neis);

    IGRAPH_CHECK(igraph_strength(graph, res, igraph_vss_all(), IGRAPH_ALL,
                                 IGRAPH_LOOPS, weights));

    for (nn = no_of_nodes - 1; nn >= 0; nn--) {
        node = VECTOR(order)[nn];

        IGRAPH_ALLOW_INTERRUPTION();

        neis1 = igraph_inclist_get(&allinc, node);
        neilen1 = igraph_vector_int_size(neis1);

        /* Mark the neighbors of the node */
        for (i = 0; i < neilen1; i++) {
            int edge = VECTOR(*neis1)[i];
            int nei = IGRAPH_OTHER(graph, edge, node);
            VECTOR(*edge1)[nei] = VECTOR(*weights)[edge];
            neis[nei] = node + 1;
        }

        for (i = 0; i < neilen1; i++) {
            long int edge = VECTOR(*neis1)[i];
            long int nei = IGRAPH_OTHER(graph, edge, node);
            igraph_real_t w = VECTOR(*weights)[edge];
            neis2 = igraph_inclist_get(&allinc, nei);
            neilen2 = igraph_vector_int_size(neis2);
            for (j = 0; j < neilen2; j++) {
                long int edge2 = VECTOR(*neis2)[j];
                long int nei2 = IGRAPH_OTHER(graph, edge2, nei);
                igraph_real_t w2 = VECTOR(*weights)[edge2];
                if (neis[nei2] == node + 1) {
                    VECTOR(*res)[node] += w2;
                    VECTOR(*res)[nei2] += w;
                    VECTOR(*res)[nei] += VECTOR(*edge1)[nei2];
                }
            }
        }
    }

    igraph_free(neis);
    igraph_inclist_destroy(&allinc);
    igraph_vector_int_destroy(&rank);
    igraph_vector_destroy(&degree);
    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
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

int igraph_local_scan_1_ecount(const igraph_t *graph, igraph_vector_t *res,
                               const igraph_vector_t *weights,
                               igraph_neimode_t mode) {

    if (igraph_is_directed(graph)) {
        if (mode != IGRAPH_ALL) {
            return igraph_i_local_scan_1_directed(graph, res, weights, mode);
        } else {
            return igraph_i_local_scan_1_directed_all(graph, res, weights);
        }
    } else {
        if (weights) {
            return igraph_i_local_scan_1_sumweights(graph, res, weights);
        } else {

#define TRIEDGES
#include "properties/triangles_template.h"
#undef TRIEDGES

        }
    }

    return 0;
}

static int igraph_i_local_scan_0_them_w(const igraph_t *us, const igraph_t *them,
                                        igraph_vector_t *res,
                                        const igraph_vector_t *weights_them,
                                        igraph_neimode_t mode) {

    igraph_t is;
    igraph_vector_t map2;
    int i, m;

    if (!weights_them) {
        IGRAPH_ERROR("Edge weights not given for weighted scan-0",
                     IGRAPH_EINVAL);
    }
    if (igraph_vector_size(weights_them) != igraph_ecount(them)) {
        IGRAPH_ERROR("Invalid weights length for scan-0", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&map2, 0);
    igraph_intersection(&is, us, them, /*map1=*/ 0, &map2);
    IGRAPH_FINALLY(igraph_destroy, &is);

    /* Rewrite the map as edge weights */
    m = igraph_vector_size(&map2);
    for (i = 0; i < m; i++) {
        VECTOR(map2)[i] = VECTOR(*weights_them)[ (int) VECTOR(map2)[i] ];
    }

    igraph_strength(&is, res, igraph_vss_all(), mode, IGRAPH_LOOPS,
                    /*weights=*/ &map2);

    igraph_destroy(&is);
    igraph_vector_destroy(&map2);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
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

int igraph_local_scan_0_them(const igraph_t *us, const igraph_t *them,
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

    igraph_intersection(&is, us, them, /*edgemap1=*/ 0, /*edgemap2=*/ 0);
    IGRAPH_FINALLY(igraph_destroy, &is);

    igraph_degree(&is, res, igraph_vss_all(), mode, IGRAPH_LOOPS);

    igraph_destroy(&is);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
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

int igraph_local_scan_1_ecount_them(const igraph_t *us, const igraph_t *them,
                                    igraph_vector_t *res,
                                    const igraph_vector_t *weights_them,
                                    igraph_neimode_t mode) {

    int no_of_nodes = igraph_vcount(us);
    igraph_adjlist_t adj_us;
    igraph_inclist_t incs_them;
    igraph_vector_int_t neis;
    int node;

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

    IGRAPH_CHECK(igraph_vector_int_init(&neis, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {
        igraph_vector_int_t *neis_us = igraph_adjlist_get(&adj_us, node);
        igraph_vector_int_t *edges1_them = igraph_inclist_get(&incs_them, node);
        int len1_us = igraph_vector_int_size(neis_us);
        int len1_them = igraph_vector_int_size(edges1_them);
        int i;

        IGRAPH_ALLOW_INTERRUPTION();

        /* Mark neighbors and self in us */
        VECTOR(neis)[node] = node + 1;
        for (i = 0; i < len1_us; i++) {
            int nei = VECTOR(*neis_us)[i];
            VECTOR(neis)[nei] = node + 1;
        }

        /* Crawl neighbors in them, first ego */
        for (i = 0; i < len1_them; i++) {
            int e = VECTOR(*edges1_them)[i];
            int nei = IGRAPH_OTHER(them, e, node);
            if (VECTOR(neis)[nei] == node + 1) {
                igraph_real_t w = weights_them ? VECTOR(*weights_them)[e] : 1;
                VECTOR(*res)[node] += w;
            }
        }
        /* Then the rest */
        for (i = 0; i < len1_us; i++) {
            int nei = VECTOR(*neis_us)[i];
            igraph_vector_int_t *edges2_them = igraph_inclist_get(&incs_them, nei);
            int j, len2_them = igraph_vector_int_size(edges2_them);
            for (j = 0; j < len2_them; j++) {
                int e2 = VECTOR(*edges2_them)[j];
                int nei2 = IGRAPH_OTHER(them, e2, nei);
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

    return 0;
}

/**
 * \function igraph_local_scan_k_ecount
 * Local scan-statistics, general function, edge count and sum of weights
 *
 * Count the number of edges or the sum the edge weights in the
 * k-neighborhood of vertices.
 *
 * \param graph The input graph
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

int igraph_local_scan_k_ecount(const igraph_t *graph, int k,
                               igraph_vector_t *res,
                               const igraph_vector_t *weights,
                               igraph_neimode_t mode) {

    int no_of_nodes = igraph_vcount(graph);
    int node;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t marked;
    igraph_inclist_t incs;

    if (k < 0) {
        IGRAPH_ERROR("k must be non-negative in k-scan", IGRAPH_EINVAL);
    }
    if (weights && igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length in k-scan", IGRAPH_EINVAL);
    }

    if (k == 0) {
        return igraph_local_scan_0(graph, res, weights, mode);
    }
    if (k == 1) {
        return igraph_local_scan_1_ecount(graph, res, weights, mode);
    }

    /* We do a BFS form each node, and simply count the number
       of edges on the way */

    IGRAPH_CHECK(igraph_dqueue_int_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &Q);
    IGRAPH_CHECK(igraph_vector_int_init(&marked, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &marked);
    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0 ; node < no_of_nodes ; node++) {
        igraph_dqueue_int_push(&Q, node);
        igraph_dqueue_int_push(&Q, 0);
        VECTOR(marked)[node] = node + 1;
        while (!igraph_dqueue_int_empty(&Q)) {
            int act = igraph_dqueue_int_pop(&Q);
            int dist = igraph_dqueue_int_pop(&Q) + 1;
            igraph_vector_int_t *edges = igraph_inclist_get(&incs, act);
            int i, edgeslen = igraph_vector_int_size(edges);
            for (i = 0; i < edgeslen; i++) {
                int edge = VECTOR(*edges)[i];
                int nei = IGRAPH_OTHER(graph, edge, act);
                if (dist <= k || VECTOR(marked)[nei] == node + 1) {
                    igraph_real_t w = weights ? VECTOR(*weights)[edge] : 1;
                    VECTOR(*res)[node] += w;
                }
                if (dist <= k && VECTOR(marked)[nei] != node + 1) {
                    igraph_dqueue_int_push(&Q, nei);
                    igraph_dqueue_int_push(&Q, dist);
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

    return 0;
}

/**
 * \function igraph_local_scan_k_ecount_them
 * Local THEM scan-statistics, general function, edge count and sum of weights
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

int igraph_local_scan_k_ecount_them(const igraph_t *us, const igraph_t *them,
                                    int k, igraph_vector_t *res,
                                    const igraph_vector_t *weights_them,
                                    igraph_neimode_t mode) {

    int no_of_nodes = igraph_vcount(us);
    int node;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t marked;
    igraph_stack_int_t ST;
    igraph_inclist_t incs_us, incs_them;

    if (igraph_vcount(them) != no_of_nodes) {
        IGRAPH_ERROR("Number of vertices must match in scan-k", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(us) != igraph_is_directed(them)) {
        IGRAPH_ERROR("Directedness must match in scan-k", IGRAPH_EINVAL);
    }
    if (k < 0) {
        IGRAPH_ERROR("k must be non-negative in k-scan", IGRAPH_EINVAL);
    }
    if (weights_them &&
        igraph_vector_size(weights_them) != igraph_ecount(them)) {
        IGRAPH_ERROR("Invalid weight vector length in k-scan (them)",
                     IGRAPH_EINVAL);
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
    IGRAPH_CHECK(igraph_vector_int_init(&marked, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &marked);
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
            int act = igraph_dqueue_int_pop(&Q);
            int dist = igraph_dqueue_int_pop(&Q) + 1;
            igraph_vector_int_t *edges = igraph_inclist_get(&incs_us, act);
            int i, edgeslen = igraph_vector_int_size(edges);
            for (i = 0; i < edgeslen; i++) {
                int edge = VECTOR(*edges)[i];
                int nei = IGRAPH_OTHER(us, edge, act);
                if (dist <= k && VECTOR(marked)[nei] != node + 1) {
                    igraph_dqueue_int_push(&Q, nei);
                    igraph_dqueue_int_push(&Q, dist);
                    VECTOR(marked)[nei] = node + 1;
                    igraph_stack_int_push(&ST, nei);
                }
            }
        }

        /* Now check the edges of all nodes in THEM */
        while (!igraph_stack_int_empty(&ST)) {
            int act = igraph_stack_int_pop(&ST);
            igraph_vector_int_t *edges = igraph_inclist_get(&incs_them, act);
            int i, edgeslen = igraph_vector_int_size(edges);
            for (i = 0; i < edgeslen; i++) {
                int edge = VECTOR(*edges)[i];
                int nei = IGRAPH_OTHER(them, edge, act);
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

    return 0;
}

/**
 * \function igraph_local_scan_neighborhood_ecount
 * Local scan-statistics with pre-calculated neighborhoods
 *
 * Count the number of edges, or sum the edge weigths in
 * neighborhoods given as a parameter.
 *
 * \param graph The graph to perform the counting/summing in.
 * \param res Initialized vector, the result is stored here.
 * \param weights Weight vector for weighted graphs, null pointer for
 *        unweighted graphs.
 * \param neighborhoods List of <code>igraph_vector_int_t</code>
 *        objects, the neighborhoods, one for each vertex in the
 *        graph.
 * \return Error code.
 */

int igraph_local_scan_neighborhood_ecount(const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vector_t *weights,
        const igraph_vector_ptr_t *neighborhoods) {

    int node, no_of_nodes = igraph_vcount(graph);
    igraph_inclist_t incs;
    igraph_vector_int_t marked;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (weights && igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length in local scan", IGRAPH_EINVAL);
    }
    if (igraph_vector_ptr_size(neighborhoods) != no_of_nodes) {
        IGRAPH_ERROR("Invalid neighborhood list length in local scan",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_init(&marked, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &marked);
    IGRAPH_CHECK(igraph_inclist_init(graph, &incs, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (node = 0; node < no_of_nodes; node++) {
        igraph_vector_int_t *nei = VECTOR(*neighborhoods)[node];
        int i, neilen = igraph_vector_int_size(nei);
        VECTOR(marked)[node] = node + 1;
        for (i = 0; i < neilen; i++) {
            int vertex = VECTOR(*nei)[i];
            if (vertex < 0 || vertex >= no_of_nodes) {
                IGRAPH_ERROR("Invalid vertex id in neighborhood list in local scan",
                             IGRAPH_EINVAL);
            }
            VECTOR(marked)[vertex] = node + 1;
        }

        for (i = 0; i < neilen; i++) {
            int vertex = VECTOR(*nei)[i];
            igraph_vector_int_t *edges = igraph_inclist_get(&incs, vertex);
            int j, edgeslen = igraph_vector_int_size(edges);
            for (j = 0; j < edgeslen; j++) {
                int edge = VECTOR(*edges)[j];
                int nei2 = IGRAPH_OTHER(graph, edge, vertex);
                if (VECTOR(marked)[nei2] == node + 1) {
                    igraph_real_t w = weights ? VECTOR(*weights)[edge] : 1;
                    VECTOR(*res)[node] += w;
                }
            }
        }
        if (!directed) {
            VECTOR(*res)[node] /= 2.0;
        }
    }

    igraph_inclist_destroy(&incs);
    igraph_vector_int_destroy(&marked);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}
