/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_transitivity.h"

#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_memory.h"
#include "igraph_motifs.h"
#include "igraph_structural.h"

#include "core/interruption.h"
#include "properties/properties_internal.h"

/**
 * \function igraph_transitivity_avglocal_undirected
 * \brief Average local transitivity (clustering coefficient).
 *
 * The transitivity measures the probability that two neighbors of a
 * vertex are connected. In case of the average local transitivity,
 * this probability is calculated for each vertex and then the average
 * is taken. Vertices with less than two neighbors require special treatment,
 * they will either be left out from the calculation or they will be considered
 * as having zero transitivity, depending on the \c mode argument.
 * Edge directions and edge multiplicities are ignored.
 *
 * </para><para>
 * Note that this measure is different from the global transitivity measure
 * (see \ref igraph_transitivity_undirected() ) as it simply takes the
 * average local transitivity across the whole network.
 *
 * </para><para>
 * Clustering coefficient is an alternative name for transitivity.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * D. J. Watts and S. Strogatz: Collective dynamics of small-world networks.
 * Nature 393(6684):440-442 (1998).
 *
 * \param graph The input graph. Edge directions and multiplicites are ignored.
 * \param res Pointer to a real variable, the result will be stored here.
 * \param mode Defines how to treat vertices with degree less than two.
 *    \c IGRAPH_TRANSITIVITY_NAN leaves them out from averaging,
 *    \c IGRAPH_TRANSITIVITY_ZERO includes them with zero transitivity.
 *    The result will be \c NaN if the mode is \c IGRAPH_TRANSITIVITY_NAN
 *    and there are no vertices with more than one neighbor.
 *
 * \return Error code.
 *
 * \sa \ref igraph_transitivity_undirected(), \ref
 * igraph_transitivity_local_undirected().
 *
 * Time complexity: O(|V|*d^2), |V| is the number of vertices in the
 * graph and d is the average degree.
 */

igraph_error_t igraph_transitivity_avglocal_undirected(const igraph_t *graph,
        igraph_real_t *res,
        igraph_transitivity_mode_t mode) {

    igraph_integer_t i, no_of_nodes = igraph_vcount(graph), nans = 0;
    igraph_real_t sum = 0.0;
    igraph_vector_t vec;

    if (no_of_nodes == 0) {
        if (mode == IGRAPH_TRANSITIVITY_ZERO) {
            *res = 0;
        } else {
            *res = IGRAPH_NAN;
        }
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&vec, no_of_nodes);

        IGRAPH_CHECK(igraph_transitivity_local_undirected(graph, &vec, igraph_vss_all(), mode));

        for (i = 0, nans = 0; i < no_of_nodes; i++) {
            if (!isnan(VECTOR(vec)[i])) {
                sum += VECTOR(vec)[i];
            } else {
                nans++;
            }
        }

        igraph_vector_destroy(&vec);
        IGRAPH_FINALLY_CLEAN(1);

        *res = sum / (no_of_nodes - nans);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_transitivity_local_undirected1(const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vs_t vids,
        igraph_transitivity_mode_t mode) {

#define TRANSIT
#include "properties/triangles_template1.h"
#undef TRANSIT

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_transitivity_local_undirected2(const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vs_t vids,
        igraph_transitivity_mode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vit_t vit;
    igraph_integer_t nodes_to_calc, affected_nodes;
    igraph_integer_t maxdegree = 0;
    igraph_integer_t i, j, k, nn;
    igraph_lazy_adjlist_t adjlist;
    igraph_vector_int_t degree;
    igraph_vector_t indexv, avids, rank, triangles;
    igraph_vector_int_t order;
    igraph_integer_t *neis;

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    IGRAPH_VECTOR_INIT_FINALLY(&indexv, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&avids, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&avids, nodes_to_calc));
    k = 0;
    for (i = 0; i < nodes_to_calc; IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t v = IGRAPH_VIT_GET(vit);
        igraph_vector_int_t *neis2;
        igraph_integer_t neilen;
        if (VECTOR(indexv)[v] == 0) {
            VECTOR(indexv)[v] = k + 1; k++;
            IGRAPH_CHECK(igraph_vector_push_back(&avids, v));
        }

        neis2 = igraph_lazy_adjlist_get(&adjlist, v);
        IGRAPH_CHECK_OOM(neis2, "Failed to query neighbors.");

        neilen = igraph_vector_int_size(neis2);
        for (j = 0; j < neilen; j++) {
            igraph_integer_t nei = VECTOR(*neis2)[j];
            if (VECTOR(indexv)[nei] == 0) {
                VECTOR(indexv)[nei] = k + 1; k++;
                IGRAPH_CHECK(igraph_vector_push_back(&avids, nei));
            }
        }
    }

    /* Degree, ordering, ranking */
    affected_nodes = igraph_vector_size(&avids);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, affected_nodes);
    for (i = 0; i < affected_nodes; i++) {
        igraph_integer_t v = VECTOR(avids)[i];
        igraph_vector_int_t *neis2;
        igraph_integer_t deg;
        neis2 = igraph_lazy_adjlist_get(&adjlist, v);
        IGRAPH_CHECK_OOM(neis2, "Failed to query neighbors.");
        VECTOR(degree)[i] = deg = igraph_vector_int_size(neis2);
        if (deg > maxdegree) {
            maxdegree = deg;
        }
    }
    IGRAPH_CHECK(igraph_vector_int_order1(&degree, &order, maxdegree + 1));
    igraph_vector_int_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_VECTOR_INIT_FINALLY(&rank, affected_nodes);
    for (i = 0; i < affected_nodes; i++) {
        VECTOR(rank)[ VECTOR(order)[i] ] = affected_nodes - i - 1;
    }

    neis = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (neis == 0) {
        IGRAPH_ERROR("Insufficient memory for local transitivity calculation.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, neis);

    IGRAPH_VECTOR_INIT_FINALLY(&triangles, affected_nodes);
    for (nn = affected_nodes - 1; nn >= 0; nn--) {
        igraph_integer_t node = VECTOR(avids) [ VECTOR(order)[nn] ];
        igraph_vector_int_t *neis1, *neis2;
        igraph_integer_t neilen1, neilen2;
        igraph_integer_t nodeindex = VECTOR(indexv)[node];
        igraph_integer_t noderank = VECTOR(rank) [nodeindex - 1];

        IGRAPH_ALLOW_INTERRUPTION();

        neis1 = igraph_lazy_adjlist_get(&adjlist, node);
        IGRAPH_CHECK_OOM(neis1, "Failed to query neighbors.");

        neilen1 = igraph_vector_int_size(neis1);
        for (i = 0; i < neilen1; i++) {
            igraph_integer_t nei = VECTOR(*neis1)[i];
            neis[nei] = node + 1;
        }
        for (i = 0; i < neilen1; i++) {
            igraph_integer_t nei = VECTOR(*neis1)[i];
            igraph_integer_t neiindex = VECTOR(indexv)[nei];
            igraph_integer_t neirank = VECTOR(rank)[neiindex - 1];

            /*       fprintf(stderr, "  nei %li (indexv %li, rank %li)\n", nei, */
            /*        neiindex, neirank); */
            if (neirank > noderank) {
                neis2 = igraph_lazy_adjlist_get(&adjlist, nei);
                IGRAPH_CHECK_OOM(neis2, "Failed to query neighbors.");
                neilen2 = igraph_vector_int_size(neis2);
                for (j = 0; j < neilen2; j++) {
                    igraph_integer_t nei2 = VECTOR(*neis2)[j];
                    igraph_integer_t nei2index = VECTOR(indexv)[nei2];
                    igraph_integer_t nei2rank = VECTOR(rank)[nei2index - 1];
                    /*    fprintf(stderr, "    triple %li %li %li\n", node, nei, nei2); */
                    if (nei2rank < neirank) {
                        continue;
                    }
                    if (neis[nei2] == node + 1) {
                        /*      fprintf(stderr, "    triangle\n"); */
                        VECTOR(triangles) [ nei2index - 1 ] += 1;
                        VECTOR(triangles) [ neiindex - 1 ] += 1;
                        VECTOR(triangles) [ nodeindex - 1 ] += 1;
                    }
                }
            }
        }
    }

    /* Ok, for all affected vertices the number of triangles were counted */

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    IGRAPH_VIT_RESET(vit);
    for (i = 0; i < nodes_to_calc; i++, IGRAPH_VIT_NEXT(vit)) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        igraph_integer_t idx = VECTOR(indexv)[node] - 1;
        igraph_vector_int_t *neis2 = igraph_lazy_adjlist_get(&adjlist, node);
        igraph_integer_t deg;

        IGRAPH_CHECK_OOM(neis2, "Failed to query neighbors.");

        deg = igraph_vector_int_size(neis2);
        if (mode == IGRAPH_TRANSITIVITY_ZERO && deg < 2) {
            VECTOR(*res)[i] = 0.0;
        } else {
            VECTOR(*res)[i] = VECTOR(triangles)[idx] / deg / (deg - 1) * 2.0;
        }
        /*     fprintf(stderr, "%f %f\n", VECTOR(triangles)[idx], triples); */
    }

    igraph_vector_destroy(&triangles);
    igraph_free(neis);
    igraph_vector_destroy(&rank);
    igraph_vector_int_destroy(&order);
    igraph_vector_destroy(&avids);
    igraph_vector_destroy(&indexv);
    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(8);

    return IGRAPH_SUCCESS;
}

/* This removes loop, multiple edges and edges that point
     "backwards" according to the rank vector. */
/* Note: Also used in scan.c */
igraph_error_t igraph_i_trans4_al_simplify(igraph_adjlist_t *al,
                                const igraph_vector_int_t *rank) {
    igraph_integer_t i;
    igraph_integer_t n = al->length;
    igraph_vector_int_t mark;

    IGRAPH_CHECK(igraph_vector_int_init(&mark, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &mark);

    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->adjs[i];
        igraph_integer_t j, l = igraph_vector_int_size(v);
        igraph_integer_t irank = VECTOR(*rank)[i];
        VECTOR(mark)[i] = i + 1;
        for (j = 0; j < l; /* nothing */) {
            igraph_integer_t e = VECTOR(*v)[j];
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

static igraph_error_t igraph_transitivity_local_undirected4(const igraph_t *graph,
        igraph_vector_t *res,
        igraph_transitivity_mode_t mode) {

#define TRANSIT 1
#include "properties/triangles_template.h"
#undef TRANSIT

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_transitivity_local_undirected
 * \brief The local transitivity (clustering coefficient) of some vertices.
 *
 * The transitivity measures the probability that two neighbors of a
 * vertex are connected. In case of the local transitivity, this
 * probability is calculated separately for each vertex.
 *
 * </para><para>
 * Note that this measure is different from the global transitivity measure
 * (see \ref igraph_transitivity_undirected() ) as it calculates a transitivity
 * value for each vertex individually.
 *
 * </para><para>
 * Clustering coefficient is an alternative name for transitivity.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * D. J. Watts and S. Strogatz: Collective dynamics of small-world networks.
 * Nature 393(6684):440-442 (1998).
 *
 * \param graph The input graph. Edge directions and multiplicities are ignored.
 * \param res Pointer to an initialized vector, the result will be
 *   stored here. It will be resized as needed.
 * \param vids Vertex set, the vertices for which the local
 *   transitivity will be calculated.
 * \param mode Defines how to treat vertices with degree less than two.
 *    \c IGRAPH_TRANSITIVITY_NAN returns \c NaN for these vertices,
 *    \c IGRAPH_TRANSITIVITY_ZERO returns zero.
 * \return Error code.
 *
 * \sa \ref igraph_transitivity_undirected(), \ref
 * igraph_transitivity_avglocal_undirected().
 *
 * Time complexity: O(n*d^2), n is the number of vertices for which
 * the transitivity is calculated, d is the average vertex degree.
 */

igraph_error_t igraph_transitivity_local_undirected(const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vs_t vids,
        igraph_transitivity_mode_t mode) {

    if (igraph_vs_is_all(&vids)) {
        return igraph_transitivity_local_undirected4(graph, res, mode);
    } else {
        igraph_vit_t vit;
        igraph_integer_t size;
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        size = IGRAPH_VIT_SIZE(vit);
        igraph_vit_destroy(&vit);
        IGRAPH_FINALLY_CLEAN(1);
        if (size < 100) {
            return igraph_transitivity_local_undirected1(graph, res, vids, mode);
        } else {
            return igraph_transitivity_local_undirected2(graph, res, vids, mode);
        }
    }
}

static igraph_error_t igraph_adjacent_triangles1(const igraph_t *graph,
                                      igraph_vector_t *res,
                                      const igraph_vs_t vids) {
# include "properties/triangles_template1.h"
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_adjacent_triangles4(const igraph_t *graph,
                                      igraph_vector_t *res) {
# include "properties/triangles_template.h"
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_adjacent_triangles
 * \brief Count the number of triangles a vertex is part of.
 *
 * \param graph The input graph. Edge directions and multiplicities are ignored.
 * \param res Initiliazed vector, the results are stored here.
 * \param vids The vertices to perform the calculation for.
 * \return Error mode.
 *
 * \sa \ref igraph_list_triangles() to list them.
 *
 * Time complexity: O(d^2 n), d is the average vertex degree of the
 * queried vertices, n is their number.
 */

igraph_error_t igraph_adjacent_triangles(const igraph_t *graph,
                              igraph_vector_t *res,
                              const igraph_vs_t vids) {
    if (igraph_vs_is_all(&vids)) {
        return igraph_adjacent_triangles4(graph, res);
    } else {
        return igraph_adjacent_triangles1(graph, res, vids);
    }
}

/**
 * \function igraph_list_triangles
 * \brief Find all triangles in a graph.
 *
 * </para><para>
 * The triangles are reported as a long list of vertex ID triplets. Use
 * the \c int variant of \ref igraph_matrix_view_from_vector() to create a
 * matrix view into the vector where each triangle is stored in a column of the
 * matrix (see the example).
 *
 * \param graph The input graph, edge directions are ignored.
 *        Multiple edges are ignored.
 * \param res Pointer to an initialized integer vector, the result
 *        is stored here, in a long list of triples of vertex IDs.
 *        Each triple is a triangle in the graph. Each triangle is
 *        listed exactly once.
 * \return Error code.
 *
 * \sa \ref igraph_transitivity_undirected() to count the triangles,
 * \ref igraph_adjacent_triangles() to count the triangles a vertex
 * participates in.
 *
 * Time complexity: O(d^2 n), d is the average degree, n is the number
 * of vertices.
 *
 * \example examples/simple/igraph_list_triangles.c
 */

igraph_error_t igraph_list_triangles(const igraph_t *graph,
                          igraph_vector_int_t *res) {
# define TRIANGLES
# include "properties/triangles_template.h"
# undef TRIANGLES
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_transitivity_undirected
 * \brief Calculates the transitivity (clustering coefficient) of a graph.
 *
 * </para><para>
 * The transitivity measures the probability that two neighbors of a
 * vertex are connected. More precisely, this is the ratio of the
 * triangles and connected triples in the graph, the result is a
 * single real number. Directed graphs are considered as undirected ones
 * and multi-edges are ignored.
 *
 * </para><para>
 * Note that this measure is different from the local transitivity measure
 * (see \ref igraph_transitivity_local_undirected() ) as it calculates a single
 * value for the whole graph.
 *
 * </para><para>
 * Clustering coefficient is an alternative name for transitivity.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * S. Wasserman and K. Faust: Social Network Analysis: Methods and
 * Applications. Cambridge: Cambridge University Press, 1994.
 *
 * \param graph The graph object. Edge directions and multiplicites are ignored.
 * \param res Pointer to a real variable, the result will be stored here.
 * \param mode Defines how to treat graphs with no connected triples.
 *   \c IGRAPH_TRANSITIVITY_NAN returns \c NaN in this case,
 *   \c IGRAPH_TRANSITIVITY_ZERO returns zero.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory for
 *         temporary data.
 *
 * \sa \ref igraph_transitivity_local_undirected(),
 * \ref igraph_transitivity_avglocal_undirected().
 *
 * Time complexity: O(|V|*d^2), |V| is the number of vertices in
 * the graph, d is the average node degree.
 *
 * \example examples/simple/igraph_transitivity.c
 */

igraph_error_t igraph_transitivity_undirected(const igraph_t *graph,
                                   igraph_real_t *res,
                                   igraph_transitivity_mode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_real_t triples = 0, triangles = 0;
    igraph_integer_t node, nn;
    igraph_integer_t maxdegree;
    igraph_integer_t *neis;
    igraph_vector_int_t order;
    igraph_vector_t rank;
    igraph_vector_int_t degree;

    igraph_adjlist_t allneis;
    igraph_vector_int_t *neis1, *neis2;
    igraph_integer_t i, j, neilen1, neilen2;

    if (no_of_nodes == 0) {
        *res = mode == IGRAPH_TRANSITIVITY_ZERO ? 0.0 : IGRAPH_NAN;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, no_of_nodes);

    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL,
                               IGRAPH_LOOPS));
    maxdegree = igraph_vector_int_max(&degree) + 1;
    IGRAPH_CHECK(igraph_vector_int_order1(&degree, &order, maxdegree));

    igraph_vector_int_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INIT_FINALLY(&rank, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(rank)[ VECTOR(order)[i] ] = no_of_nodes - i - 1;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    neis = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (! neis) {
        IGRAPH_ERROR("Insufficient memory for undirected global transitivity.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, neis);

    for (nn = no_of_nodes - 1; nn >= 0; nn--) {
        node = VECTOR(order)[nn];

        IGRAPH_ALLOW_INTERRUPTION();

        neis1 = igraph_adjlist_get(&allneis, node);
        neilen1 = igraph_vector_int_size(neis1);
        triples += (igraph_real_t)neilen1 * (neilen1 - 1);
        /* Mark the neighbors of 'node' */
        for (i = 0; i < neilen1; i++) {
            igraph_integer_t nei = VECTOR(*neis1)[i];
            neis[nei] = node + 1;
        }
        for (i = 0; i < neilen1; i++) {
            igraph_integer_t nei = VECTOR(*neis1)[i];
            /* If 'nei' is not ready yet */
            if (VECTOR(rank)[nei] > VECTOR(rank)[node]) {
                neis2 = igraph_adjlist_get(&allneis, nei);
                neilen2 = igraph_vector_int_size(neis2);
                for (j = 0; j < neilen2; j++) {
                    igraph_integer_t nei2 = VECTOR(*neis2)[j];
                    if (neis[nei2] == node + 1) {
                        triangles += 1.0;
                    }
                }
            }
        }
    }

    IGRAPH_FREE(neis);
    igraph_adjlist_destroy(&allneis);
    igraph_vector_destroy(&rank);
    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(4);

    if (triples == 0 && mode == IGRAPH_TRANSITIVITY_ZERO) {
        *res = 0;
    } else {
        *res = triangles / triples * 2.0;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_transitivity_barrat1(
        const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vs_t vids,
        const igraph_vector_t *weights,
        igraph_transitivity_mode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vit_t vit;
    igraph_integer_t nodes_to_calc;
    igraph_vector_int_t *adj1, *adj2;
    igraph_vector_int_t neis;
    igraph_vector_t actw;
    igraph_lazy_inclist_t incident;
    igraph_integer_t i;
    igraph_vector_t strength;

    /* Precondition: weight vector is not null, its length equals the number of
     * edges, and the graph has at least one vertex. The graph must not have
     * multi-edges. These must be ensured by the caller.
     */

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_vector_int_init(&neis, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

    IGRAPH_VECTOR_INIT_FINALLY(&actw, no_of_nodes);

    IGRAPH_VECTOR_INIT_FINALLY(&strength, 0);
    IGRAPH_CHECK(igraph_strength(graph, &strength, igraph_vss_all(), IGRAPH_ALL,
                                 IGRAPH_LOOPS, weights));

    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &incident, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &incident);

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        igraph_integer_t adjlen1, adjlen2, j, k;
        igraph_real_t triples, triangles;

        IGRAPH_ALLOW_INTERRUPTION();

        adj1 = igraph_lazy_inclist_get(&incident, node);
        IGRAPH_CHECK_OOM(adj1, "Failed to query incident edges.");
        adjlen1 = igraph_vector_int_size(adj1);
        /* Mark the neighbors of the node */
        for (j = 0; j < adjlen1; j++) {
            igraph_integer_t edge = VECTOR(*adj1)[j];
            igraph_integer_t nei = IGRAPH_OTHER(graph, edge, node);
            VECTOR(neis)[nei] = i + 1;
            VECTOR(actw)[nei] = VECTOR(*weights)[edge];
        }
        triples = VECTOR(strength)[node] * (adjlen1 - 1);
        triangles = 0.0;

        for (j = 0; j < adjlen1; j++) {
            igraph_integer_t edge1 = VECTOR(*adj1)[j];
            igraph_real_t weight1 = VECTOR(*weights)[edge1];
            igraph_integer_t v = IGRAPH_OTHER(graph, edge1, node);
            adj2 = igraph_lazy_inclist_get(&incident, v);
            IGRAPH_CHECK_OOM(adj2, "Failed to query incident edges.");
            adjlen2 = igraph_vector_int_size(adj2);
            for (k = 0; k < adjlen2; k++) {
                igraph_integer_t edge2 = VECTOR(*adj2)[k];
                igraph_integer_t v2 = IGRAPH_OTHER(graph, edge2, v);
                if (VECTOR(neis)[v2] == i + 1) {
                    triangles += (VECTOR(actw)[v2] + weight1) / 2.0;
                }
            }
        }
        if (mode == IGRAPH_TRANSITIVITY_ZERO && triples == 0) {
            VECTOR(*res)[i] = 0.0;
        } else {
            VECTOR(*res)[i] = triangles / triples;
        }
    }

    igraph_lazy_inclist_destroy(&incident);
    igraph_vector_destroy(&strength);
    igraph_vector_destroy(&actw);
    igraph_vector_int_destroy(&neis);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_transitivity_barrat4(
        const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vector_t *weights,
        igraph_transitivity_mode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t order;
    igraph_vector_int_t degree;
    igraph_vector_t strength;
    igraph_vector_t rank;
    igraph_integer_t maxdegree;
    igraph_inclist_t incident;
    igraph_vector_int_t neis;
    igraph_vector_int_t *adj1, *adj2;
    igraph_vector_t actw;
    igraph_integer_t i, nn;

    /* Precondition: weight vector is not null, its length equals the number of
     * edges, and the graph has at least one vertex. The graph must not have
     * multi-edges. These must be ensured by the caller.
     */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&strength, no_of_nodes);

    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL,
                               IGRAPH_LOOPS));
    maxdegree = igraph_vector_int_max(&degree) + 1;
    IGRAPH_CHECK(igraph_vector_int_order1(&degree, &order, maxdegree));

    IGRAPH_CHECK(igraph_strength(graph, &strength, igraph_vss_all(), IGRAPH_ALL,
                                 IGRAPH_LOOPS, weights));

    IGRAPH_VECTOR_INIT_FINALLY(&rank, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(rank)[ VECTOR(order)[i] ] = no_of_nodes - i - 1;
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &incident, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incident);

    IGRAPH_CHECK(igraph_vector_int_init(&neis, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

    IGRAPH_VECTOR_INIT_FINALLY(&actw, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    for (nn = no_of_nodes - 1; nn >= 0; nn--) {
        igraph_integer_t adjlen1, adjlen2;
        igraph_real_t triples;
        igraph_integer_t node = VECTOR(order)[nn];

        IGRAPH_ALLOW_INTERRUPTION();

        adj1 = igraph_inclist_get(&incident, node);
        adjlen1 = igraph_vector_int_size(adj1);
        triples = VECTOR(strength)[node] * (adjlen1 - 1) / 2.0;
        /* Mark the neighbors of the node */
        for (i = 0; i < adjlen1; i++) {
            igraph_integer_t edge = VECTOR(*adj1)[i];
            igraph_integer_t nei = IGRAPH_OTHER(graph, edge, node);
            VECTOR(neis)[nei] = node + 1;
            VECTOR(actw)[nei] = VECTOR(*weights)[edge];
        }

        for (i = 0; i < adjlen1; i++) {
            igraph_integer_t edge1 = VECTOR(*adj1)[i];
            igraph_real_t weight1 = VECTOR(*weights)[edge1];
            igraph_integer_t nei = IGRAPH_OTHER(graph, edge1, node);
            igraph_integer_t j;
            if (VECTOR(rank)[nei] > VECTOR(rank)[node]) {
                adj2 = igraph_inclist_get(&incident, nei);
                adjlen2 = igraph_vector_int_size(adj2);
                for (j = 0; j < adjlen2; j++) {
                    igraph_integer_t edge2 = VECTOR(*adj2)[j];
                    igraph_real_t weight2 = VECTOR(*weights)[edge2];
                    igraph_integer_t nei2 = IGRAPH_OTHER(graph, edge2, nei);
                    if (VECTOR(rank)[nei2] < VECTOR(rank)[nei]) {
                        continue;
                    }
                    if (VECTOR(neis)[nei2] == node + 1) {
                        VECTOR(*res)[nei2] += (VECTOR(actw)[nei2] + weight2) / 2.0;
                        VECTOR(*res)[nei] += (weight1 + weight2) / 2.0;
                        VECTOR(*res)[node] += (VECTOR(actw)[nei2] + weight1) / 2.0;
                    }
                }
            }
        }

        if (mode == IGRAPH_TRANSITIVITY_ZERO && triples == 0) {
            VECTOR(*res)[node] = 0.0;
        } else {
            VECTOR(*res)[node] /= triples;
        }
    }

    igraph_vector_destroy(&actw);
    igraph_vector_int_destroy(&neis);
    igraph_inclist_destroy(&incident);
    igraph_vector_destroy(&rank);
    igraph_vector_int_destroy(&degree);
    igraph_vector_destroy(&strength);
    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_transitivity_barrat
 * \brief Weighted local transitivity of some vertices, as defined by A. Barrat.
 *
 * This is a local transitivity, i.e. a vertex-level index. For a
 * given vertex \c i, from all triangles in which it participates we
 * consider the weight of the edges incident on \c i. The transitivity
 * is the sum of these weights divided by twice the strength of the
 * vertex (see \ref igraph_strength()) and the degree of the vertex
 * minus one. See equation (5) in Alain Barrat, Marc Barthelemy, Romualdo
 * Pastor-Satorras, Alessandro Vespignani: The architecture of complex
 * weighted networks, Proc. Natl. Acad. Sci. USA 101, 3747 (2004) at
 * https://doi.org/10.1073/pnas.0400087101 for the exact formula.
 *
 * \param graph The input graph. Edge directions are ignored for
 *   directed graphs. Note that the function does \em not work for
 *   non-simple graphs.
 * \param res Pointer to an initialized vector, the result will be
 *   stored here. It will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param weights Edge weights. If this is a null pointer, then a
 *   warning is given and \ref igraph_transitivity_local_undirected()
 *   is called.
 * \param mode Defines how to treat vertices with zero strength.
 *   \c IGRAPH_TRANSITIVITY_NAN says that the transitivity of these
 *   vertices is \c NaN, \c IGRAPH_TRANSITIVITY_ZERO says it is zero.
 *
 * \return Error code.
 *
 * Time complexity: O(|V|*d^2), |V| is the number of vertices in
 * the graph, d is the average node degree.
 *
 * \sa \ref igraph_transitivity_undirected(), \ref
 * igraph_transitivity_local_undirected() and \ref
 * igraph_transitivity_avglocal_undirected() for other kinds of
 * (non-weighted) transitivity.
 */

igraph_error_t igraph_transitivity_barrat(const igraph_t *graph,
                               igraph_vector_t *res,
                               const igraph_vs_t vids,
                               const igraph_vector_t *weights,
                               igraph_transitivity_mode_t mode) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t has_multiple;

    /* Handle fallback to unweighted version and common cases */
    if (!weights) {
        if (no_of_edges != 0) {
            IGRAPH_WARNING("No weights given for Barrat's transitivity, unweighted version is used.");
        }
        return igraph_transitivity_local_undirected(graph, res, vids, mode);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Edge weight vector length (%" IGRAPH_PRId ") not equal to "
                      "number of edges (%" IGRAPH_PRId ").", IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (no_of_nodes == 0) {
        igraph_vector_clear(res);
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_has_multiple(graph, &has_multiple));
    if (! has_multiple && igraph_is_directed(graph)) {
        /* When the graph is directed, mutual edges are effectively multi-edges as we
         * are ignoring edge directions. */
        IGRAPH_CHECK(igraph_has_mutual(graph, &has_multiple, false));
    }
    if (has_multiple) {
        IGRAPH_ERROR("Barrat's weighted transitivity measure works only if the graph has no multi-edges.",
                     IGRAPH_EINVAL);
    }

    /* Preconditions validated, now we can call the real implementation */

    if (igraph_vs_is_all(&vids)) {
        return igraph_i_transitivity_barrat4(graph, res, weights, mode);
    } else {
        return igraph_i_transitivity_barrat1(graph, res, vids, weights, mode);
    }
}
