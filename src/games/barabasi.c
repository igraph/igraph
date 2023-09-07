/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_conversion.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_psumtree.h"
#include "igraph_random.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

/* Attraction function for barabasi_game.
 * We special-case power == 0 to ensure that 0^0 is computed as 1 instead of NaN. */
static igraph_real_t attraction(igraph_real_t degree, igraph_real_t power, igraph_real_t A) {
    return ( power == 0 ? 1.0 : pow(degree, power) ) + A;
}

static igraph_error_t igraph_i_barabasi_game_bag(igraph_t *graph, igraph_integer_t n,
                                      igraph_integer_t m,
                                      const igraph_vector_int_t *outseq,
                                      igraph_bool_t outpref,
                                      igraph_bool_t directed,
                                      const igraph_t *start_from);

static igraph_error_t igraph_i_barabasi_game_psumtree_multiple(igraph_t *graph,
                                                    igraph_integer_t n,
                                                    igraph_real_t power,
                                                    igraph_integer_t m,
                                                    const igraph_vector_int_t *outseq,
                                                    igraph_bool_t outpref,
                                                    igraph_real_t A,
                                                    igraph_bool_t directed,
                                                    const igraph_t *start_from);

static igraph_error_t igraph_i_barabasi_game_psumtree(igraph_t *graph,
                                           igraph_integer_t n,
                                           igraph_real_t power,
                                           igraph_integer_t m,
                                           const igraph_vector_int_t *outseq,
                                           igraph_bool_t outpref,
                                           igraph_real_t A,
                                           igraph_bool_t directed,
                                           const igraph_t *start_from);

static igraph_error_t igraph_i_barabasi_game_bag(igraph_t *graph, igraph_integer_t n,
                                      igraph_integer_t m,
                                      const igraph_vector_int_t *outseq,
                                      igraph_bool_t outpref,
                                      igraph_bool_t directed,
                                      const igraph_t *start_from) {

    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_neighbors = m;
    igraph_integer_t *bag;
    igraph_integer_t bagp = 0;
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t resp;
    igraph_integer_t i, j, k;
    igraph_integer_t bagsize, start_nodes, start_edges, new_edges, no_of_edges;

    if (!directed) {
        outpref = true;
    }

    start_nodes = start_from ? igraph_vcount(start_from) : 1;
    start_edges = start_from ? igraph_ecount(start_from) : 0;
    if (outseq) {
        if (igraph_vector_int_size(outseq) > 1) {
            IGRAPH_CHECK(igraph_i_safe_vector_int_sum(outseq, &new_edges));
            new_edges -= VECTOR(*outseq)[0];
        } else {
            new_edges = 0;
        }
    } else {
        IGRAPH_SAFE_MULT(no_of_nodes - start_nodes, no_of_neighbors, &new_edges);
    }
    IGRAPH_SAFE_ADD(start_edges, new_edges, &no_of_edges);
    /* To ensure the size of the edges vector will not overflow. */
    if (no_of_edges > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Overflow in number of edges.", IGRAPH_EOVERFLOW);
    }
    resp = start_edges * 2;
    bagsize = no_of_nodes;
    IGRAPH_SAFE_ADD(bagsize, no_of_edges, &bagsize);
    if (outpref) {
        IGRAPH_SAFE_ADD(bagsize, no_of_edges, &bagsize);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);

    bag = IGRAPH_CALLOC(bagsize, igraph_integer_t);
    if (bag == 0) {
        IGRAPH_ERROR("barabasi_game failed", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, bag);

    /* The first node(s) in the bag */
    if (start_from) {
        igraph_vector_int_t deg;
        igraph_integer_t ii, jj, sn = igraph_vcount(start_from);
        igraph_neimode_t mm = outpref ? IGRAPH_ALL : IGRAPH_IN;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&deg, sn);
        IGRAPH_CHECK(igraph_degree(start_from, &deg, igraph_vss_all(), mm,
                                   IGRAPH_LOOPS));
        for (ii = 0; ii < sn; ii++) {
            igraph_integer_t d = VECTOR(deg)[ii];
            for (jj = 0; jj <= d; jj++) {
                bag[bagp++] = ii;
            }
        }

        igraph_vector_int_destroy(&deg);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        bag[bagp++] = 0;
    }

    /* Initialize the edges vector */
    if (start_from) {
        IGRAPH_CHECK(igraph_get_edgelist(start_from, &edges, /* bycol= */ false));
        igraph_vector_int_resize(&edges, no_of_edges * 2);
    }

    RNG_BEGIN();

    /* and the others */

    for (i = (start_from ? start_nodes : 1), k = (start_from ? 0 : 1);
         i < no_of_nodes; i++, k++) {

        IGRAPH_ALLOW_INTERRUPTION();

        /* draw edges */
        if (outseq) {
            no_of_neighbors = VECTOR(*outseq)[k];
        }
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_integer_t to = bag[RNG_INTEGER(0, bagp - 1)];
            VECTOR(edges)[resp++] = i;
            VECTOR(edges)[resp++] = to;
        }
        /* update bag */
        bag[bagp++] = i;
        for (j = 0; j < no_of_neighbors; j++) {
            bag[bagp++] = VECTOR(edges)[resp - 2 * j - 1];
            if (outpref) {
                bag[bagp++] = i;
            }
        }
    }

    RNG_END();

    IGRAPH_FREE(bag);
    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_barabasi_game_psumtree_multiple(igraph_t *graph,
                                                    igraph_integer_t n,
                                                    igraph_real_t power,
                                                    igraph_integer_t m,
                                                    const igraph_vector_int_t *outseq,
                                                    igraph_bool_t outpref,
                                                    igraph_real_t A,
                                                    igraph_bool_t directed,
                                                    const igraph_t *start_from) {

    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_neighbors = m;
    igraph_vector_int_t edges;
    igraph_integer_t i, j, k;
    igraph_psumtree_t sumtree;
    igraph_integer_t edgeptr = 0;
    igraph_vector_int_t degree;
    igraph_integer_t start_nodes, start_edges, new_edges, no_of_edges;

    if (!directed) {
        outpref = true;
    }

    start_nodes = start_from ? igraph_vcount(start_from) : 1;
    start_edges = start_from ? igraph_ecount(start_from) : 0;
    if (outseq) {
        if (igraph_vector_int_size(outseq) > 1) {
            IGRAPH_CHECK(igraph_i_safe_vector_int_sum(outseq, &new_edges));
            new_edges -= VECTOR(*outseq)[0];
        } else {
            new_edges = 0;
        }
    } else {
        IGRAPH_SAFE_MULT(no_of_nodes - start_nodes, no_of_neighbors, &new_edges);
    }
    IGRAPH_SAFE_ADD(start_edges, new_edges, &no_of_edges);
    /* To ensure the size of the edges vector will not overflow. */
    if (no_of_edges > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Overflow in number of edges.", IGRAPH_EOVERFLOW);
    }
    edgeptr = start_edges * 2;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, no_of_nodes);

    /* First node(s): */
    if (start_from) {
        igraph_integer_t ii, sn = igraph_vcount(start_from);
        igraph_neimode_t mm = outpref ? IGRAPH_ALL : IGRAPH_IN;
        IGRAPH_CHECK(igraph_degree(start_from, &degree, igraph_vss_all(), mm,
                                   IGRAPH_LOOPS));
        IGRAPH_CHECK(igraph_vector_int_resize(&degree,  no_of_nodes));
        for (ii = 0; ii < sn; ii++) {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, ii, attraction(VECTOR(degree)[ii], power, A)));
        }
    } else {
        /* Any weight may be used for the first node. In the first step, it will be connected to
         * with certainty, after which its weight will be set appropriately. */
        IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, 1.0));
    }

    /* Initialize the edges vector */
    if (start_from) {
        IGRAPH_CHECK(igraph_get_edgelist(start_from, &edges, /* bycol= */ false));
        igraph_vector_int_resize(&edges, no_of_edges * 2);
    }

    RNG_BEGIN();

    /* And the rest: */
    for (i = (start_from ? start_nodes : 1), k = (start_from ? 0 : 1);
         i < no_of_nodes; i++, k++) {
        igraph_real_t sum = igraph_psumtree_sum(&sumtree);
        igraph_integer_t to;

        IGRAPH_ALLOW_INTERRUPTION();

        if (outseq) {
            no_of_neighbors = VECTOR(*outseq)[k];
        }
        for (j = 0; j < no_of_neighbors; j++) {
            if (sum == 0) {
                /* If none of the so-far added nodes have positive weights,
                 * we choose one uniformly to connect to. */
                to = RNG_INTEGER(0, i-1);
            } else {
                igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            }
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
        }
        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_integer_t nn = VECTOR(edges)[edgeptr - 2 * j - 1];
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, nn, attraction(VECTOR(degree)[nn], power, A)));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, attraction(VECTOR(degree)[i], power, A)));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, attraction(0, power, A)));
        }
    }

    RNG_END();

    igraph_psumtree_destroy(&sumtree);
    igraph_vector_int_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_barabasi_game_psumtree(igraph_t *graph,
                                           igraph_integer_t n,
                                           igraph_real_t power,
                                           igraph_integer_t m,
                                           const igraph_vector_int_t *outseq,
                                           igraph_bool_t outpref,
                                           igraph_real_t A,
                                           igraph_bool_t directed,
                                           const igraph_t *start_from) {

    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_neighbors = m;
    igraph_vector_int_t edges;
    igraph_integer_t i, j, k;
    igraph_psumtree_t sumtree;
    igraph_integer_t edgeptr = 0;
    igraph_vector_int_t degree;
    igraph_integer_t start_nodes, start_edges, new_edges, no_of_edges;

    if (!directed) {
        outpref = true;
    }

    start_nodes = start_from ? igraph_vcount(start_from) : 1;
    start_edges = start_from ? igraph_ecount(start_from) : 0;
    if (outseq) {
        if (igraph_vector_int_size(outseq) > 1) {
            IGRAPH_CHECK(igraph_i_safe_vector_int_sum(outseq, &new_edges));
            new_edges -= VECTOR(*outseq)[0];
        } else {
            new_edges = 0;
        }
    } else {
        IGRAPH_SAFE_MULT(no_of_nodes - start_nodes, no_of_neighbors, &new_edges);
    }
    IGRAPH_SAFE_ADD(start_edges, new_edges, &no_of_edges);
    /* To ensure the size of the edges vector will not overflow. */
    if (no_of_edges > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Overflow in number of edges.", IGRAPH_EOVERFLOW);
    }
    edgeptr = start_edges * 2;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, no_of_nodes);

    RNG_BEGIN();

    /* First node(s): */
    if (start_from) {
        igraph_integer_t ii, sn = igraph_vcount(start_from);
        igraph_neimode_t mm = outpref ? IGRAPH_ALL : IGRAPH_IN;
        IGRAPH_CHECK(igraph_degree(start_from, &degree, igraph_vss_all(), mm,
                                   IGRAPH_LOOPS));
        IGRAPH_CHECK(igraph_vector_int_resize(&degree,  no_of_nodes));
        for (ii = 0; ii < sn; ii++) {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, ii, attraction(VECTOR(degree)[ii], power, A)));
        }
    } else {
        /* Any weight may be used for the first node. In the first step, it will be connected to
         * with certainty, after which its weight will be set appropriately. */
        IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, 1.0));
    }

    /* Initialize the edges vector */
    if (start_from) {
        IGRAPH_CHECK(igraph_get_edgelist(start_from, &edges, /* bycol= */ false));
    }

    /* And the rest: */
    for (i = (start_from ? start_nodes : 1), k = (start_from ? 0 : 1);
         i < no_of_nodes; i++, k++) {
        igraph_real_t sum;
        igraph_integer_t to;

        IGRAPH_ALLOW_INTERRUPTION();

        if (outseq) {
            no_of_neighbors = VECTOR(*outseq)[k];
        }
        if (no_of_neighbors >= i) {
            /* All existing vertices are cited */
            for (to = 0; to < i; to++) {
                VECTOR(degree)[to]++;
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                edgeptr += 2;
                IGRAPH_CHECK(igraph_psumtree_update(&sumtree, to, attraction(VECTOR(degree)[to], power, A)));
            }
        } else {
            for (j = 0; j < no_of_neighbors; j++) {
                sum = igraph_psumtree_sum(&sumtree);
                if (sum == 0) {
                    /* If none of the so-far added nodes have positive weights,
                     * we choose one uniformly to connect to. */
                    to = RNG_INTEGER(0, i-1);
                } else {
                    igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
                }
                VECTOR(degree)[to]++;
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                edgeptr += 2;
                IGRAPH_CHECK(igraph_psumtree_update(&sumtree, to, 0.0));
            }
            /* update probabilities */
            for (j = 0; j < no_of_neighbors; j++) {
                igraph_integer_t nn = VECTOR(edges)[edgeptr - 2 * j - 1];
                IGRAPH_CHECK(igraph_psumtree_update(&sumtree, nn, attraction(VECTOR(degree)[nn], power, A)));
            }
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors > i ? i : no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, attraction(VECTOR(degree)[i], power, A)));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, attraction(0, power, A)));
        }
    }

    RNG_END();

    igraph_psumtree_destroy(&sumtree);
    igraph_vector_int_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_barabasi_game
 * \brief Generates a graph based on the Barab&aacute;si-Albert model.
 *
 * This function implements several variants of the preferential attachment
 * process, including linear and non-linear varieties of the Barabási-Albert
 * and Price models. The graph construction starts with a single vertex,
 * or an existing graph given by the \p start_from parameter. Then new vertices
 * are added one at a time. Each new vertex connects to \p m existing vertices,
 * choosing them with probabilities proportional to
 *
 * </para><para>
 * <code>d^power + A</code>,
 *
 * </para><para>
 * where \c d is the in- or total degree of the existing vertex (controlled
 * by the \p outpref argument), while \p power and \p A are given by
 * parameters. The <emphasis>constant attractiveness</emphasis> \p A
 * is used to ensure that vertices with zero in-degree can also be
 * connected to with non-zero probability.
 *
 * </para><para>
 * Barabási, A.-L. and Albert R. 1999. Emergence of scaling in
 * random networks, Science, 286 509--512.
 * https://doi.org/10.1126/science.286.5439.509
 *
 * </para><para>
 * de Solla Price, D. J. 1965. Networks of Scientific Papers, Science,
 * 149 510--515.
 * https://doi.org/10.1126/science.149.3683.510
 *
 * \param graph An uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param power Power of the preferential attachment. In the classic preferential
 *        attachment model <code>power=1</code>. Other values allow for
 *        sampling from a non-linear preferential attachment model.
 *        Negative values are only allowed when no zero-degree vertices
 *        are present during the construction process, i.e. when
 *        the starting graph has no isolated vertices and \p outpref
 *        is set to \c true.
 * \param m The number of outgoing edges generated for each
 *        vertex. Only used when \p outseq is \c NULL.
 * \param outseq Gives the (out-)degrees of the vertices. If this is
 *        constant, this can be a \c NULL pointer or an empty vector.
 *        In this case \p m contains the constant out-degree.
 *        The very first vertex has by definition no outgoing edges,
 *        so the first number in this vector is ignored.
 * \param outpref Boolean, if true not only the in- but also the out-degree
 *        of a vertex increases its citation probability. I.e., the
 *        citation probability is determined by the total degree of
 *        the vertices. Ignored and assumed to be true if the graph
 *        being generated is undirected.
 * \param A The constant attractiveness of vertices. When \p outpref
 *        is set to \c false, it should be positive to ensure that
 *        zero in-degree vertices can be connected to as well.
 * \param directed Boolean, whether to generate a directed graph.
 *        When set to \c false, outpref is assumed to be \c true.
 * \param algo The algorithm to use to generate the network. Possible
 *        values:
 *        \clist
 *        \cli IGRAPH_BARABASI_BAG
 *          This is the algorithm that was previously (before version
 *          0.6) solely implemented in igraph. It works by putting the
 *          IDs of the vertices into a bag (multiset, really), exactly
 *          as many times as their (in-)degree, plus once more. Then
 *          the required number of cited vertices are drawn from the
 *          bag, with replacement. This method might generate multiple
 *          edges. It only works if power=1 and A=1.
 *        \cli IGRAPH_BARABASI_PSUMTREE
 *          This algorithm uses a partial prefix-sum tree to generate
 *          the graph. It does not generate multiple edges and
 *          works for any power and A values.
 *        \cli IGRAPH_BARABASI_PSUMTREE_MULTIPLE
 *          This algorithm also uses a partial prefix-sum tree to
 *          generate the graph. The difference is, that now multiple
 *          edges are allowed. This method was implemented under the
 *          name \c igraph_nonlinear_barabasi_game before version 0.6.
 *        \endclist
 * \param start_from Either a \c NULL pointer, or a graph. In the former
 *        case, the starting configuration is a clique of size \p m.
 *        In the latter case, the graph is a starting configuration.
 *        The graph must be non-empty, i.e. it must have at least one
 *        vertex. If a graph is supplied here and the \p outseq
 *        argument is also given, then \p outseq should only contain
 *        information on the vertices that are not in the \p
 *        start_from graph.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid \p n, \p m, \p A or \p outseq parameter.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges.
 *
 * \example examples/simple/igraph_barabasi_game.c
 * \example examples/simple/igraph_barabasi_game2.c
 */
igraph_error_t igraph_barabasi_game(igraph_t *graph, igraph_integer_t n,
                         igraph_real_t power,
                         igraph_integer_t m,
                         const igraph_vector_int_t *outseq,
                         igraph_bool_t outpref,
                         igraph_real_t A,
                         igraph_bool_t directed,
                         igraph_barabasi_algorithm_t algo,
                         const igraph_t *start_from) {

    igraph_integer_t start_nodes = start_from ? igraph_vcount(start_from) : 0;
    igraph_integer_t newn = start_from ? n - start_nodes : n;

    /* Fix obscure parameterizations */
    if (outseq && igraph_vector_int_empty(outseq)) {
        outseq = NULL;
    }
    if (!directed) {
        outpref = true;
    }

    /* Check arguments */
    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    } else if (newn < 0) {
        IGRAPH_ERROR("Starting graph has too many vertices.", IGRAPH_EINVAL);
    }
    if (start_from && start_nodes == 0) {
        IGRAPH_ERROR("Cannot start from an empty graph.", IGRAPH_EINVAL);
    }
    if (outseq && igraph_vector_int_size(outseq) != newn) {
        IGRAPH_ERROR("Invalid out-degree sequence length.", IGRAPH_EINVAL);
    }
    if (!outseq && m < 0) {
        IGRAPH_ERROR("Number of edges added per step must not be negative.", IGRAPH_EINVAL);
    }
    if (outseq && igraph_vector_int_min(outseq) < 0) {
        IGRAPH_ERROR("Negative out-degree in sequence.", IGRAPH_EINVAL);
    }
    if (!outpref && A <= 0) {
        IGRAPH_ERROR("Constant attractiveness (A) must be positive.",
                     IGRAPH_EINVAL);
    }
    if (outpref && A < 0) {
        IGRAPH_ERROR("Constant attractiveness (A) must be non-negative.",
                     IGRAPH_EINVAL);
    }
    if (algo == IGRAPH_BARABASI_BAG) {
        if (power != 1) {
            IGRAPH_ERROR("Power must be one for bag algorithm.", IGRAPH_EINVAL);
        }
        if (A != 1) {
            IGRAPH_ERROR("Constant attractiveness (A) must be one for bag algorithm.",
                         IGRAPH_EINVAL);
        }
    }
    if (start_from && directed != igraph_is_directed(start_from)) {
        IGRAPH_WARNING("Directedness of the start graph and the output graph mismatch.");
    }
    if (start_from && !igraph_is_directed(start_from) && !outpref) {
        IGRAPH_ERROR("`outpref' must be true if starting from an undirected graph.",
                     IGRAPH_EINVAL);
    }

    if (n == 0) {
        return igraph_empty(graph, 0, directed);
    }

    switch (algo) {
    case IGRAPH_BARABASI_BAG:
        return igraph_i_barabasi_game_bag(graph, n, m, outseq, outpref, directed, start_from);

    case IGRAPH_BARABASI_PSUMTREE:
        return igraph_i_barabasi_game_psumtree(graph, n, power, m, outseq,
                                               outpref, A, directed, start_from);
    case IGRAPH_BARABASI_PSUMTREE_MULTIPLE:
        return igraph_i_barabasi_game_psumtree_multiple(graph, n, power, m,
                                                        outseq, outpref, A,
                                                        directed, start_from);
    default:
        IGRAPH_ERROR("Invalid algorithm for Barabasi game.", IGRAPH_EINVAL);
    }
}

/* Attraction function for barabasi_aging_game.
 * We special-case deg_exp == 0 to ensure that 0^0 is computed as 1 instead of NaN. */
static igraph_real_t attraction_aging(
        igraph_real_t deg, igraph_real_t age,
        igraph_real_t deg_exp, igraph_real_t age_exp,
        igraph_real_t deg_A, igraph_real_t age_A,
        igraph_real_t deg_coef, igraph_real_t age_coef) {

    igraph_real_t dp = deg_exp == 0 ? 1.0 : pow(deg, deg_exp);
    igraph_real_t ap = pow(age, age_exp);
    return (deg_coef * dp + deg_A) * (age_coef * ap + age_A);
}

/**
 * \function igraph_barabasi_aging_game
 * \brief Preferential attachment with aging of vertices.
 *
 * </para><para>
 * This game starts with one vertex (if \p nodes > 0). In each step
 * a new node is added, and it is connected to \p m existing nodes.
 * Existing nodes to connect to are chosen with probability dependent
 * on their (in-)degree (\c k) and age (\c l).
 * The degree-dependent part is
 * <code>deg_coef * k^pa_exp + zero_deg_appeal</code>,
 * while the age-dependent part is
 * <code>age_coef * l^aging_exp + zero_age_appeal</code>,
 * which are multiplied to obtain the final weight.
 *
 * </para><para>
 * The age \c l is based on the number of vertices in the
 * network and the \p aging_bins argument: the age of a node
 * is incremented by 1 after each
 * <code>floor(nodes / aging_bins) + 1</code>
 * time steps.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the graph.
 * \param m The number of edges to add in each time step.
 *        Ignored if \p outseq is a non-zero length vector.
 * \param outseq The number of edges to add in each time step. If it
 *        is \c NULL or a zero-length vector then it is ignored
 *        and the \p m argument is used instead.
 * \param outpref Logical constant, whether the edges
 *        initiated by a vertex contribute to the probability to gain
 *        a new edge.
 * \param pa_exp The exponent of the preferential attachment, a small
 *        positive number usually, the value 1 yields the classic
 *        linear preferential attachment.
 * \param aging_exp The exponent of the aging, this is a negative
 *        number usually.
 * \param aging_bins Integer constant, the number of age bins to use.
 * \param zero_deg_appeal The degree dependent part of the
 *        attractiveness of the zero degree vertices.
 * \param zero_age_appeal The age dependent part of the attractiveness
 *        of the vertices of age zero. This parameter is usually zero.
 * \param deg_coef The coefficient for the degree.
 * \param age_coef The coefficient for the age.
 * \param directed Logical constant, whether to generate a directed
 *        graph.
 * \return Error code.
 *
 * Time complexity: O((|V|+|V|/aging_bins)*log(|V|)+|E|). |V| is the number
 * of vertices, |E| the number of edges.
 */
igraph_error_t igraph_barabasi_aging_game(igraph_t *graph,
                               igraph_integer_t nodes,
                               igraph_integer_t m,
                               const igraph_vector_int_t *outseq,
                               igraph_bool_t outpref,
                               igraph_real_t pa_exp,
                               igraph_real_t aging_exp,
                               igraph_integer_t aging_bins,
                               igraph_real_t zero_deg_appeal,
                               igraph_real_t zero_age_appeal,
                               igraph_real_t deg_coef,
                               igraph_real_t age_coef,
                               igraph_bool_t directed) {
    igraph_integer_t no_of_nodes = nodes;
    igraph_integer_t no_of_neighbors = m;
    igraph_integer_t binwidth;
    igraph_integer_t no_of_edges;
    igraph_vector_int_t edges;
    igraph_integer_t i, j, k;
    igraph_psumtree_t sumtree;
    igraph_integer_t edgeptr = 0;
    igraph_vector_int_t degree;

    if (no_of_nodes < 0) {
        IGRAPH_ERRORF("Number of nodes must not be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, no_of_nodes);
    }
    if (outseq != 0 && igraph_vector_int_size(outseq) != 0 && igraph_vector_int_size(outseq) != no_of_nodes) {
        IGRAPH_ERRORF("The length of the out-degree sequence (%" IGRAPH_PRId ") does not agree with the number of nodes (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_int_size(outseq), no_of_nodes);
    }
    if ( (outseq == 0 || igraph_vector_int_size(outseq) == 0) && m < 0) {
        IGRAPH_ERRORF("The number of edges per time step must not be negative, got %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL,
                      m);
    }
    if (aging_bins <= 0) {
        IGRAPH_ERRORF("Number of aging bins must be positive, got %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL,
                      aging_bins);
    }
    if (deg_coef < 0) {
        IGRAPH_ERRORF("Degree coefficient must be non-negative, got %g.",
                      IGRAPH_EINVAL,
                      deg_coef);
    }
    if (age_coef < 0) {
        IGRAPH_ERRORF("Age coefficient must be non-negative, got %g.",
                      IGRAPH_EINVAL,
                      deg_coef);
    }

    if (zero_deg_appeal < 0) {
        IGRAPH_ERRORF("Zero degree appeal must be non-negative, got %g.",
                      IGRAPH_EINVAL,
                      zero_deg_appeal);
    }
    if (zero_age_appeal < 0) {
        IGRAPH_ERRORF("Zero age appeal must be non-negative, got %g.",
                      IGRAPH_EINVAL,
                      zero_age_appeal);
    }

    if (no_of_nodes == 0) {
         return igraph_empty(graph, 0, directed);
    }

    binwidth = no_of_nodes / aging_bins + 1;

    if (outseq == 0 || igraph_vector_int_size(outseq) == 0) {
        no_of_neighbors = m;
        IGRAPH_SAFE_MULT(no_of_nodes - 1, no_of_neighbors, &no_of_edges);
    } else {
        IGRAPH_CHECK(igraph_i_safe_vector_int_sum(outseq, &no_of_edges));
        no_of_edges -= VECTOR(*outseq)[0];
    }
    /* To ensure the size of the edges vector will not overflow. */
    if (no_of_edges > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Overflow in number of edges.", IGRAPH_EOVERFLOW);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, no_of_nodes);

    RNG_BEGIN();

    /* First node: */
    /* Any weight may be used for the first node. In the first step, it will be connected to
     * with certainty, after which its weight will be set appropriately. */
    IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, 1.0));

    /* And the rest: */
    for (i = 1; i < no_of_nodes; i++) {
        igraph_real_t sum;
        igraph_integer_t to;

        IGRAPH_ALLOW_INTERRUPTION();

        if (outseq != 0 && igraph_vector_int_size(outseq) != 0) {
            no_of_neighbors = VECTOR(*outseq)[i];
        }
        sum = igraph_psumtree_sum(&sumtree);
        for (j = 0; j < no_of_neighbors; j++) {
            if (sum == 0) {
                /* If none of the so-far added nodes have positive weights,
                 * we choose one uniformly to connect to. */
                to = RNG_INTEGER(0, i-1);
            } else {
                igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            }
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
        }
        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_integer_t n = VECTOR(edges)[edgeptr - 2 * j - 1];
            igraph_integer_t age = (i - n) / binwidth;
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, n,
                attraction_aging(VECTOR(degree)[n], age+1,
                                 pa_exp, aging_exp,
                                 zero_deg_appeal, zero_age_appeal,
                                 deg_coef, age_coef)
            ));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, i,
                 attraction_aging(VECTOR(degree)[i], 1,
                                  pa_exp, aging_exp,
                                  zero_deg_appeal, zero_age_appeal,
                                  deg_coef, age_coef)
            ));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, i,
                 attraction_aging(0, 1,
                                  pa_exp, aging_exp,
                                  zero_deg_appeal, zero_age_appeal,
                                  deg_coef, age_coef)
            ));
        }

        /* aging */
        for (k = 1; binwidth * k <= i; k++) {
            igraph_integer_t shnode = i - binwidth * k;
            igraph_integer_t deg = VECTOR(degree)[shnode];
            igraph_integer_t age = (i - shnode) / binwidth;
            /* igraph_real_t old=igraph_psumtree_get(&sumtree, shnode); */
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, shnode,
                 attraction_aging(deg, age + 2,
                                  pa_exp, aging_exp,
                                  zero_deg_appeal, zero_age_appeal,
                                  deg_coef, age_coef)
                ));
        }
    }

    RNG_END();

    igraph_vector_int_destroy(&degree);
    igraph_psumtree_destroy(&sumtree);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
