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

static int igraph_i_barabasi_game_bag(igraph_t *graph, igraph_integer_t n,
                                      igraph_integer_t m,
                                      const igraph_vector_t *outseq,
                                      igraph_bool_t outpref,
                                      igraph_bool_t directed,
                                      const igraph_t *start_from);

static int igraph_i_barabasi_game_psumtree_multiple(igraph_t *graph,
                                                    igraph_integer_t n,
                                                    igraph_real_t power,
                                                    igraph_integer_t m,
                                                    const igraph_vector_t *outseq,
                                                    igraph_bool_t outpref,
                                                    igraph_real_t A,
                                                    igraph_bool_t directed,
                                                    const igraph_t *start_from);

static int igraph_i_barabasi_game_psumtree(igraph_t *graph,
                                           igraph_integer_t n,
                                           igraph_real_t power,
                                           igraph_integer_t m,
                                           const igraph_vector_t *outseq,
                                           igraph_bool_t outpref,
                                           igraph_real_t A,
                                           igraph_bool_t directed,
                                           const igraph_t *start_from);

static int igraph_i_barabasi_game_bag(igraph_t *graph, igraph_integer_t n,
                                      igraph_integer_t m,
                                      const igraph_vector_t *outseq,
                                      igraph_bool_t outpref,
                                      igraph_bool_t directed,
                                      const igraph_t *start_from) {

    long int no_of_nodes = n;
    long int no_of_neighbors = m;
    long int *bag;
    long int bagp = 0;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int resp;
    long int i, j, k;
    long int bagsize, start_nodes, start_edges, new_edges, no_of_edges;

    if (!directed) {
        outpref = 1;
    }

    start_nodes = start_from ? igraph_vcount(start_from) : 1;
    start_edges = start_from ? igraph_ecount(start_from) : 0;
    if (outseq) {
        if (igraph_vector_size(outseq) > 1) {
            new_edges = (long int) (igraph_vector_sum(outseq) - VECTOR(*outseq)[0]);
        } else {
            new_edges = 0;
        }
    } else {
        new_edges = (no_of_nodes - start_nodes) * no_of_neighbors;
    }
    no_of_edges = start_edges + new_edges;
    resp = start_edges * 2;
    bagsize = no_of_nodes + no_of_edges + (outpref ? no_of_edges : 0);

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    bag = IGRAPH_CALLOC(bagsize, long int);
    if (bag == 0) {
        IGRAPH_ERROR("barabasi_game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, bag);

    /* The first node(s) in the bag */
    if (start_from) {
        igraph_vector_t deg;
        long int ii, jj, sn = igraph_vcount(start_from);
        igraph_neimode_t mm = outpref ? IGRAPH_ALL : IGRAPH_IN;

        IGRAPH_VECTOR_INIT_FINALLY(&deg, sn);
        IGRAPH_CHECK(igraph_degree(start_from, &deg, igraph_vss_all(), mm,
                                   IGRAPH_LOOPS));
        for (ii = 0; ii < sn; ii++) {
            long int d = (long int) VECTOR(deg)[ii];
            for (jj = 0; jj <= d; jj++) {
                bag[bagp++] = ii;
            }
        }

        igraph_vector_destroy(&deg);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        bag[bagp++] = 0;
    }

    /* Initialize the edges vector */
    if (start_from) {
        IGRAPH_CHECK(igraph_get_edgelist(start_from, &edges, /* bycol= */ 0));
        igraph_vector_resize(&edges, no_of_edges * 2);
    }

    RNG_BEGIN();

    /* and the others */

    for (i = (start_from ? start_nodes : 1), k = (start_from ? 0 : 1);
         i < no_of_nodes; i++, k++) {

        IGRAPH_ALLOW_INTERRUPTION();

        /* draw edges */
        if (outseq) {
            no_of_neighbors = (long int) VECTOR(*outseq)[k];
        }
        for (j = 0; j < no_of_neighbors; j++) {
            long int to = bag[RNG_INTEGER(0, bagp - 1)];
            VECTOR(edges)[resp++] = i;
            VECTOR(edges)[resp++] = to;
        }
        /* update bag */
        bag[bagp++] = i;
        for (j = 0; j < no_of_neighbors; j++) {
            bag[bagp++] = (long int) VECTOR(edges)[resp - 2 * j - 1];
            if (outpref) {
                bag[bagp++] = i;
            }
        }
    }

    RNG_END();

    IGRAPH_FREE(bag);
    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

static int igraph_i_barabasi_game_psumtree_multiple(igraph_t *graph,
                                                    igraph_integer_t n,
                                                    igraph_real_t power,
                                                    igraph_integer_t m,
                                                    const igraph_vector_t *outseq,
                                                    igraph_bool_t outpref,
                                                    igraph_real_t A,
                                                    igraph_bool_t directed,
                                                    const igraph_t *start_from) {

    long int no_of_nodes = n;
    long int no_of_neighbors = m;
    igraph_vector_t edges;
    long int i, j, k;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;
    long int start_nodes, start_edges, new_edges, no_of_edges;

    if (!directed) {
        outpref = 1;
    }

    start_nodes = start_from ? igraph_vcount(start_from) : 1;
    start_edges = start_from ? igraph_ecount(start_from) : 0;
    if (outseq) {
        if (igraph_vector_size(outseq) > 1) {
            new_edges = (long int) (igraph_vector_sum(outseq) - VECTOR(*outseq)[0]);
        } else {
            new_edges = 0;
        }
    } else {
        new_edges = (no_of_nodes - start_nodes) * no_of_neighbors;
    }
    no_of_edges = start_edges + new_edges;
    edgeptr = start_edges * 2;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

    /* first node(s) */
    if (start_from) {
        long int ii, sn = igraph_vcount(start_from);
        igraph_neimode_t mm = outpref ? IGRAPH_ALL : IGRAPH_IN;
        IGRAPH_CHECK(igraph_degree(start_from, &degree, igraph_vss_all(), mm,
                                   IGRAPH_LOOPS));
        IGRAPH_CHECK(igraph_vector_resize(&degree,  no_of_nodes));
        for (ii = 0; ii < sn; ii++) {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, ii, pow(VECTOR(degree)[ii], power) + A));
        }
    } else {
        IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, A));
    }

    /* Initialize the edges vector */
    if (start_from) {
        IGRAPH_CHECK(igraph_get_edgelist(start_from, &edges, /* bycol= */ 0));
        igraph_vector_resize(&edges, no_of_edges * 2);
    }

    RNG_BEGIN();

    /* and the rest */
    for (i = (start_from ? start_nodes : 1), k = (start_from ? 0 : 1);
         i < no_of_nodes; i++, k++) {
        igraph_real_t sum = igraph_psumtree_sum(&sumtree);
        long int to;

        IGRAPH_ALLOW_INTERRUPTION();

        if (outseq) {
            no_of_neighbors = (long int) VECTOR(*outseq)[k];
        }
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
        }
        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            long int nn = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, nn, pow(VECTOR(degree)[nn], power) + A));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, pow(VECTOR(degree)[i], power) + A));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, A));
        }
    }

    RNG_END();

    igraph_psumtree_destroy(&sumtree);
    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_barabasi_game_psumtree(igraph_t *graph,
                                           igraph_integer_t n,
                                           igraph_real_t power,
                                           igraph_integer_t m,
                                           const igraph_vector_t *outseq,
                                           igraph_bool_t outpref,
                                           igraph_real_t A,
                                           igraph_bool_t directed,
                                           const igraph_t *start_from) {

    long int no_of_nodes = n;
    long int no_of_neighbors = m;
    igraph_vector_t edges;
    long int i, j, k;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;
    long int start_nodes, start_edges, new_edges, no_of_edges;

    if (!directed) {
        outpref = 1;
    }

    start_nodes = start_from ? igraph_vcount(start_from) : 1;
    start_edges = start_from ? igraph_ecount(start_from) : 0;
    if (outseq) {
        if (igraph_vector_size(outseq) > 1) {
            new_edges = (long int) (igraph_vector_sum(outseq) - VECTOR(*outseq)[0]);
        } else {
            new_edges = 0;
        }
    } else {
        new_edges = (no_of_nodes - start_nodes) * no_of_neighbors;
    }
    no_of_edges = start_edges + new_edges;
    edgeptr = start_edges * 2;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

    RNG_BEGIN();

    /* first node(s) */
    if (start_from) {
        long int ii, sn = igraph_vcount(start_from);
        igraph_neimode_t mm = outpref ? IGRAPH_ALL : IGRAPH_IN;
        IGRAPH_CHECK(igraph_degree(start_from, &degree, igraph_vss_all(), mm,
                                   IGRAPH_LOOPS));
        IGRAPH_CHECK(igraph_vector_resize(&degree,  no_of_nodes));
        for (ii = 0; ii < sn; ii++) {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, ii, pow(VECTOR(degree)[ii], power) + A));
        }
    } else {
        IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, A));
    }

    /* Initialize the edges vector */
    if (start_from) {
        IGRAPH_CHECK(igraph_get_edgelist(start_from, &edges, /* bycol= */ 0));
    }

    /* and the rest */
    for (i = (start_from ? start_nodes : 1), k = (start_from ? 0 : 1);
         i < no_of_nodes; i++, k++) {
        igraph_real_t sum;
        long int to;

        IGRAPH_ALLOW_INTERRUPTION();

        if (outseq) {
            no_of_neighbors = (long int) VECTOR(*outseq)[k];
        }
        if (no_of_neighbors >= i) {
            /* All existing vertices are cited */
            for (to = 0; to < i; to++) {
                VECTOR(degree)[to]++;
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
                edgeptr += 2;
                IGRAPH_CHECK(igraph_psumtree_update(&sumtree, to, pow(VECTOR(degree)[to], power) + A));
            }
        } else {
            for (j = 0; j < no_of_neighbors; j++) {
                sum = igraph_psumtree_sum(&sumtree);
                igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
                VECTOR(degree)[to]++;
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
                edgeptr += 2;
                IGRAPH_CHECK(igraph_psumtree_update(&sumtree, to, 0.0));
            }
            /* update probabilities */
            for (j = 0; j < no_of_neighbors; j++) {
                long int nn = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
                IGRAPH_CHECK(igraph_psumtree_update(&sumtree, nn, pow(VECTOR(degree)[nn], power) + A));
            }
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors > i ? i : no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, pow(VECTOR(degree)[i], power) + A));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, A));
        }
    }

    RNG_END();

    igraph_psumtree_destroy(&sumtree);
    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \ingroup generators
 * \function igraph_barabasi_game
 * \brief Generates a graph based on the Barab&aacute;si-Albert model.
 *
 * \param graph An uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param power Power of the preferential attachment. The probability
 *        that a vertex is cited is proportional to d^power+A, where
 *        d is its degree (see also the \p outpref argument), power
 *        and A are given by arguments. In the classic preferential
 *        attachment model power=1.
 * \param m The number of outgoing edges generated for each
 *        vertex. (Only if \p outseq is \c NULL.)
 * \param outseq Gives the (out-)degrees of the vertices. If this is
 *        constant, this can be a NULL pointer or an empty (but
 *        initialized!) vector, in this case \p m contains
 *        the constant out-degree. The very first vertex has by definition
 *        no outgoing edges, so the first number in this vector is
 *        ignored.
 * \param outpref Boolean, if true not only the in- but also the out-degree
 *        of a vertex increases its citation probability. I.e., the
 *        citation probability is determined by the total degree of
 *        the vertices. Ignored and assumed to be true if the graph
 *        being generated is undirected.
 * \param A The probability that a vertex is cited is proportional to
 *        d^power+A, where d is its degree (see also the \p outpref
 *        argument), power and A are given by arguments. In the
 *        previous versions of the function this parameter was
 *        implicitly set to one.
 * \param directed Boolean, whether to generate a directed graph.
 * \param algo The algorithm to use to generate the network. Possible
 *        values:
 *        \clist
 *        \cli IGRAPH_BARABASI_BAG
 *          This is the algorithm that was previously (before version
 *          0.6) solely implemented in igraph. It works by putting the
 *          ids of the vertices into a bag (multiset, really), exactly
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
 * \param start_from Either a null pointer, or a graph. In the former
 *        case, the starting configuration is a clique of size \p m.
 *        In the latter case, the graph is a starting configuration.
 *        The graph must be non-empty, i.e. it must have at least one
 *        vertex. If a graph is supplied here and the \p outseq
 *        argument is also given, then \p outseq should only contain
 *        information on the vertices that are not in the \p
 *        start_from graph.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid \p n,
 *         \p m or \p outseq parameter.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges.
 *
 * \example examples/simple/igraph_barabasi_game.c
 * \example examples/simple/igraph_barabasi_game2.c
 */
int igraph_barabasi_game(igraph_t *graph, igraph_integer_t n,
                         igraph_real_t power,
                         igraph_integer_t m,
                         const igraph_vector_t *outseq,
                         igraph_bool_t outpref,
                         igraph_real_t A,
                         igraph_bool_t directed,
                         igraph_barabasi_algorithm_t algo,
                         const igraph_t *start_from) {

    long int start_nodes = start_from ? igraph_vcount(start_from) : 0;
    long int newn = start_from ? n - start_nodes : n;

    /* Fix obscure parameterizations */
    if (outseq && igraph_vector_size(outseq) == 0) {
        outseq = 0;
    }
    if (!directed) {
        outpref = 1;
    }

    /* Check arguments */

    if (algo != IGRAPH_BARABASI_BAG &&
        algo != IGRAPH_BARABASI_PSUMTREE &&
        algo != IGRAPH_BARABASI_PSUMTREE_MULTIPLE) {
        IGRAPH_ERROR("Invalid algorithm", IGRAPH_EINVAL);
    }
    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    } else if (newn < 0) {
        IGRAPH_ERROR("Starting graph has too many vertices.", IGRAPH_EINVAL);
    }
    if (start_from && start_nodes == 0) {
        IGRAPH_ERROR("Cannot start from an empty graph.", IGRAPH_EINVAL);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 &&
        igraph_vector_size(outseq) != newn) {
        IGRAPH_ERROR("Invalid out-degree sequence length.", IGRAPH_EINVAL);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
        IGRAPH_ERROR("Number of edges added per step must not be negative.", IGRAPH_EINVAL);
    }
    if (outseq && igraph_vector_min(outseq) < 0) {
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
            IGRAPH_ERROR("Power must be one for 'bag' algorithm.", IGRAPH_EINVAL);
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

    if (algo == IGRAPH_BARABASI_BAG) {
        return igraph_i_barabasi_game_bag(graph, n, m, outseq, outpref, directed,
                                          start_from);
    } else if (algo == IGRAPH_BARABASI_PSUMTREE) {
        return igraph_i_barabasi_game_psumtree(graph, n, power, m, outseq,
                                               outpref, A, directed, start_from);
    } else if (algo == IGRAPH_BARABASI_PSUMTREE_MULTIPLE) {
        return igraph_i_barabasi_game_psumtree_multiple(graph, n, power, m,
                outseq, outpref, A,
                directed, start_from);
    }

    return 0;
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
 * which are summed to obtain the final weight.
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
int igraph_barabasi_aging_game(igraph_t *graph,
                               igraph_integer_t nodes,
                               igraph_integer_t m,
                               const igraph_vector_t *outseq,
                               igraph_bool_t outpref,
                               igraph_real_t pa_exp,
                               igraph_real_t aging_exp,
                               igraph_integer_t aging_bins,
                               igraph_real_t zero_deg_appeal,
                               igraph_real_t zero_age_appeal,
                               igraph_real_t deg_coef,
                               igraph_real_t age_coef,
                               igraph_bool_t directed) {
    long int no_of_nodes = nodes;
    long int no_of_neighbors = m;
    long int binwidth;
    long int no_of_edges;
    igraph_vector_t edges;
    long int i, j, k;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;

    if (no_of_nodes < 0) {
        IGRAPH_ERRORF("Number of nodes must not be negative, got %ld.", IGRAPH_EINVAL, no_of_nodes);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != no_of_nodes) {
        IGRAPH_ERRORF("The length of the out-degree sequence (%ld) does not agree with the number of nodes (%ld).",
                       IGRAPH_EINVAL,
                       igraph_vector_size(outseq), no_of_nodes);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
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

    if (outseq == 0 || igraph_vector_size(outseq) == 0) {
        no_of_neighbors = m;
        no_of_edges = (no_of_nodes - 1) * no_of_neighbors;
    } else {
        no_of_edges = 0;
        for (i = 1; i < igraph_vector_size(outseq); i++) {
            no_of_edges += VECTOR(*outseq)[i];
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

    RNG_BEGIN();

    /* first node */
    IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, zero_deg_appeal * (1 + zero_age_appeal)));

    /* and the rest */
    for (i = 1; i < no_of_nodes; i++) {
        igraph_real_t sum;
        long int to;

        IGRAPH_ALLOW_INTERRUPTION();

        if (outseq != 0 && igraph_vector_size(outseq) != 0) {
            no_of_neighbors = (long int) VECTOR(*outseq)[i];
        }
        sum = igraph_psumtree_sum(&sumtree);
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
        }
        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            long int n = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
            long int age = (i - n) / binwidth;
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, n,
                (deg_coef * pow(VECTOR(degree)[n], pa_exp) + zero_deg_appeal) *
                (age_coef * pow(age + 1, aging_exp) + zero_age_appeal)
            ));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, i,
                (zero_age_appeal + 1) * (deg_coef * pow(VECTOR(degree)[i], pa_exp) + zero_deg_appeal)
            ));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, i, (1 + zero_age_appeal) * zero_deg_appeal
            ));
        }

        /* aging */
        for (k = 1; i - binwidth * k + 1 >= 1; k++) {
            long int shnode = i - binwidth * k;
            long int deg = (long int) VECTOR(degree)[shnode];
            long int age = (i - shnode) / binwidth;
            /* igraph_real_t old=igraph_psumtree_get(&sumtree, shnode); */
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, shnode,
                (deg_coef * pow(deg, pa_exp) + zero_deg_appeal) *
                (age_coef * pow(age + 2, aging_exp) + zero_age_appeal)
                ));
        }
    }

    RNG_END();

    igraph_vector_destroy(&degree);
    igraph_psumtree_destroy(&sumtree);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
