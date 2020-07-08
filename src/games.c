/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph R library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_interface.h"
#include "igraph_games.h"
#include "igraph_random.h"
#include "igraph_memory.h"
#include "igraph_interrupt_internal.h"
#include "igraph_attributes.h"
#include "igraph_constructors.h"
#include "igraph_nongraph.h"
#include "igraph_conversion.h"
#include "igraph_psumtree.h"
#include "igraph_dqueue.h"
#include "igraph_adjlist.h"
#include "igraph_iterators.h"
#include "igraph_progress.h"
#include "igraph_topology.h"
#include "igraph_types_internal.h"
#include "config.h"

#include <math.h>

typedef struct {
    long int no;
    igraph_psumtree_t *sumtrees;
} igraph_i_citing_cited_type_game_struct_t;

static void igraph_i_citing_cited_type_game_free (
        igraph_i_citing_cited_type_game_struct_t *s);
/**
 * \section about_games
 *
 * <para>Games are randomized graph generators. Randomization means that
 * they generate a different graph every time you call them. </para>
 */

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

    bag = igraph_Calloc(bagsize, long int);
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

    igraph_Free(bag);
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
            igraph_psumtree_update(&sumtree, ii, pow(VECTOR(degree)[ii], power) + A);
        }
    } else {
        igraph_psumtree_update(&sumtree, 0, A);
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
            igraph_psumtree_update(&sumtree, nn,
                                   pow(VECTOR(degree)[nn], power) + A);
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            igraph_psumtree_update(&sumtree, i,
                                   pow(VECTOR(degree)[i], power) + A);
        } else {
            igraph_psumtree_update(&sumtree, i, A);
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
            igraph_psumtree_update(&sumtree, ii, pow(VECTOR(degree)[ii], power) + A);
        }
    } else {
        igraph_psumtree_update(&sumtree, 0, A);
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
        if (outseq) {
            no_of_neighbors = (long int) VECTOR(*outseq)[k];
        }
        if (no_of_neighbors >= i) {
            /* All existing vertices are cited */
            for (to = 0; to < i; to++) {
                VECTOR(degree)[to]++;
                igraph_vector_push_back(&edges, i);
                igraph_vector_push_back(&edges, to);
                edgeptr += 2;
                igraph_psumtree_update(&sumtree, to, pow(VECTOR(degree)[to], power) + A);
            }
        } else {
            for (j = 0; j < no_of_neighbors; j++) {
                sum = igraph_psumtree_sum(&sumtree);
                igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
                VECTOR(degree)[to]++;
                igraph_vector_push_back(&edges, i);
                igraph_vector_push_back(&edges, to);
                edgeptr += 2;
                igraph_psumtree_update(&sumtree, to, 0.0);
            }
            /* update probabilities */
            for (j = 0; j < no_of_neighbors; j++) {
                long int nn = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
                igraph_psumtree_update(&sumtree, nn,
                                       pow(VECTOR(degree)[nn], power) + A);
            }
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors > i ? i : no_of_neighbors;
            igraph_psumtree_update(&sumtree, i,
                                   pow(VECTOR(degree)[i], power) + A);
        } else {
            igraph_psumtree_update(&sumtree, i, A);
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
 *        of a vertex increases its citation probability. Ie. the
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
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    } else if (newn < 0) {
        IGRAPH_ERROR("Starting graph has too many vertices", IGRAPH_EINVAL);
    }
    if (start_from && start_nodes == 0) {
        IGRAPH_ERROR("Cannot start from an empty graph", IGRAPH_EINVAL);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 &&
        igraph_vector_size(outseq) != newn) {
        IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
        IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
    }
    if (outseq && igraph_vector_min(outseq) < 0) {
        IGRAPH_ERROR("Negative out degree in sequence", IGRAPH_EINVAL);
    }
    if (!outpref && A <= 0) {
        IGRAPH_ERROR("Constant attractiveness (A) must be positive",
                     IGRAPH_EINVAL);
    }
    if (outpref && A < 0) {
        IGRAPH_ERROR("Constant attractiveness (A) must be non-negative",
                     IGRAPH_EINVAL);
    }
    if (algo == IGRAPH_BARABASI_BAG) {
        if (power != 1) {
            IGRAPH_ERROR("Power must be one for 'bag' algorithm", IGRAPH_EINVAL);
        }
        if (A != 1) {
            IGRAPH_ERROR("Constant attractiveness (A) must be one for bag algorithm",
                         IGRAPH_EINVAL);
        }
    }
    if (start_from && directed != igraph_is_directed(start_from)) {
        IGRAPH_WARNING("Directedness of the start graph and the output graph"
                       " mismatch");
    }
    if (start_from && !igraph_is_directed(start_from) && !outpref) {
        IGRAPH_ERROR("`outpref' must be true if starting from an undirected "
                     "graph", IGRAPH_EINVAL);
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
 * \ingroup internal
 */

int igraph_erdos_renyi_game_gnp(igraph_t *graph, igraph_integer_t n, igraph_real_t p,
                                igraph_bool_t directed, igraph_bool_t loops) {

    long int no_of_nodes = n;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int retval = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (p < 0.0 || p > 1.0) {
        IGRAPH_ERROR("Invalid probability given", IGRAPH_EINVAL);
    }

    if (p == 0.0 || no_of_nodes <= 1) {
        IGRAPH_CHECK(retval = igraph_empty(graph, n, directed));
    } else if (p == 1.0) {
        IGRAPH_CHECK(retval = igraph_full(graph, n, directed, loops));
    } else {

        long int i;
        double maxedges = n, last;
        if (directed && loops) {
            maxedges *= n;
        } else if (directed && !loops) {
            maxedges *= (n - 1);
        } else if (!directed && loops) {
            maxedges *= (n + 1) / 2.0;
        } else {
            maxedges *= (n - 1) / 2.0;
        }

        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&s, (long int) (maxedges * p * 1.1)));

        RNG_BEGIN();

        last = RNG_GEOM(p);
        while (last < maxedges) {
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            last += RNG_GEOM(p);
            last += 1;
        }

        RNG_END();

        IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s) * 2));

        if (directed && loops) {
            for (i = 0; i < igraph_vector_size(&s); i++) {
                long int to = (long int) floor(VECTOR(s)[i] / no_of_nodes);
                long int from = (long int) (VECTOR(s)[i] - ((igraph_real_t)to) * no_of_nodes);
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        } else if (directed && !loops) {
            for (i = 0; i < igraph_vector_size(&s); i++) {
                long int to = (long int) floor(VECTOR(s)[i] / no_of_nodes);
                long int from = (long int) (VECTOR(s)[i] - ((igraph_real_t)to) * no_of_nodes);
                if (from == to) {
                    to = no_of_nodes - 1;
                }
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        } else if (!directed && loops) {
            for (i = 0; i < igraph_vector_size(&s); i++) {
                long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2);
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        } else { /* !directed && !loops */
            for (i = 0; i < igraph_vector_size(&s); i++) {
                long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2);
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        }

        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_CHECK(retval = igraph_create(graph, &edges, n, directed));
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return retval;
}

int igraph_erdos_renyi_game_gnm(igraph_t *graph, igraph_integer_t n, igraph_real_t m,
                                igraph_bool_t directed, igraph_bool_t loops) {

    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_edges = (igraph_integer_t) m;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int retval = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (m < 0) {
        IGRAPH_ERROR("Invalid number of edges", IGRAPH_EINVAL);
    }

    if (m == 0.0 || no_of_nodes <= 1) {
        IGRAPH_CHECK(retval = igraph_empty(graph, n, directed));
    } else {

        long int i;
        double maxedges = n;
        if (directed && loops) {
            maxedges *= n;
        } else if (directed && !loops) {
            maxedges *= (n - 1);
        } else if (!directed && loops) {
            maxedges *= (n + 1) / 2.0;
        } else {
            maxedges *= (n - 1) / 2.0;
        }

        if (no_of_edges > maxedges) {
            IGRAPH_ERROR("Invalid number (too large) of edges", IGRAPH_EINVAL);
        }

        if (maxedges == no_of_edges) {
            retval = igraph_full(graph, n, directed, loops);
        } else {

            long int slen;

            IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
            IGRAPH_CHECK(igraph_random_sample(&s, 0, maxedges - 1,
                                              (igraph_integer_t) no_of_edges));

            IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
            IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s) * 2));

            slen = igraph_vector_size(&s);
            if (directed && loops) {
                for (i = 0; i < slen; i++) {
                    long int to = (long int) floor(VECTOR(s)[i] / no_of_nodes);
                    long int from = (long int) (VECTOR(s)[i] - ((igraph_real_t)to) * no_of_nodes);
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            } else if (directed && !loops) {
                for (i = 0; i < slen; i++) {
                    long int from = (long int) floor(VECTOR(s)[i] / (no_of_nodes - 1));
                    long int to = (long int) (VECTOR(s)[i] - ((igraph_real_t)from) * (no_of_nodes - 1));
                    if (from == to) {
                        to = no_of_nodes - 1;
                    }
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            } else if (!directed && loops) {
                for (i = 0; i < slen; i++) {
                    long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                    long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2);
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            } else { /* !directed && !loops */
                for (i = 0; i < slen; i++) {
                    long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                    long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2);
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            }

            igraph_vector_destroy(&s);
            IGRAPH_FINALLY_CLEAN(1);
            retval = igraph_create(graph, &edges, n, directed);
            igraph_vector_destroy(&edges);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    return retval;
}

/**
 * \ingroup generators
 * \function igraph_erdos_renyi_game
 * \brief Generates a random (Erdos-Renyi) graph.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param type The type of the random graph, possible values:
 *        \clist
 *        \cli IGRAPH_ERDOS_RENYI_GNM
 *          G(n,m) graph,
 *          m edges are
 *          selected uniformly randomly in a graph with
 *          n vertices.
 *        \cli IGRAPH_ERDOS_RENYI_GNP
 *          G(n,p) graph,
 *          every possible edge is included in the graph with
 *          probability p.
 *        \endclist
 * \param n The number of vertices in the graph.
 * \param p_or_m This is the p parameter for
 *        G(n,p) graphs and the
 *        m
 *        parameter for G(n,m) graphs.
 * \param directed Logical, whether to generate a directed graph.
 * \param loops Logical, whether to generate loops (self) edges.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid
 *         \p type, \p n,
 *         \p p or \p m
 *          parameter.
 *         \c IGRAPH_ENOMEM: there is not enough
 *         memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_barabasi_game(), \ref igraph_growing_random_game()
 *
 * \example examples/simple/igraph_erdos_renyi_game.c
 */

int igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
                            igraph_integer_t n, igraph_real_t p_or_m,
                            igraph_bool_t directed, igraph_bool_t loops) {
    int retval = 0;
    if (type == IGRAPH_ERDOS_RENYI_GNP) {
        retval = igraph_erdos_renyi_game_gnp(graph, n, p_or_m, directed, loops);
    } else if (type == IGRAPH_ERDOS_RENYI_GNM) {
        retval = igraph_erdos_renyi_game_gnm(graph, n, p_or_m, directed, loops);
    } else {
        IGRAPH_ERROR("Invalid type", IGRAPH_EINVAL);
    }

    return retval;
}

int igraph_degree_sequence_game_simple(igraph_t *graph,
                                       const igraph_vector_t *out_seq,
                                       const igraph_vector_t *in_seq);

int igraph_degree_sequence_game_simple(igraph_t *graph,
                                       const igraph_vector_t *out_seq,
                                       const igraph_vector_t *in_seq) {

    long int outsum = 0, insum = 0;
    igraph_bool_t directed = (in_seq != 0 && igraph_vector_size(in_seq) != 0);
    igraph_bool_t degseq_ok;
    long int no_of_nodes, no_of_edges;
    long int *bag1 = 0, *bag2 = 0;
    long int bagp1 = 0, bagp2 = 0;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int i, j;

    IGRAPH_CHECK(igraph_is_degree_sequence(out_seq, in_seq, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR(in_seq ? "No directed graph can realize the given degree sequences" :
                     "No undirected graph can realize the given degree sequence", IGRAPH_EINVAL);
    }

    outsum = (long int) igraph_vector_sum(out_seq);
    if (directed) {
        insum = (long int) igraph_vector_sum(in_seq);
    }

    no_of_nodes = igraph_vector_size(out_seq);
    no_of_edges = directed ? outsum : outsum / 2;

    bag1 = igraph_Calloc(outsum, long int);
    if (bag1 == 0) {
        IGRAPH_ERROR("degree sequence game (simple)", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, bag1);

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < VECTOR(*out_seq)[i]; j++) {
            bag1[bagp1++] = i;
        }
    }
    if (directed) {
        bag2 = igraph_Calloc(insum, long int);
        if (bag2 == 0) {
            IGRAPH_ERROR("degree sequence game (simple)", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, bag2);
        for (i = 0; i < no_of_nodes; i++) {
            for (j = 0; j < VECTOR(*in_seq)[i]; j++) {
                bag2[bagp2++] = i;
            }
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    RNG_BEGIN();

    if (directed) {
        for (i = 0; i < no_of_edges; i++) {
            long int from = RNG_INTEGER(0, bagp1 - 1);
            long int to = RNG_INTEGER(0, bagp2 - 1);
            igraph_vector_push_back(&edges, bag1[from]); /* safe, already reserved */
            igraph_vector_push_back(&edges, bag2[to]);   /* ditto */
            bag1[from] = bag1[bagp1 - 1];
            bag2[to] = bag2[bagp2 - 1];
            bagp1--; bagp2--;
        }
    } else {
        for (i = 0; i < no_of_edges; i++) {
            long int from = RNG_INTEGER(0, bagp1 - 1);
            long int to;
            igraph_vector_push_back(&edges, bag1[from]); /* safe, already reserved */
            bag1[from] = bag1[bagp1 - 1];
            bagp1--;
            to = RNG_INTEGER(0, bagp1 - 1);
            igraph_vector_push_back(&edges, bag1[to]);   /* ditto */
            bag1[to] = bag1[bagp1 - 1];
            bagp1--;
        }
    }

    RNG_END();

    igraph_Free(bag1);
    IGRAPH_FINALLY_CLEAN(1);
    if (directed) {
        igraph_Free(bag2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

int igraph_degree_sequence_game_no_multiple_undirected(
    igraph_t *graph, const igraph_vector_t *seq) {

    igraph_vector_t stubs = IGRAPH_VECTOR_NULL;
    igraph_vector_int_t *neis;
    igraph_vector_t residual_degrees = IGRAPH_VECTOR_NULL;
    igraph_set_t incomplete_vertices;
    igraph_adjlist_t al;
    igraph_bool_t finished, failed;
    igraph_integer_t from, to, dummy;
    long int i, j, k;
    long int no_of_nodes, outsum = 0;
    igraph_bool_t degseq_ok;

    IGRAPH_CHECK(igraph_is_graphical_degree_sequence(seq, 0, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR("No simple undirected graph can realize the given degree sequence",
                     IGRAPH_EINVAL);
    }

    outsum = (long int) igraph_vector_sum(seq);
    no_of_nodes = igraph_vector_size(seq);

    /* Allocate required data structures */
    IGRAPH_CHECK(igraph_adjlist_init_empty(&al, (igraph_integer_t) no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    IGRAPH_VECTOR_INIT_FINALLY(&stubs, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&stubs, outsum));
    IGRAPH_VECTOR_INIT_FINALLY(&residual_degrees, no_of_nodes);
    IGRAPH_CHECK(igraph_set_init(&incomplete_vertices, 0));
    IGRAPH_FINALLY(igraph_set_destroy, &incomplete_vertices);

    /* Start the RNG */
    RNG_BEGIN();

    /* Outer loop; this will try to construct a graph several times from scratch
     * until it finally succeeds. */
    finished = 0;
    while (!finished) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* Be optimistic :) */
        failed = 0;

        /* Clear the adjacency list to get rid of the previous attempt (if any) */
        igraph_adjlist_clear(&al);

        /* Initialize the residual degrees from the degree sequence */
        IGRAPH_CHECK(igraph_vector_update(&residual_degrees, seq));

        /* While there are some unconnected stubs left... */
        while (!finished && !failed) {
            /* Construct the initial stub vector */
            igraph_vector_clear(&stubs);
            for (i = 0; i < no_of_nodes; i++) {
                for (j = 0; j < VECTOR(residual_degrees)[i]; j++) {
                    igraph_vector_push_back(&stubs, i);
                }
            }

            /* Clear the skipped stub counters and the set of incomplete vertices */
            igraph_vector_null(&residual_degrees);
            igraph_set_clear(&incomplete_vertices);

            /* Shuffle the stubs in-place */
            igraph_vector_shuffle(&stubs);

            /* Connect the stubs where possible */
            k = igraph_vector_size(&stubs);
            for (i = 0; i < k; ) {
                from = (igraph_integer_t) VECTOR(stubs)[i++];
                to = (igraph_integer_t) VECTOR(stubs)[i++];

                if (from > to) {
                    dummy = from; from = to; to = dummy;
                }

                neis = igraph_adjlist_get(&al, from);
                if (from == to || igraph_vector_int_binsearch(neis, to, &j)) {
                    /* Edge exists already */
                    VECTOR(residual_degrees)[from]++;
                    VECTOR(residual_degrees)[to]++;
                    IGRAPH_CHECK(igraph_set_add(&incomplete_vertices, from));
                    IGRAPH_CHECK(igraph_set_add(&incomplete_vertices, to));
                } else {
                    /* Insert the edge */
                    IGRAPH_CHECK(igraph_vector_int_insert(neis, j, to));
                }
            }

            finished = igraph_set_empty(&incomplete_vertices);

            if (!finished) {
                /* We are not done yet; check if the remaining stubs are feasible. This
                 * is done by enumerating all possible pairs and checking whether at
                 * least one feasible pair is found. */
                i = 0;
                failed = 1;
                while (failed && igraph_set_iterate(&incomplete_vertices, &i, &from)) {
                    j = 0;
                    while (igraph_set_iterate(&incomplete_vertices, &j, &to)) {
                        if (from == to) {
                            /* This is used to ensure that each pair is checked once only */
                            break;
                        }
                        if (from > to) {
                            dummy = from; from = to; to = dummy;
                        }
                        neis = igraph_adjlist_get(&al, from);
                        if (!igraph_vector_int_binsearch(neis, to, 0)) {
                            /* Found a suitable pair, so we can continue */
                            failed = 0;
                            break;
                        }
                    }
                }
            }
        }
    }

    /* Finish the RNG */
    RNG_END();

    /* Clean up */
    igraph_set_destroy(&incomplete_vertices);
    igraph_vector_destroy(&residual_degrees);
    igraph_vector_destroy(&stubs);
    IGRAPH_FINALLY_CLEAN(3);

    /* Create the graph. We cannot use IGRAPH_ALL here for undirected graphs
     * because we did not add edges in both directions in the adjacency list.
     * We will use igraph_to_undirected in an extra step. */
    IGRAPH_CHECK(igraph_adjlist(graph, &al, IGRAPH_OUT, 1));
    IGRAPH_CHECK(igraph_to_undirected(graph, IGRAPH_TO_UNDIRECTED_EACH, 0));

    /* Clear the adjacency list */
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

int igraph_degree_sequence_game_no_multiple_directed(igraph_t *graph,
        const igraph_vector_t *out_seq, const igraph_vector_t *in_seq) {
    igraph_adjlist_t al;
    igraph_bool_t deg_seq_ok, failed, finished;
    igraph_vector_t in_stubs = IGRAPH_VECTOR_NULL;
    igraph_vector_t out_stubs = IGRAPH_VECTOR_NULL;
    igraph_vector_int_t *neis;
    igraph_vector_t residual_in_degrees = IGRAPH_VECTOR_NULL;
    igraph_vector_t residual_out_degrees = IGRAPH_VECTOR_NULL;
    igraph_set_t incomplete_in_vertices;
    igraph_set_t incomplete_out_vertices;
    igraph_integer_t from, to;
    long int i, j, k;
    long int no_of_nodes, outsum;

    IGRAPH_CHECK(igraph_is_graphical_degree_sequence(out_seq, in_seq, &deg_seq_ok));
    if (!deg_seq_ok) {
        IGRAPH_ERROR("No simple directed graph can realize the given degree sequence",
                     IGRAPH_EINVAL);
    }

    outsum = (long int) igraph_vector_sum(out_seq);
    no_of_nodes = igraph_vector_size(out_seq);

    /* Allocate required data structures */
    IGRAPH_CHECK(igraph_adjlist_init_empty(&al, (igraph_integer_t) no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    IGRAPH_VECTOR_INIT_FINALLY(&out_stubs, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&out_stubs, outsum));
    IGRAPH_VECTOR_INIT_FINALLY(&in_stubs, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&in_stubs, outsum));
    IGRAPH_VECTOR_INIT_FINALLY(&residual_out_degrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&residual_in_degrees, no_of_nodes);
    IGRAPH_CHECK(igraph_set_init(&incomplete_out_vertices, 0));
    IGRAPH_FINALLY(igraph_set_destroy, &incomplete_out_vertices);
    IGRAPH_CHECK(igraph_set_init(&incomplete_in_vertices, 0));
    IGRAPH_FINALLY(igraph_set_destroy, &incomplete_in_vertices);

    /* Start the RNG */
    RNG_BEGIN();

    /* Outer loop; this will try to construct a graph several times from scratch
     * until it finally succeeds. */
    finished = 0;
    while (!finished) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* Be optimistic :) */
        failed = 0;

        /* Clear the adjacency list to get rid of the previous attempt (if any) */
        igraph_adjlist_clear(&al);

        /* Initialize the residual degrees from the degree sequences */
        IGRAPH_CHECK(igraph_vector_update(&residual_out_degrees, out_seq));
        IGRAPH_CHECK(igraph_vector_update(&residual_in_degrees, in_seq));

        /* While there are some unconnected stubs left... */
        while (!finished && !failed) {
            /* Construct the initial stub vectors */
            igraph_vector_clear(&out_stubs);
            igraph_vector_clear(&in_stubs);
            for (i = 0; i < no_of_nodes; i++) {
                for (j = 0; j < VECTOR(residual_out_degrees)[i]; j++) {
                    igraph_vector_push_back(&out_stubs, i);
                }
                for (j = 0; j < VECTOR(residual_in_degrees)[i]; j++) {
                    igraph_vector_push_back(&in_stubs, i);
                }
            }

            /* Clear the skipped stub counters and the set of incomplete vertices */
            igraph_vector_null(&residual_out_degrees);
            igraph_vector_null(&residual_in_degrees);
            igraph_set_clear(&incomplete_out_vertices);
            igraph_set_clear(&incomplete_in_vertices);
            outsum = 0;

            /* Shuffle the out-stubs in-place */
            igraph_vector_shuffle(&out_stubs);

            /* Connect the stubs where possible */
            k = igraph_vector_size(&out_stubs);
            for (i = 0; i < k; i++) {
                from = (igraph_integer_t) VECTOR(out_stubs)[i];
                to = (igraph_integer_t) VECTOR(in_stubs)[i];

                neis = igraph_adjlist_get(&al, from);
                if (from == to || igraph_vector_int_binsearch(neis, to, &j)) {
                    /* Edge exists already */
                    VECTOR(residual_out_degrees)[from]++;
                    VECTOR(residual_in_degrees)[to]++;
                    IGRAPH_CHECK(igraph_set_add(&incomplete_out_vertices, from));
                    IGRAPH_CHECK(igraph_set_add(&incomplete_in_vertices, to));
                } else {
                    /* Insert the edge */
                    IGRAPH_CHECK(igraph_vector_int_insert(neis, j, to));
                }
            }

            /* Are we finished? */
            finished = igraph_set_empty(&incomplete_out_vertices);

            if (!finished) {
                /* We are not done yet; check if the remaining stubs are feasible. This
                 * is done by enumerating all possible pairs and checking whether at
                 * least one feasible pair is found. */
                i = 0;
                failed = 1;
                while (failed && igraph_set_iterate(&incomplete_out_vertices, &i, &from)) {
                    j = 0;
                    while (igraph_set_iterate(&incomplete_in_vertices, &j, &to)) {
                        neis = igraph_adjlist_get(&al, from);
                        if (from != to && !igraph_vector_int_binsearch(neis, to, 0)) {
                            /* Found a suitable pair, so we can continue */
                            failed = 0;
                            break;
                        }
                    }
                }
            }
        }
    }

    /* Finish the RNG */
    RNG_END();

    /* Clean up */
    igraph_set_destroy(&incomplete_in_vertices);
    igraph_set_destroy(&incomplete_out_vertices);
    igraph_vector_destroy(&residual_in_degrees);
    igraph_vector_destroy(&residual_out_degrees);
    igraph_vector_destroy(&in_stubs);
    igraph_vector_destroy(&out_stubs);
    IGRAPH_FINALLY_CLEAN(6);

    /* Create the graph */
    IGRAPH_CHECK(igraph_adjlist(graph, &al, IGRAPH_OUT, 1));

    /* Clear the adjacency list */
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* swap two elements of a vector_int */
#define SWAP_INT_ELEM(vec, i, j) \
    { \
        igraph_integer_t temp; \
        temp = VECTOR(vec)[i]; \
        VECTOR(vec)[i] = VECTOR(vec)[j]; \
        VECTOR(vec)[j] = temp; \
    }

int igraph_degree_sequence_game_no_multiple_undirected_uniform(igraph_t *graph, const igraph_vector_t *degseq) {
    igraph_vector_int_t stubs;
    igraph_vector_t edges;
    igraph_bool_t degseq_ok;
    igraph_vector_ptr_t adjlist;
    long i, j;
    long vcount, ecount, stub_count;

    IGRAPH_CHECK(igraph_is_graphical_degree_sequence(degseq, NULL, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR("No simple undirected graph can realize the given degree sequence", IGRAPH_EINVAL);
    }

    stub_count = (long) igraph_vector_sum(degseq);
    ecount = stub_count / 2;
    vcount = igraph_vector_size(degseq);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&stubs, stub_count);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, stub_count);

    /* Fill stubs vector. */
    {
        long k = 0;
        for (i = 0; i < vcount; ++i) {
            long deg = (long) VECTOR(*degseq)[i];
            for (j = 0; j < deg; ++j) {
                VECTOR(stubs)[k++] = i;
            }
        }
    }

    /* Build an adjacency list in terms of sets; used to check for multi-edges. */
    IGRAPH_CHECK(igraph_vector_ptr_init(&adjlist, vcount));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&adjlist, igraph_set_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &adjlist);
    for (i = 0; i < vcount; ++i) {
        igraph_set_t *set = igraph_malloc(sizeof(igraph_set_t));
        if (! set) {
            IGRAPH_ERROR("Out of memory", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_set_init(set, 0));
        VECTOR(adjlist)[i] = set;
        IGRAPH_CHECK(igraph_set_reserve(set, (long) VECTOR(*degseq)[i]));
    }    

    RNG_BEGIN();

    for (;;) {
        igraph_bool_t success = 1;

        /* Shuffle stubs vector with Fisher-Yates and check for self-loops and multi-edges as we go. */
        for (i = 0; i < ecount; ++i) {
            long k, from, to;

            k = RNG_INTEGER(2*i, stub_count-1);
            SWAP_INT_ELEM(stubs, 2*i, k);

            k = RNG_INTEGER(2*i+1, stub_count-1);
            SWAP_INT_ELEM(stubs, 2*i+1, k);

            from = VECTOR(stubs)[2*i];
            to   = VECTOR(stubs)[2*i+1];

            /* self-loop, fail */
            if (from == to) {
                success = 0;
                break;
            }

            /* multi-edge, fail */
            if (igraph_set_contains((igraph_set_t *) VECTOR(adjlist)[to], from)) {
                success = 0;
                break;
            }

            /* sets are already reserved */
            igraph_set_add((igraph_set_t *) VECTOR(adjlist)[to], from);
            igraph_set_add((igraph_set_t *) VECTOR(adjlist)[from], to);

            /* register edge */
            VECTOR(edges)[2 * i]   = from;
            VECTOR(edges)[2 * i + 1] = to;
        }

        if (success) {
            break;
        }

        /* Clear adjacency list. */
        for (j = 0; j < vcount; ++j) {
            igraph_set_clear((igraph_set_t *) VECTOR(adjlist)[j]);
        }

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    igraph_vector_ptr_destroy_all(&adjlist);
    igraph_vector_int_destroy(&stubs);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, vcount, /* directed = */ 0));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


int igraph_degree_sequence_game_no_multiple_directed_uniform(
    igraph_t *graph, const igraph_vector_t *out_deg, const igraph_vector_t *in_deg) {
    igraph_vector_int_t out_stubs, in_stubs;
    igraph_vector_t edges;
    igraph_bool_t degseq_ok;
    igraph_vector_ptr_t adjlist;
    long i, j;
    long vcount, ecount;

    IGRAPH_CHECK(igraph_is_graphical_degree_sequence(out_deg, in_deg, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR("No simple directed graph can realize the given degree sequence", IGRAPH_EINVAL);
    }

    ecount = (long) igraph_vector_sum(out_deg);
    vcount = igraph_vector_size(out_deg);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&out_stubs, ecount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&in_stubs, ecount);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * ecount);

    /* Fill in- and out-stubs vectors. */
    {
        long k = 0, l = 0;
        for (i = 0; i < vcount; ++i) {
            long dout, din;

            dout = (long) VECTOR(*out_deg)[i];
            for (j = 0; j < dout; ++j) {
                VECTOR(out_stubs)[k++] = i;
            }

            din  = (long) VECTOR(*in_deg)[i];
            for (j = 0; j < din; ++j) {
                VECTOR(in_stubs)[l++] = i;
            }
        }
    }

    /* Build an adjacency list in terms of sets; used to check for multi-edges. */
    IGRAPH_CHECK(igraph_vector_ptr_init(&adjlist, vcount));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&adjlist, igraph_set_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &adjlist);
    for (i = 0; i < vcount; ++i) {
        igraph_set_t *set = igraph_malloc(sizeof(igraph_set_t));
        if (! set) {
            IGRAPH_ERROR("Out of memory", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_set_init(set, 0));
        VECTOR(adjlist)[i] = set;
        IGRAPH_CHECK(igraph_set_reserve(set, (long) VECTOR(*out_deg)[i]));
    }    

    RNG_BEGIN();

    for (;;) {
        igraph_bool_t success = 1;

        /* Shuffle out-stubs vector with Fisher-Yates and check for self-loops and multi-edges as we go. */
        for (i = 0; i < ecount; ++i) {
            long k, from, to;
            igraph_set_t *set;

            k = RNG_INTEGER(i, ecount-1);
            SWAP_INT_ELEM(out_stubs, i, k);

            from = VECTOR(out_stubs)[i];
            to   = VECTOR(in_stubs)[i];

            /* self-loop, fail */
            if (to == from) {
                success = 0;
                break;
            }

            /* multi-edge, fail */
            set = (igraph_set_t *) VECTOR(adjlist)[from];
            if (igraph_set_contains(set, to)) {
                success = 0;
                break;
            }

            /* sets are already reserved */
            igraph_set_add(set, to);

            /* register edge */
            VECTOR(edges)[2 * i]   = from;
            VECTOR(edges)[2 * i + 1] = to;
        }

        if (success) {
            break;
        }

        /* Clear adjacency list. */
        for (j = 0; j < vcount; ++j) {
            igraph_set_clear((igraph_set_t *) VECTOR(adjlist)[j]);
        }

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    igraph_vector_ptr_destroy_all(&adjlist);
    igraph_vector_int_destroy(&out_stubs);
    igraph_vector_int_destroy(&in_stubs);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, vcount, /* directed = */ 1));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

#undef SWAP_INT_ELEM


/* This is in gengraph_mr-connected.cpp */

int igraph_degree_sequence_game_vl(igraph_t *graph,
                                   const igraph_vector_t *out_seq,
                                   const igraph_vector_t *in_seq);
/**
 * \ingroup generators
 * \function igraph_degree_sequence_game
 * \brief Generates a random graph with a given degree sequence
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param out_deg The degree sequence for an undirected graph (if
 *        \p in_seq is \c NULL or of length zero), or the out-degree
 *        sequence of a directed graph (if \p in_deq is not
 *        of length zero).
 * \param in_deg It is either a zero-length vector or
 *        \c NULL (if an undirected
 *        graph is generated), or the in-degree sequence.
 * \param method The method to generate the graph. Possible values:
 *        \clist
 *          \cli IGRAPH_DEGSEQ_SIMPLE
 *          This method implements the configuration model.
 *          For undirected graphs, it puts all vertex IDs in a bag
 *          such that the multiplicity of a vertex in the bag is the same as
 *          its degree. Then it draws pairs from the bag until the bag becomes
 *          empty. This method may generate both loop (self) edges and multiple
 *          edges. For directed graphs, the algorithm is basically the same,
 *          but two separate bags are used for the in- and out-degrees.
 *          Undirected graphs are generated with probability proportional to
 *          <code>(\prod_{i&lt;j} A_{ij} ! \prod_i A_{ii} !!)^{-1}</code>,
 *          where \c A denotes the adjacency matrix and <code>!!</code> denotes
 *          the double factorial. Here \c A is assumed to have twice the number of
 *          self-loops on its diagonal.
 *          The corresponding  expression for directed graphs is
 *          <code>(\prod_{i,j} A_{ij}!)^{-1}</code>.
 *          Thus the probability of all simple graphs (which only have 0s and 1s
 *          in the adjacency matrix) is the same, while that of
 *          non-simple ones depends on their edge and self-loop multiplicities.
 *          \cli IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE
 *          This method generates simple graphs.
 *          It is similar to \c IGRAPH_DEGSEQ_SIMPLE
 *          but tries to avoid multiple and loop edges and restarts the
 *          generation from scratch if it gets stuck. It can generate all simple
 *          realizations of a degree sequence, but it is not guaranteed
 *          to sample them uniformly. This method is relatively fast and it will
 *          eventually succeed if the provided degree sequence is graphical,
 *          but there is no upper bound on the number of iterations.
 *          \cli IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM
 *          This method is identical to \c IGRAPH_DEGSEQ_SIMPLE, but if the
 *          generated graph is not simple, it rejects it and re-starts the
 *          generation. It generates all simple graphs with the same probability.
 *          \cli IGRAPH_DEGSEQ_VL
 *          This method samples undirected connected graphs approximately
 *          uniformly. It is a Monte Carlo method based on degree-preserving
 *          edge swaps.
 *          This generator should be favoured if undirected and connected
 *          graphs are to be generated and execution time is not a concern.
 *          igraph uses the original implementation of Fabien Viger; for the algorithm,
 *          see https://www-complexnetworks.lip6.fr/~latapy/FV/generation.html
 *          and the paper https://arxiv.org/abs/cs/0502085
 *        \endclist
 * \return Error code:
 *          \c IGRAPH_ENOMEM: there is not enough
 *           memory to perform the operation.
 *          \c IGRAPH_EINVAL: invalid method parameter, or
 *           invalid in- and/or out-degree vectors. The degree vectors
 *           should be non-negative, \p out_deg should sum
 *           up to an even integer for undirected graphs; the length
 *           and sum of \p out_deg and
 *           \p in_deg
 *           should match for directed graphs.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number of edges
 *                  for \c IGRAPH_DEGSEQ_SIMPLE. The time complexity of the
 *                  other modes is not known.
 *
 * \sa \ref igraph_barabasi_game(), \ref igraph_erdos_renyi_game(),
 *     \ref igraph_is_degree_sequence(),
 *     \ref igraph_is_graphical_degree_sequence()
 *
 * \example examples/simple/igraph_degree_sequence_game.c
 */

int igraph_degree_sequence_game(igraph_t *graph, const igraph_vector_t *out_deg,
                                const igraph_vector_t *in_deg,
                                igraph_degseq_t method) {
    if (in_deg && igraph_vector_empty(in_deg) && !igraph_vector_empty(out_deg)) {
        in_deg = 0;
    }

    switch (method) {
    case IGRAPH_DEGSEQ_SIMPLE:
        return igraph_degree_sequence_game_simple(graph, out_deg, in_deg);

    case IGRAPH_DEGSEQ_VL:
        return igraph_degree_sequence_game_vl(graph, out_deg, in_deg);

    case IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE:
        if (in_deg == 0) {
            return igraph_degree_sequence_game_no_multiple_undirected(graph, out_deg);
        } else {
            return igraph_degree_sequence_game_no_multiple_directed(graph, out_deg, in_deg);
        }

    case IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM:
        if (in_deg == 0) {
            return igraph_degree_sequence_game_no_multiple_undirected_uniform(graph, out_deg);
        } else {
            return igraph_degree_sequence_game_no_multiple_directed_uniform(graph, out_deg, in_deg);
        }

    default:
        IGRAPH_ERROR("Invalid degree sequence game method", IGRAPH_EINVAL);
    }
}

/**
 * \ingroup generators
 * \function igraph_growing_random_game
 * \brief Generates a growing random graph.
 *
 * </para><para>
 * This function simulates a growing random graph. In each discrete
 * time step a new vertex is added and a number of new edges are also
 * added. These graphs are known to be different from standard (not
 * growing) random graphs.
 * \param graph Uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param m The number of edges to add in a time step (ie. after
 *        adding a vertex).
 * \param directed Boolean, whether to generate a directed graph.
 * \param citation Boolean, if \c TRUE, the edges always
 *        originate from the most recently added vertex.
 * \return Error code:
 *          \c IGRAPH_EINVAL: invalid
 *          \p n or \p m
 *          parameter.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges.
 *
 * \example examples/simple/igraph_growing_random_game.c
 */
int igraph_growing_random_game(igraph_t *graph, igraph_integer_t n,
                               igraph_integer_t m, igraph_bool_t directed,
                               igraph_bool_t citation) {

    long int no_of_nodes = n;
    long int no_of_neighbors = m;
    long int no_of_edges;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;

    long int resp = 0;

    long int i, j;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (m < 0) {
        IGRAPH_ERROR("Invalid number of edges per step (m)", IGRAPH_EINVAL);
    }

    no_of_edges = (no_of_nodes - 1) * no_of_neighbors;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    RNG_BEGIN();

    for (i = 1; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_neighbors; j++) {
            if (citation) {
                long int to = RNG_INTEGER(0, i - 1);
                VECTOR(edges)[resp++] = i;
                VECTOR(edges)[resp++] = to;
            } else {
                long int from = RNG_INTEGER(0, i);
                long int to = RNG_INTEGER(1, i);
                VECTOR(edges)[resp++] = from;
                VECTOR(edges)[resp++] = to;
            }
        }
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_callaway_traits_game
 * \brief Simulate a growing network with vertex types.
 *
 * </para><para>
 * The different types of vertices prefer to connect other types of
 * vertices with a given probability.</para><para>
 *
 * </para><para>
 * The simulation goes like this: in each discrete time step a new
 * vertex is added to the graph. The type of this vertex is generated
 * based on \p type_dist. Then two vertices are selected uniformly
 * randomly from the graph. The probability that they will be
 * connected depends on the types of these vertices and is taken from
 * \p pref_matrix. Then another two vertices are selected and this is
 * repeated \p edges_per_step times in each time step.
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of nodes in the graph.
 * \param types Number of node types.
 * \param edges_per_step The number of edges to be add per time step.
 * \param type_dist Vector giving the distribution of the vertex
 * types.
 * \param pref_matrix Matrix giving the connection probabilities for
 * the vertex types.
 * \param directed Logical, whether to generate a directed graph.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|e*log(|V|)), |V| is the number of vertices, e
 * is \p edges_per_step.
 */

int igraph_callaway_traits_game (igraph_t *graph, igraph_integer_t nodes,
                                 igraph_integer_t types, igraph_integer_t edges_per_step,
                                 igraph_vector_t *type_dist,
                                 igraph_matrix_t *pref_matrix,
                                 igraph_bool_t directed) {
    long int i, j;
    igraph_vector_t edges;
    igraph_vector_t cumdist;
    igraph_real_t maxcum;
    igraph_vector_t nodetypes;

    /* TODO: parameter checks */

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&nodetypes, nodes);

    VECTOR(cumdist)[0] = 0;
    for (i = 0; i < types; i++) {
        VECTOR(cumdist)[i + 1] = VECTOR(cumdist)[i] + VECTOR(*type_dist)[i];
    }
    maxcum = igraph_vector_tail(&cumdist);

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        igraph_real_t uni = RNG_UNIF(0, maxcum);
        long int type;
        igraph_vector_binsearch(&cumdist, uni, &type);
        VECTOR(nodetypes)[i] = type - 1;
    }

    for (i = 1; i < nodes; i++) {
        for (j = 0; j < edges_per_step; j++) {
            long int node1 = RNG_INTEGER(0, i);
            long int node2 = RNG_INTEGER(0, i);
            long int type1 = (long int) VECTOR(nodetypes)[node1];
            long int type2 = (long int) VECTOR(nodetypes)[node2];
            /*    printf("unif: %f, %f, types: %li, %li\n", uni1, uni2, type1, type2); */
            if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, node1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, node2));
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&nodetypes);
    igraph_vector_destroy(&cumdist);
    IGRAPH_FINALLY_CLEAN(2);
    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_establishment_game
 * \brief Generates a graph with a simple growing model with vertex types.
 *
 * </para><para>
 * The simulation goes like this: a single vertex is added at each
 * time step. This new vertex tries to connect to \p k vertices in the
 * graph. The probability that such a connection is realized depends
 * on the types of the vertices involved.
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param k The number of connections tried in each time step.
 * \param type_dist Vector giving the distribution of vertex types.
 * \param pref_matrix Matrix giving the connection probabilities for
 * different vertex types.
 * \param directed Logical, whether to generate a directed graph.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|*k*log(|V|)), |V| is the number of vertices
 * and k is the \p k parameter.
 */

int igraph_establishment_game(igraph_t *graph, igraph_integer_t nodes,
                              igraph_integer_t types, igraph_integer_t k,
                              igraph_vector_t *type_dist,
                              igraph_matrix_t *pref_matrix,
                              igraph_bool_t directed) {

    long int i, j;
    igraph_vector_t edges;
    igraph_vector_t cumdist;
    igraph_vector_t potneis;
    igraph_real_t maxcum;
    igraph_vector_t nodetypes;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&potneis, k);
    IGRAPH_VECTOR_INIT_FINALLY(&nodetypes, nodes);

    VECTOR(cumdist)[0] = 0;
    for (i = 0; i < types; i++) {
        VECTOR(cumdist)[i + 1] = VECTOR(cumdist)[i] + VECTOR(*type_dist)[i];
    }
    maxcum = igraph_vector_tail(&cumdist);

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        igraph_real_t uni = RNG_UNIF(0, maxcum);
        long int type;
        igraph_vector_binsearch(&cumdist, uni, &type);
        VECTOR(nodetypes)[i] = type - 1;
    }

    for (i = k; i < nodes; i++) {
        long int type1 = (long int) VECTOR(nodetypes)[i];
        igraph_random_sample(&potneis, 0, i - 1, k);
        for (j = 0; j < k; j++) {
            long int type2 = (long int) VECTOR(nodetypes)[(long int)VECTOR(potneis)[j]];
            if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(potneis)[j]));
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&nodetypes);
    igraph_vector_destroy(&potneis);
    igraph_vector_destroy(&cumdist);
    IGRAPH_FINALLY_CLEAN(3);
    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_recent_degree_game
 * \brief Stochastic graph generator based on the number of incident edges a node has gained recently
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph, this is the same as
 *        the number of time steps.
 * \param power The exponent, the probability that a node gains a
 *        new edge is proportional to the number of edges it has
 *        gained recently (in the last \p window time steps) to \p
 *        power.
 * \param window Integer constant, the size of the time window to use
 *        to count the number of recent edges.
 * \param m Integer constant, the number of edges to add per time
 *        step if the \p outseq parameter is a null pointer or a
 *        zero-length vector.
 * \param outseq The number of edges to add in each time step. This
 *        argument is ignored if it is a null pointer or a zero length
 *        vector, is this case the constant \p m parameter is used.
 * \param outpref Logical constant, if true the edges originated by a
 *        vertex also count as recent incident edges. It is false in
 *        most cases.
 * \param zero_appeal Constant giving the attractiveness of the
 *        vertices which haven't gained any edge recently.
 * \param directed Logical constant, whether to generate a directed
 *        graph.
 * \return Error code.
 *
 * Time complexity: O(|V|*log(|V|)+|E|), |V| is the number of
 * vertices, |E| is the number of edges in the graph.
 *
 */

int igraph_recent_degree_game(igraph_t *graph, igraph_integer_t n,
                              igraph_real_t power,
                              igraph_integer_t window,
                              igraph_integer_t m,
                              const igraph_vector_t *outseq,
                              igraph_bool_t outpref,
                              igraph_real_t zero_appeal,
                              igraph_bool_t directed) {

    long int no_of_nodes = n;
    long int no_of_neighbors = m;
    long int no_of_edges;
    igraph_vector_t edges;
    long int i, j;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;
    long int time_window = window;
    igraph_dqueue_t history;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != n) {
        IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
        IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
    }

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
    IGRAPH_CHECK(igraph_dqueue_init(&history,
                                    time_window * (no_of_neighbors + 1) + 10));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &history);

    RNG_BEGIN();

    /* first node */
    igraph_psumtree_update(&sumtree, 0, zero_appeal);
    igraph_dqueue_push(&history, -1);

    /* and the rest */
    for (i = 1; i < no_of_nodes; i++) {
        igraph_real_t sum;
        long int to;
        if (outseq != 0 && igraph_vector_size(outseq) != 0) {
            no_of_neighbors = (long int) VECTOR(*outseq)[i];
        }

        if (i >= time_window) {
            while ((j = (long int) igraph_dqueue_pop(&history)) != -1) {
                VECTOR(degree)[j] -= 1;
                igraph_psumtree_update(&sumtree, j,
                                       pow(VECTOR(degree)[j], power) + zero_appeal);
            }
        }

        sum = igraph_psumtree_sum(&sumtree);
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
            igraph_dqueue_push(&history, to);
        }
        igraph_dqueue_push(&history, -1);

        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            long int nn = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
            igraph_psumtree_update(&sumtree, nn,
                                   pow(VECTOR(degree)[nn], power) + zero_appeal);
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            igraph_psumtree_update(&sumtree, i,
                                   pow(VECTOR(degree)[i], power) + zero_appeal);
        } else {
            igraph_psumtree_update(&sumtree, i, zero_appeal);
        }
    }

    RNG_END();

    igraph_dqueue_destroy(&history);
    igraph_psumtree_destroy(&sumtree);
    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_barabasi_aging_game
 * \brief Preferential attachment with aging of vertices
 *
 * </para><para>
 * In this game, the probability that a node gains a new edge is
 * given by its (in-)degree (k) and age (l). This probability has a
 * degree dependent component multiplied by an age dependent
 * component. The degree dependent part is: \p deg_coef times k to the
 * power of \p pa_exp plus \p zero_deg_appeal; and the age dependent
 * part is \p age_coef times l to the power of \p aging_exp plus \p
 * zero_age_appeal.
 *
 * </para><para>
 * The age is based on the number of vertices in the
 * network and the \p aging_bin argument: vertices grew one unit older
 * after each \p aging_bin vertices added to the network.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the graph.
 * \param m The number of edges to add in each time step. If the \p
 *        outseq argument is not a null vector and not a zero-length
 *        vector.
 * \param outseq The number of edges to add in each time step. If it
 *        is a null pointer or a zero-length vector then it is ignored
 *        and the \p m argument is used instead.
 * \param outpref Logical constant, whether the edges
 *        initiated by a vertex contribute to the probability to gain
 *        a new edge.
 * \param pa_exp The exponent of the preferential attachment, a small
 *        positive number usually, the value 1 yields the classic
 *        linear preferential attachment.
 * \param aging_exp The exponent of the aging, this is a negative
 *        number usually.
 * \param aging_bin Integer constant, the number of vertices to add
 *        before vertices in the network grew one unit older.
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
 * Time complexity: O((|V|+|V|/aging_bin)*log(|V|)+|E|). |V| is the number
 * of vertices, |E| the number of edges.
 */

int igraph_barabasi_aging_game(igraph_t *graph,
                               igraph_integer_t nodes,
                               igraph_integer_t m,
                               const igraph_vector_t *outseq,
                               igraph_bool_t outpref,
                               igraph_real_t pa_exp,
                               igraph_real_t aging_exp,
                               igraph_integer_t aging_bin,
                               igraph_real_t zero_deg_appeal,
                               igraph_real_t zero_age_appeal,
                               igraph_real_t deg_coef,
                               igraph_real_t age_coef,
                               igraph_bool_t directed) {
    long int no_of_nodes = nodes;
    long int no_of_neighbors = m;
    long int binwidth = nodes / aging_bin + 1;
    long int no_of_edges;
    igraph_vector_t edges;
    long int i, j, k;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;

    if (no_of_nodes < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != no_of_nodes) {
        IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
        IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
    }
    if (aging_bin <= 0) {
        IGRAPH_ERROR("Invalid aging bin", IGRAPH_EINVAL);
    }

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
    igraph_psumtree_update(&sumtree, 0, zero_deg_appeal * (1 + zero_age_appeal));

    /* and the rest */
    for (i = 1; i < no_of_nodes; i++) {
        igraph_real_t sum;
        long int to;
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
            igraph_psumtree_update(&sumtree, n,
                                   (deg_coef * pow(VECTOR(degree)[n], pa_exp)
                                    + zero_deg_appeal)*
                                   (age_coef * pow(age + 1, aging_exp) + zero_age_appeal));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            igraph_psumtree_update(&sumtree, i, (zero_age_appeal + 1)*
                                   (deg_coef * pow(VECTOR(degree)[i], pa_exp)
                                    + zero_deg_appeal));
        } else {
            igraph_psumtree_update(&sumtree, i, (1 + zero_age_appeal)*zero_deg_appeal);
        }

        /* aging */
        for (k = 1; i - binwidth * k + 1 >= 1; k++) {
            long int shnode = i - binwidth * k;
            long int deg = (long int) VECTOR(degree)[shnode];
            long int age = (i - shnode) / binwidth;
            /* igraph_real_t old=igraph_psumtree_get(&sumtree, shnode); */
            igraph_psumtree_update(&sumtree, shnode,
                                   (deg_coef * pow(deg, pa_exp) + zero_deg_appeal) *
                                   (age_coef * pow(age + 2, aging_exp) + zero_age_appeal));
        }
    }

    RNG_END();

    igraph_vector_destroy(&degree);
    igraph_psumtree_destroy(&sumtree);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_recent_degree_aging_game
 * \brief Preferential attachment based on the number of edges gained recently, with aging of vertices
 *
 * </para><para>
 * This game is very similar to \ref igraph_barabasi_aging_game(),
 * except that instead of the total number of incident edges the
 * number of edges gained in the last \p time_window time steps are
 * counted.
 *
 * </para><para>The degree dependent part of the attractiveness is
 * given by k to the power of \p pa_exp plus \p zero_appeal; the age
 * dependent part is l to the power to \p aging_exp.
 * k is the number of edges gained in the last \p time_window time
 * steps, l is the age of the vertex.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the graph.
 * \param m The number of edges to add in each time step. If the \p
 *        outseq argument is not a null vector or a zero-length vector
 *        then it is ignored.
 * \param outseq Vector giving the number of edges to add in each time
 *        step. If it is a null pointer or a zero-length vector then
 *        it is ignored and the \p m argument is used.
 * \param outpref Logical constant, if true the edges initiated by a
 *        vertex are also counted. Normally it is false.
 * \param pa_exp The exponent for the preferential attachment.
 * \param aging_exp The exponent for the aging, normally it is
 *        negative: old vertices gain edges with less probability.
 * \param aging_bin Integer constant, gives the scale of the aging.
 *        The age of the vertices is incremented by one after every \p
 *        aging_bin vertex added.
 * \param time_window The time window to use to count the number of
 *        incident edges for the vertices.
 * \param zero_appeal The degree dependent part of the attractiveness
 *        for zero degree vertices.
 * \param directed Logical constant, whether to create a directed
 *        graph.
 * \return Error code.
 *
 * Time complexity: O((|V|+|V|/aging_bin)*log(|V|)+|E|). |V| is the number
 * of vertices, |E| the number of edges.
 */

int igraph_recent_degree_aging_game(igraph_t *graph,
                                    igraph_integer_t nodes,
                                    igraph_integer_t m,
                                    const igraph_vector_t *outseq,
                                    igraph_bool_t outpref,
                                    igraph_real_t pa_exp,
                                    igraph_real_t aging_exp,
                                    igraph_integer_t aging_bin,
                                    igraph_integer_t time_window,
                                    igraph_real_t zero_appeal,
                                    igraph_bool_t directed) {

    long int no_of_nodes = nodes;
    long int no_of_neighbors = m;
    long int binwidth = nodes / aging_bin + 1;
    long int no_of_edges;
    igraph_vector_t edges;
    long int i, j, k;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;
    igraph_dqueue_t history;

    if (no_of_nodes < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != no_of_nodes) {
        IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
        IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
    }
    if (aging_bin <= 0) {
        IGRAPH_ERROR("Invalid aging bin", IGRAPH_EINVAL);
    }

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
    IGRAPH_CHECK(igraph_dqueue_init(&history,
                                    time_window * (no_of_neighbors + 1) + 10));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &history);

    RNG_BEGIN();

    /* first node */
    igraph_psumtree_update(&sumtree, 0, zero_appeal);
    igraph_dqueue_push(&history, -1);

    /* and the rest */
    for (i = 1; i < no_of_nodes; i++) {
        igraph_real_t sum;
        long int to;
        if (outseq != 0 && igraph_vector_size(outseq) != 0) {
            no_of_neighbors = (long int) VECTOR(*outseq)[i];
        }

        if (i >= time_window) {
            while ((j = (long int) igraph_dqueue_pop(&history)) != -1) {
                long int age = (i - j) / binwidth;
                VECTOR(degree)[j] -= 1;
                igraph_psumtree_update(&sumtree, j,
                                       (pow(VECTOR(degree)[j], pa_exp) + zero_appeal)*
                                       pow(age + 1, aging_exp));
            }
        }

        sum = igraph_psumtree_sum(&sumtree);
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
            igraph_dqueue_push(&history, to);
        }
        igraph_dqueue_push(&history, -1);

        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            long int n = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
            long int age = (i - n) / binwidth;
            igraph_psumtree_update(&sumtree, n,
                                   (pow(VECTOR(degree)[n], pa_exp) + zero_appeal)*
                                   pow(age + 1, aging_exp));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            igraph_psumtree_update(&sumtree, i,
                                   pow(VECTOR(degree)[i], pa_exp) + zero_appeal);
        } else {
            igraph_psumtree_update(&sumtree, i, zero_appeal);
        }

        /* aging */
        for (k = 1; i - binwidth * k + 1 >= 1; k++) {
            long int shnode = i - binwidth * k;
            long int deg = (long int) VECTOR(degree)[shnode];
            long int age = (i - shnode) / binwidth;
            igraph_psumtree_update(&sumtree, shnode,
                                   (pow(deg, pa_exp) + zero_appeal) *
                                   pow(age + 2, aging_exp));
        }
    }

    RNG_END();

    igraph_dqueue_destroy(&history);
    igraph_vector_destroy(&degree);
    igraph_psumtree_destroy(&sumtree);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_grg_game
 * \brief Generating geometric random graphs.
 *
 * A geometric random graph is created by dropping points (=vertices)
 * randomly to the unit square and then connecting all those pairs
 * which are less than \c radius apart in Euclidean norm.
 *
 * </para><para>
 * Original code contributed by Keith Briggs, thanks Keith.
 * \param graph Pointer to an uninitialized graph object,
 * \param nodes The number of vertices in the graph.
 * \param radius The radius within which the vertices will be connected.
 * \param torus Logical constant, if true periodic boundary conditions
 *        will be used, ie. the vertices are assumed to be on a torus
 *        instead of a square.
 * \return Error code.
 *
 * Time complexity: TODO, less than O(|V|^2+|E|).
 *
 * \example examples/simple/igraph_grg_game.c
 */

int igraph_grg_game(igraph_t *graph, igraph_integer_t nodes,
                    igraph_real_t radius, igraph_bool_t torus,
                    igraph_vector_t *x, igraph_vector_t *y) {

    long int i;
    igraph_vector_t myx, myy, *xx = &myx, *yy = &myy, edges;
    igraph_real_t r2 = radius * radius;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes));

    if (x) {
        xx = x;
        IGRAPH_CHECK(igraph_vector_resize(xx, nodes));
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(xx, nodes);
    }
    if (y) {
        yy = y;
        IGRAPH_CHECK(igraph_vector_resize(yy, nodes));
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(yy, nodes);
    }

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        VECTOR(*xx)[i] = RNG_UNIF01();
        VECTOR(*yy)[i] = RNG_UNIF01();
    }

    RNG_END();

    igraph_vector_sort(xx);

    if (!torus) {
        for (i = 0; i < nodes; i++) {
            igraph_real_t xx1 = VECTOR(*xx)[i];
            igraph_real_t yy1 = VECTOR(*yy)[i];
            long int j = i + 1;
            igraph_real_t dx, dy;
            while ( j < nodes && (dx = VECTOR(*xx)[j] - xx1) < radius) {
                dy = VECTOR(*yy)[j] - yy1;
                if (dx * dx + dy * dy < r2) {
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                }
                j++;
            }
        }
    } else {
        for (i = 0; i < nodes; i++) {
            igraph_real_t xx1 = VECTOR(*xx)[i];
            igraph_real_t yy1 = VECTOR(*yy)[i];
            long int j = i + 1;
            igraph_real_t dx, dy;
            while ( j < nodes && (dx = VECTOR(*xx)[j] - xx1) < radius) {
                dy = fabs(VECTOR(*yy)[j] - yy1);
                if (dx > 0.5) {
                    dx = 1 - dx;
                }
                if (dy > 0.5) {
                    dy = 1 - dy;
                }
                if (dx * dx + dy * dy < r2) {
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                }
                j++;
            }
            if (j == nodes) {
                j = 0;
                while (j < i && (dx = 1 - xx1 + VECTOR(*xx)[j]) < radius &&
                       xx1 - VECTOR(*xx)[j] >= radius) {
                    dy = fabs(VECTOR(*yy)[j] - yy1);
                    if (dy > 0.5) {
                        dy = 1 - dy;
                    }
                    if (dx * dx + dy * dy < r2) {
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                    }
                    j++;
                }
            }
        }
    }

    if (!y) {
        igraph_vector_destroy(yy);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (!x) {
        igraph_vector_destroy(xx);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, IGRAPH_UNDIRECTED));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}


static void igraph_i_preference_game_free_vids_by_type(igraph_vector_ptr_t *vecs) {
    int i = 0, n;
    igraph_vector_t *v;

    n = (int) igraph_vector_ptr_size(vecs);
    for (i = 0; i < n; i++) {
        v = (igraph_vector_t*)VECTOR(*vecs)[i];
        if (v) {
            igraph_vector_destroy(v);
        }
    }
    igraph_vector_ptr_destroy_all(vecs);
}

/**
 * \function igraph_preference_game
 * \brief Generates a graph with vertex types and connection preferences
 *
 * </para><para>
 * This is practically the nongrowing variant of \ref
 * igraph_establishment_game. A given number of vertices are
 * generated. Every vertex is assigned to a vertex type according to
 * the given type probabilities. Finally, every
 * vertex pair is evaluated and an edge is created between them with a
 * probability depending on the types of the vertices involved.
 *
 * </para><para>
 * In other words, this function generates a graph according to a
 * block-model. Vertices are divided into groups (or blocks), and
 * the probability the two vertices are connected depends on their
 * groups only.
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param type_dist Vector giving the distribution of vertex types. If
 *   \c NULL, all vertex types will have equal probability. See also the
 *   \c fixed_sizes argument.
 * \param fixed_sizes Boolean. If true, then the number of vertices with a
 *   given vertex type is fixed and the \c type_dist argument gives these
 *   numbers for each vertex type. If true, and \c type_dist is \c NULL,
 *   then the function tries to make vertex groups of the same size. If this
 *   is not possible, then some groups will have an extra vertex.
 * \param pref_matrix Matrix giving the connection probabilities for
 *   different vertex types. This should be symmetric if the requested
 *   graph is undirected.
 * \param node_type_vec A vector where the individual generated vertex types
 *   will be stored. If \c NULL , the vertex types won't be saved.
 * \param directed Logical, whether to generate a directed graph. If undirected
 *   graphs are requested, only the lower left triangle of the preference
 *   matrix is considered.
 * \param loops Logical, whether loop edges are allowed.
 * \return Error code.
 *
 * Added in version 0.3.</para><para>
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa igraph_establishment_game()
 *
 * \example examples/simple/igraph_preference_game.c
 */

int igraph_preference_game(igraph_t *graph, igraph_integer_t nodes,
                           igraph_integer_t types,
                           const igraph_vector_t *type_dist,
                           igraph_bool_t fixed_sizes,
                           const igraph_matrix_t *pref_matrix,
                           igraph_vector_t *node_type_vec,
                           igraph_bool_t directed,
                           igraph_bool_t loops) {

    long int i, j;
    igraph_vector_t edges, s;
    igraph_vector_t* nodetypes;
    igraph_vector_ptr_t vids_by_type;
    igraph_real_t maxcum, maxedges;

    if (types < 1) {
        IGRAPH_ERROR("types must be >= 1", IGRAPH_EINVAL);
    }
    if (nodes < 0) {
        IGRAPH_ERROR("nodes must be >= 0", IGRAPH_EINVAL);
    }
    if (type_dist && igraph_vector_size(type_dist) != types) {
        if (igraph_vector_size(type_dist) > types) {
            IGRAPH_WARNING("length of type_dist > types, type_dist will be trimmed");
        } else {
            IGRAPH_ERROR("type_dist vector too short", IGRAPH_EINVAL);
        }
    }
    if (igraph_matrix_nrow(pref_matrix) < types ||
        igraph_matrix_ncol(pref_matrix) < types) {
        IGRAPH_ERROR("pref_matrix too small", IGRAPH_EINVAL);
    }

    if (fixed_sizes && type_dist) {
        if (igraph_vector_sum(type_dist) != nodes) {
            IGRAPH_ERROR("Invalid group sizes, their sum must match the number"
                         " of vertices", IGRAPH_EINVAL);
        }
    }

    if (node_type_vec) {
        IGRAPH_CHECK(igraph_vector_resize(node_type_vec, nodes));
        nodetypes = node_type_vec;
    } else {
        nodetypes = igraph_Calloc(1, igraph_vector_t);
        if (nodetypes == 0) {
            IGRAPH_ERROR("preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, nodetypes);
        IGRAPH_VECTOR_INIT_FINALLY(nodetypes, nodes);
    }

    IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_type, types));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_type);
    for (i = 0; i < types; i++) {
        VECTOR(vids_by_type)[i] = igraph_Calloc(1, igraph_vector_t);
        if (VECTOR(vids_by_type)[i] == 0) {
            IGRAPH_ERROR("preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_type)[i], 0));
    }
    IGRAPH_FINALLY_CLEAN(1);   /* removing igraph_vector_ptr_destroy_all */
    IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_type);

    RNG_BEGIN();

    if (!fixed_sizes) {

        igraph_vector_t cumdist;
        IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types + 1);

        VECTOR(cumdist)[0] = 0;
        if (type_dist) {
            for (i = 0; i < types; i++) {
                VECTOR(cumdist)[i + 1] = VECTOR(cumdist)[i] + VECTOR(*type_dist)[i];
            }
        } else {
            for (i = 0; i < types; i++) {
                VECTOR(cumdist)[i + 1] = i + 1;
            }
        }
        maxcum = igraph_vector_tail(&cumdist);

        for (i = 0; i < nodes; i++) {
            long int type1;
            igraph_real_t uni1 = RNG_UNIF(0, maxcum);
            igraph_vector_binsearch(&cumdist, uni1, &type1);
            VECTOR(*nodetypes)[i] = type1 - 1;
            IGRAPH_CHECK(igraph_vector_push_back(
                             (igraph_vector_t*)VECTOR(vids_by_type)[type1 - 1], i));
        }

        igraph_vector_destroy(&cumdist);
        IGRAPH_FINALLY_CLEAN(1);

    } else {

        int an = 0;
        if (type_dist) {
            for (i = 0; i < types; i++) {
                int no = (int) VECTOR(*type_dist)[i];
                igraph_vector_t *v = VECTOR(vids_by_type)[i];
                for (j = 0; j < no && an < nodes; j++) {
                    VECTOR(*nodetypes)[an] = i;
                    IGRAPH_CHECK(igraph_vector_push_back(v, an));
                    an++;
                }
            }
        } else {
            int fixno = (int) ceil( (double)nodes / types);
            for (i = 0; i < types; i++) {
                igraph_vector_t *v = VECTOR(vids_by_type)[i];
                for (j = 0; j < fixno && an < nodes; j++) {
                    VECTOR(*nodetypes)[an++] = i;
                    IGRAPH_CHECK(igraph_vector_push_back(v, an));
                    an++;
                }
            }
        }

    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);

    for (i = 0; i < types; i++) {
        for (j = 0; j < types; j++) {
            /* Generating the random subgraph between vertices of type i and j */
            long int k, l;
            igraph_real_t p, last;
            igraph_vector_t *v1, *v2;
            long int v1_size, v2_size;

            IGRAPH_ALLOW_INTERRUPTION();

            v1 = (igraph_vector_t*)VECTOR(vids_by_type)[i];
            v2 = (igraph_vector_t*)VECTOR(vids_by_type)[j];
            v1_size = igraph_vector_size(v1);
            v2_size = igraph_vector_size(v2);

            p = MATRIX(*pref_matrix, i, j);
            igraph_vector_clear(&s);
            if (i != j) {
                /* The two vertex sets are disjoint, this is the easier case */
                if (i > j && !directed) {
                    continue;
                }
                maxedges = v1_size * v2_size;
            } else {
                if (directed && loops) {
                    maxedges = v1_size * v1_size;
                } else if (directed && !loops) {
                    maxedges = v1_size * (v1_size - 1);
                } else if (!directed && loops) {
                    maxedges = v1_size * (v1_size + 1) / 2;
                } else {
                    maxedges = v1_size * (v1_size - 1) / 2;
                }
            }

            IGRAPH_CHECK(igraph_vector_reserve(&s, (long int) (maxedges * p * 1.1)));

            last = RNG_GEOM(p);
            while (last < maxedges) {
                IGRAPH_CHECK(igraph_vector_push_back(&s, last));
                last += RNG_GEOM(p);
                last += 1;
            }
            l = igraph_vector_size(&s);

            IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&edges) + l * 2));

            if (i != j) {
                /* Generating the subgraph between vertices of type i and j */
                for (k = 0; k < l; k++) {
                    long int to = (long int) floor(VECTOR(s)[k] / v1_size);
                    long int from = (long int) (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                    igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
                }
            } else {
                /* Generating the subgraph among vertices of type i */
                if (directed && loops) {
                    for (k = 0; k < l; k++) {
                        long int to = (long int) floor(VECTOR(s)[k] / v1_size);
                        long int from = (long int) (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else if (directed && !loops) {
                    for (k = 0; k < l; k++) {
                        long int to = (long int) floor(VECTOR(s)[k] / v1_size);
                        long int from = (long int) (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                        if (from == to) {
                            to = v1_size - 1;
                        }
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else if (!directed && loops) {
                    for (k = 0; k < l; k++) {
                        long int to = (long int) floor((sqrt(8 * VECTOR(s)[k] + 1) - 1) / 2);
                        long int from = (long int) (VECTOR(s)[k] - (((igraph_real_t)to) * (to + 1)) / 2);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else {
                    for (k = 0; k < l; k++) {
                        long int to = (long int) floor((sqrt(8 * VECTOR(s)[k] + 1) + 1) / 2);
                        long int from = (long int) (VECTOR(s)[k] - (((igraph_real_t)to) * (to - 1)) / 2);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                }
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&s);
    igraph_i_preference_game_free_vids_by_type(&vids_by_type);
    IGRAPH_FINALLY_CLEAN(2);

    if (node_type_vec == 0) {
        igraph_vector_destroy(nodetypes);
        igraph_Free(nodetypes);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_asymmetric_preference_game
 * \brief Generates a graph with asymmetric vertex types and connection preferences
 *
 * </para><para>
 * This is the asymmetric variant of \ref igraph_preference_game() .
 * A given number of vertices are generated. Every vertex is assigned to an
 * "incoming" and an "outgoing" vertex type according to the given joint
 * type probabilities. Finally, every vertex pair is evaluated and a
 * directed edge is created between them with a probability depending on the
 * "outgoing" type of the source vertex and the "incoming" type of the target
 * vertex.
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param type_dist_matrix Matrix giving the joint distribution of vertex types.
 *   If null, incoming and outgoing vertex types are independent and uniformly
 *   distributed.
 * \param pref_matrix Matrix giving the connection probabilities for
 *   different vertex types.
 * \param node_type_in_vec A vector where the individual generated "incoming"
 *   vertex types will be stored. If NULL, the vertex types won't be saved.
 * \param node_type_out_vec A vector where the individual generated "outgoing"
 *   vertex types will be stored. If NULL, the vertex types won't be saved.
 * \param loops Logical, whether loop edges are allowed.
 * \return Error code.
 *
 * Added in version 0.3.</para><para>
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_preference_game()
 */

int igraph_asymmetric_preference_game(igraph_t *graph, igraph_integer_t nodes,
                                      igraph_integer_t types,
                                      igraph_matrix_t *type_dist_matrix,
                                      igraph_matrix_t *pref_matrix,
                                      igraph_vector_t *node_type_in_vec,
                                      igraph_vector_t *node_type_out_vec,
                                      igraph_bool_t loops) {

    long int i, j, k;
    igraph_vector_t edges, cumdist, s, intersect;
    igraph_vector_t *nodetypes_in;
    igraph_vector_t *nodetypes_out;
    igraph_vector_ptr_t vids_by_intype, vids_by_outtype;
    igraph_real_t maxcum, maxedges;

    if (types < 1) {
        IGRAPH_ERROR("types must be >= 1", IGRAPH_EINVAL);
    }
    if (nodes < 0) {
        IGRAPH_ERROR("nodes must be >= 0", IGRAPH_EINVAL);
    }
    if (type_dist_matrix) {
        if (igraph_matrix_nrow(type_dist_matrix) < types ||
            igraph_matrix_ncol(type_dist_matrix) < types) {
            IGRAPH_ERROR("type_dist_matrix too small", IGRAPH_EINVAL);
        } else if (igraph_matrix_nrow(type_dist_matrix) > types ||
                   igraph_matrix_ncol(type_dist_matrix) > types) {
            IGRAPH_WARNING("type_dist_matrix will be trimmed");
        }
    }
    if (igraph_matrix_nrow(pref_matrix) < types ||
        igraph_matrix_ncol(pref_matrix) < types) {
        IGRAPH_ERROR("pref_matrix too small", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types * types + 1);

    if (node_type_in_vec) {
        nodetypes_in = node_type_in_vec;
        IGRAPH_CHECK(igraph_vector_resize(nodetypes_in, nodes));
    } else {
        nodetypes_in = igraph_Calloc(1, igraph_vector_t);
        if (nodetypes_in == 0) {
            IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_VECTOR_INIT_FINALLY(nodetypes_in, nodes);
    }

    if (node_type_out_vec) {
        nodetypes_out = node_type_out_vec;
        IGRAPH_CHECK(igraph_vector_resize(nodetypes_out, nodes));
    } else {
        nodetypes_out = igraph_Calloc(1, igraph_vector_t);
        if (nodetypes_out == 0) {
            IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_VECTOR_INIT_FINALLY(nodetypes_out, nodes);
    }

    IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_intype, types));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_intype);
    IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_outtype, types));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_outtype);
    for (i = 0; i < types; i++) {
        VECTOR(vids_by_intype)[i] = igraph_Calloc(1, igraph_vector_t);
        VECTOR(vids_by_outtype)[i] = igraph_Calloc(1, igraph_vector_t);
        if (VECTOR(vids_by_intype)[i] == 0 || VECTOR(vids_by_outtype)[i] == 0) {
            IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_intype)[i], 0));
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_outtype)[i], 0));
    }
    IGRAPH_FINALLY_CLEAN(2);   /* removing igraph_vector_ptr_destroy_all */
    IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_intype);
    IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_outtype);

    VECTOR(cumdist)[0] = 0;
    if (type_dist_matrix) {
        for (i = 0, k = 0; i < types; i++) {
            for (j = 0; j < types; j++, k++) {
                VECTOR(cumdist)[k + 1] = VECTOR(cumdist)[k] + MATRIX(*type_dist_matrix, i, j);
            }
        }
    } else {
        for (i = 0; i < types * types; i++) {
            VECTOR(cumdist)[i + 1] = i + 1;
        }
    }
    maxcum = igraph_vector_tail(&cumdist);

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        long int type1, type2;
        igraph_real_t uni1 = RNG_UNIF(0, maxcum);
        igraph_vector_binsearch(&cumdist, uni1, &type1);
        type2 = (type1 - 1) % (int)types;
        type1 = (type1 - 1) / (int)types;
        VECTOR(*nodetypes_in)[i] = type1;
        VECTOR(*nodetypes_out)[i] = type2;
        IGRAPH_CHECK(igraph_vector_push_back(
                         (igraph_vector_t*)VECTOR(vids_by_intype)[type1], i));
        IGRAPH_CHECK(igraph_vector_push_back(
                         (igraph_vector_t*)VECTOR(vids_by_outtype)[type2], i));
    }

    igraph_vector_destroy(&cumdist);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&intersect, 0);
    for (i = 0; i < types; i++) {
        for (j = 0; j < types; j++) {
            long int kk, l, c;
            igraph_real_t p, last;
            igraph_vector_t *v1, *v2;
            long int v1_size, v2_size;

            IGRAPH_ALLOW_INTERRUPTION();

            v1 = (igraph_vector_t*)VECTOR(vids_by_outtype)[i];
            v2 = (igraph_vector_t*)VECTOR(vids_by_intype)[j];
            v1_size = igraph_vector_size(v1);
            v2_size = igraph_vector_size(v2);

            maxedges = v1_size * v2_size;
            if (!loops) {
                IGRAPH_CHECK(igraph_vector_intersect_sorted(v1, v2, &intersect));
                c = igraph_vector_size(&intersect);
                maxedges -= c;
            }

            p = MATRIX(*pref_matrix, i, j);
            igraph_vector_clear(&s);
            IGRAPH_CHECK(igraph_vector_reserve(&s, (long int) (maxedges * p * 1.1)));

            last = RNG_GEOM(p);
            while (last < maxedges) {
                IGRAPH_CHECK(igraph_vector_push_back(&s, last));
                last += RNG_GEOM(p);
                last += 1;
            }
            l = igraph_vector_size(&s);

            IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&edges) + l * 2));

            if (!loops && c > 0) {
                for (kk = 0; kk < l; kk++) {
                    long int to = (long int) floor(VECTOR(s)[kk] / v1_size);
                    long int from = (long int) (VECTOR(s)[kk] - ((igraph_real_t)to) * v1_size);
                    if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
                        /* remap loop edges */
                        to = v2_size - 1;
                        igraph_vector_binsearch(&intersect, VECTOR(*v1)[from], &c);
                        from = v1_size - 1;
                        if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
                            from--;
                        }
                        while (c > 0) {
                            c--; from--;
                            if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
                                from--;
                            }
                        }
                    }
                    igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
                }
            } else {
                for (kk = 0; kk < l; kk++) {
                    long int to = (long int) floor(VECTOR(s)[kk] / v1_size);
                    long int from = (long int) (VECTOR(s)[kk] - ((igraph_real_t)to) * v1_size);
                    igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
                }
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&s);
    igraph_vector_destroy(&intersect);
    igraph_i_preference_game_free_vids_by_type(&vids_by_intype);
    igraph_i_preference_game_free_vids_by_type(&vids_by_outtype);
    IGRAPH_FINALLY_CLEAN(4);

    if (node_type_out_vec == 0) {
        igraph_vector_destroy(nodetypes_out);
        igraph_Free(nodetypes_out);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (node_type_in_vec == 0) {
        igraph_vector_destroy(nodetypes_in);
        igraph_Free(nodetypes_in);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, 1));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

static int igraph_i_rewire_edges_no_multiple(igraph_t *graph, igraph_real_t prob,
                                             igraph_bool_t loops,
                                             igraph_vector_t *edges) {

    int no_verts = igraph_vcount(graph);
    int no_edges = igraph_ecount(graph);
    igraph_vector_t eorder, tmp;
    igraph_vector_int_t first, next, prev, marked;
    int i, to_rewire, last_other = -1;

    /* Create our special graph representation */

# define ADD_STUB(vertex, stub) do {                \
        if (VECTOR(first)[(vertex)]) {              \
            VECTOR(prev)[(int) VECTOR(first)[(vertex)]-1]=(stub)+1;   \
        }                               \
        VECTOR(next)[(stub)]=VECTOR(first)[(vertex)];       \
        VECTOR(prev)[(stub)]=0;                 \
        VECTOR(first)[(vertex)]=(stub)+1;               \
    } while (0)

# define DEL_STUB(vertex, stub) do {                    \
        if (VECTOR(next)[(stub)]) {                     \
            VECTOR(prev)[VECTOR(next)[(stub)]-1]=VECTOR(prev)[(stub)];    \
        }                                   \
        if (VECTOR(prev)[(stub)]) {                     \
            VECTOR(next)[VECTOR(prev)[(stub)]-1]=VECTOR(next)[(stub)];    \
        } else {                                \
            VECTOR(first)[(vertex)]=VECTOR(next)[(stub)];         \
        }                                   \
    } while (0)

# define MARK_NEIGHBORS(vertex) do {                \
        int xxx_ =VECTOR(first)[(vertex)];              \
        while (xxx_) {                      \
            int o= (int) VECTOR(*edges)[xxx_ % 2 ? xxx_ : xxx_-2];    \
            VECTOR(marked)[o]=other+1;                \
            xxx_=VECTOR(next)[xxx_-1];                \
        }                               \
    } while (0)

    IGRAPH_CHECK(igraph_vector_int_init(&first, no_verts));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &first);
    IGRAPH_CHECK(igraph_vector_int_init(&next, no_edges * 2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &next);
    IGRAPH_CHECK(igraph_vector_int_init(&prev, no_edges * 2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &prev);
    IGRAPH_CHECK(igraph_get_edgelist(graph, edges, /*bycol=*/ 0));
    IGRAPH_VECTOR_INIT_FINALLY(&eorder, no_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_edges);
    for (i = 0; i < no_edges; i++) {
        int idx1 = 2 * i, idx2 = idx1 + 1,
            from = (int) VECTOR(*edges)[idx1], to = (int) VECTOR(*edges)[idx2];
        VECTOR(tmp)[i] = from;
        ADD_STUB(from, idx1);
        ADD_STUB(to, idx2);
    }
    IGRAPH_CHECK(igraph_vector_order1(&tmp, &eorder, no_verts));
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_vector_int_init(&marked, no_verts));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &marked);

    /* Rewire the stubs, part I */

    to_rewire = (int) RNG_GEOM(prob);
    while (to_rewire < no_edges) {
        int stub = (int) (2 * VECTOR(eorder)[to_rewire] + 1);
        int v = (int) VECTOR(*edges)[stub];
        int ostub = stub - 1;
        int other = (int) VECTOR(*edges)[ostub];
        int pot;
        if (last_other != other) {
            MARK_NEIGHBORS(other);
        }
        /* Do the rewiring */
        do {
            if (loops) {
                pot = (int) RNG_INTEGER(0, no_verts - 1);
            } else {
                pot = (int) RNG_INTEGER(0, no_verts - 2);
                pot = pot != other ? pot : no_verts - 1;
            }
        } while (VECTOR(marked)[pot] == other + 1 && pot != v);

        if (pot != v) {
            DEL_STUB(v, stub);
            ADD_STUB(pot, stub);
            VECTOR(marked)[v] = 0;
            VECTOR(marked)[pot] = other + 1;
            VECTOR(*edges)[stub] = pot;
        }

        to_rewire += RNG_GEOM(prob) + 1;
        last_other = other;
    }

    /* Create the new index, from the potentially rewired stubs */

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_edges);
    for (i = 0; i < no_edges; i++) {
        VECTOR(tmp)[i] = VECTOR(*edges)[2 * i + 1];
    }
    IGRAPH_CHECK(igraph_vector_order1(&tmp, &eorder, no_verts));
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    /* Rewire the stubs, part II */

    igraph_vector_int_null(&marked);
    last_other = -1;

    to_rewire = (int) RNG_GEOM(prob);
    while (to_rewire < no_edges) {
        int stub = (int) (2 * VECTOR(eorder)[to_rewire]);
        int v = (int) VECTOR(*edges)[stub];
        int ostub = stub + 1;
        int other = (int) VECTOR(*edges)[ostub];
        int pot;
        if (last_other != other) {
            MARK_NEIGHBORS(other);
        }
        /* Do the rewiring */
        do {
            if (loops) {
                pot = (int) RNG_INTEGER(0, no_verts - 1);
            } else {
                pot = (int) RNG_INTEGER(0, no_verts - 2);
                pot = pot != other ? pot : no_verts - 1;
            }
        } while (VECTOR(marked)[pot] == other + 1 && pot != v);
        if (pot != v) {
            DEL_STUB(v, stub);
            ADD_STUB(pot, stub);
            VECTOR(marked)[v] = 0;
            VECTOR(marked)[pot] = other + 1;
            VECTOR(*edges)[stub] = pot;
        }

        to_rewire += RNG_GEOM(prob) + 1;
        last_other = other;
    }

    igraph_vector_int_destroy(&marked);
    igraph_vector_int_destroy(&prev);
    igraph_vector_int_destroy(&next);
    igraph_vector_int_destroy(&first);
    igraph_vector_destroy(&eorder);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

#undef ADD_STUB
#undef DEL_STUB
#undef MARK_NEIGHBORS

/**
 * \function igraph_rewire_edges
 * \brief Rewire the edges of a graph with constant probability
 *
 * This function rewires the edges of a graph with a constant
 * probability. More precisely each end point of each edge is rewired
 * to a uniformly randomly chosen vertex with constant probability \p
 * prob.
 *
 * </para><para> Note that this function modifies the input \p graph,
 * call \ref igraph_copy() if you want to keep it.
 *
 * \param graph The input graph, this will be rewired, it can be
 *    directed or undirected.
 * \param prob The rewiring probability a constant between zero and
 *    one (inclusive).
 * \param loops Boolean, whether loop edges are allowed in the new
 *    graph, or not.
 * \param multiple Boolean, whether multiple edges are allowed in the
 *    new graph.
 * \return Error code.
 *
 * \sa \ref igraph_watts_strogatz_game() uses this function for the
 * rewiring.
 *
 * Time complexity: O(|V|+|E|).
 */

int igraph_rewire_edges(igraph_t *graph, igraph_real_t prob,
                        igraph_bool_t loops, igraph_bool_t multiple) {

    igraph_t newgraph;
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_nodes = igraph_vcount(graph);
    long int endpoints = no_of_edges * 2;
    long int to_rewire;
    igraph_vector_t edges;

    if (prob < 0 || prob > 1) {
        IGRAPH_ERROR("Rewiring probability should be between zero and one",
                     IGRAPH_EINVAL);
    }

    if (prob == 0) {
        /* This is easy, just leave things as they are */
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, endpoints);

    RNG_BEGIN();

    if (prob != 0 && no_of_edges > 0) {
        if (multiple) {
            /* If multiple edges are allowed, then there is an easy and fast
            method. Each endpoint of an edge is rewired with probability p,
             so the "skips" between the really rewired endpoints follow a
             geometric distribution. */
            IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
            to_rewire = (long int) RNG_GEOM(prob);
            while (to_rewire < endpoints) {
                if (loops) {
                    VECTOR(edges)[to_rewire] = RNG_INTEGER(0, no_of_nodes - 1);
                } else {
                    long int opos = to_rewire % 2 ? to_rewire - 1 : to_rewire + 1;
                    long int nei = (long int) VECTOR(edges)[opos];
                    long int r = RNG_INTEGER(0, no_of_nodes - 2);
                    VECTOR(edges)[ to_rewire ] = (r != nei ? r : no_of_nodes - 1);
                }
                to_rewire += RNG_GEOM(prob) + 1;
            }

        } else {
            IGRAPH_CHECK(igraph_i_rewire_edges_no_multiple(graph, prob, loops,
                         &edges));
        }
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(&newgraph, &edges, (igraph_integer_t) no_of_nodes,
                               igraph_is_directed(graph)));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &newgraph);
    IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
    IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 1);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = newgraph;

    return 0;
}

/**
 * \function igraph_rewire_directed_edges
 * \brief Rewire the chosen endpoint of directed edges
 *
 * This function rewires either the start or end of directed edges in a graph
 * with a constant probability. Correspondingly, either the in-degree sequence
 * or the out-degree sequence of the graph will be preserved.
 *
 * </para><para> Note that this function modifies the input \p graph,
 * call \ref igraph_copy() if you want to keep it.
 *
 * \param graph The input graph, this will be rewired, it can be
 *    directed or undirected. If it is directed, \ref igraph_rewire_edges()
 *    will be called.
 * \param prob The rewiring probability, a constant between zero and
 *    one (inclusive).
 * \param loops Boolean, whether loop edges are allowed in the new
 *    graph, or not.
 * \param mode The endpoints of directed edges to rewire. It is ignored for
 *    undirected graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          rewire the end of each directed edge
 *        \cli IGRAPH_IN
 *          rewire the start of each directed edge
 *        \cli IGRAPH_ALL
 *          rewire both endpoints of each edge
 *        \endclist
 * \return Error code.
 *
 * \sa \ref igraph_rewire_edges(), \ref igraph_rewire()
 *
 * Time complexity: O(|E|).
 */

int igraph_rewire_directed_edges(igraph_t *graph, igraph_real_t prob,
                                 igraph_bool_t loops, igraph_neimode_t mode) {

    if (prob < 0 || prob > 1) {
        IGRAPH_ERROR("Rewiring probability should be between zero and one",
                     IGRAPH_EINVAL);
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    if (prob == 0) {
        return IGRAPH_SUCCESS;
    }

    if (igraph_is_directed(graph) && mode != IGRAPH_ALL) {
        igraph_t newgraph;
        long int no_of_edges = igraph_ecount(graph);
        long int no_of_nodes = igraph_vcount(graph);
        long int to_rewire;
        long int offset;
        igraph_vector_t edges;

        IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * no_of_edges);

        switch (mode) {
        case IGRAPH_IN:
            offset = 0;
            break;
        case IGRAPH_OUT:
            offset = 1;
            break;
        case IGRAPH_ALL:
            break; /* suppress compiler warning */
        }

        IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));

        RNG_BEGIN();

        to_rewire = RNG_GEOM(prob);
        while (to_rewire < no_of_edges) {
            if (loops) {
                VECTOR(edges)[2 * to_rewire + offset] = RNG_INTEGER(0, no_of_nodes - 1);
            } else {
                long int nei = (long int) VECTOR(edges)[2 * to_rewire + (1 - offset)];
                long int r = RNG_INTEGER(0, no_of_nodes - 2);
                VECTOR(edges)[2 * to_rewire + offset] = (r != nei ? r : no_of_nodes - 1);
            }
            to_rewire += RNG_GEOM(prob) + 1;
        }

        RNG_END();

        IGRAPH_CHECK(igraph_create(&newgraph, &edges, (igraph_integer_t) no_of_nodes,
                                   igraph_is_directed(graph)));
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);

        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 1);
        IGRAPH_FINALLY_CLEAN(1);
        igraph_destroy(graph);
        *graph = newgraph;

    } else {
        IGRAPH_CHECK(igraph_rewire_edges(graph, prob, loops, /* multiple = */ 0));
    }

    return 0;
}

/**
 * \function igraph_watts_strogatz_game
 * \brief The Watts-Strogatz small-world model
 *
 * This function generates a graph according to the Watts-Strogatz
 * model of small-world networks. The graph is obtained by creating a
 * circular undirected lattice and then rewire the edges randomly with
 * a constant probability.
 *
 * </para><para>See also: Duncan J Watts and Steven H Strogatz:
 * Collective dynamics of <quote>small world</quote> networks, Nature
 * 393, 440-442, 1998.
 * \param graph The graph to initialize.
 * \param dim The dimension of the lattice.
 * \param size The size of the lattice along each dimension.
 * \param nei The size of the neighborhood for each vertex. This is
 *    the same as the \p nei argument of \ref
 *    igraph_connect_neighborhood().
 * \param p The rewiring probability. A real number between zero and
 *   one (inclusive).
 * \param loops Logical, whether to generate loop edges.
 * \param multiple Logical, whether to allow multiple edges in the
 *   generated graph.
 * \return Error code.
 *
 * \sa \ref igraph_lattice(), \ref igraph_connect_neighborhood() and
 * \ref igraph_rewire_edges() can be used if more flexibility is
 * needed, eg. a different type of lattice.
 *
 * Time complexity: O(|V|*d^o+|E|), |V| and |E| are the number of
 * vertices and edges, d is the average degree, o is the \p nei
 * argument.
 */

int igraph_watts_strogatz_game(igraph_t *graph, igraph_integer_t dim,
                               igraph_integer_t size, igraph_integer_t nei,
                               igraph_real_t p, igraph_bool_t loops,
                               igraph_bool_t multiple) {

    igraph_vector_t dimvector;
    long int i;

    if (dim < 1) {
        IGRAPH_ERROR("WS game: dimension should be at least one", IGRAPH_EINVAL);
    }
    if (size < 1) {
        IGRAPH_ERROR("WS game: lattice size should be at least one",
                     IGRAPH_EINVAL);
    }
    if (p < 0 || p > 1) {
        IGRAPH_ERROR("WS game: rewiring probability should be between 0 and 1",
                     IGRAPH_EINVAL);
    }

    /* Create the lattice first */

    IGRAPH_VECTOR_INIT_FINALLY(&dimvector, dim);
    for (i = 0; i < dim; i++) {
        VECTOR(dimvector)[i] = size;
    }

    IGRAPH_CHECK(igraph_lattice(graph, &dimvector, nei, IGRAPH_UNDIRECTED,
                                0 /* mutual */, 1 /* circular */));
    igraph_vector_destroy(&dimvector);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_destroy, graph);

    /* Rewire the edges then */

    IGRAPH_CHECK(igraph_rewire_edges(graph, p, loops, multiple));

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_lastcit_game
 * \brief Simulate citation network, based on time passed since the last citation.
 *
 * This is a quite special stochastic graph generator, it models an
 * evolving graph. In each time step a single vertex is added to the
 * network and it cites a number of other vertices (as specified by
 * the \p edges_per_step argument). The cited vertices are selected
 * based on the last time they were cited. Time is measured by the
 * addition of vertices and it is binned into \p pagebins bins.
 * So if the current time step is \c t and the last citation to a
 * given \c i vertex was made in time step \c t0, then \c
 * (t-t0)/binwidth is calculated where binwidth is \c nodes/pagebins+1,
 * in the last expression '/' denotes integer division, so the
 * fraction part is omitted.
 *
 * </para><para>
 * The \p preference argument specifies the preferences for the
 * citation lags, ie. its first elements contains the attractivity
 * of the very recently cited vertices, etc. The last element is
 * special, it contains the attractivity of the vertices which were
 * never cited. This element should be bigger than zero.
 *
 * </para><para>
 * Note that this function generates networks with multiple edges if
 * \p edges_per_step is bigger than one, call \ref igraph_simplify()
 * on the result to get rid of these edges.
 * \param graph Pointer to an uninitialized graph object, the result
 *     will be stored here.
 * \param node The number of vertices in the network.
 * \param edges_per_node The number of edges to add in each time
 *     step.
 * \param pagebins The number of age bins to use.
 * \param preference Pointer to an initialized vector of length
 *     \c pagebins+1. This contains the `attractivity' of the various
 *     age bins, the last element is the attractivity of the vertices
 *     which were never cited, and it should be greater than zero.
 *     It is a good idea to have all positive values in this vector.
 * \param directed Logical constant, whether to create directed
 *      networks.
 * \return Error code.
 *
 * \sa \ref igraph_barabasi_aging_game().
 *
 * Time complexity: O(|V|*a+|E|*log|V|), |V| is the number of vertices,
 * |E| is the total number of edges, a is the \p pagebins parameter.
 */

int igraph_lastcit_game(igraph_t *graph,
                        igraph_integer_t nodes, igraph_integer_t edges_per_node,
                        igraph_integer_t pagebins,
                        const igraph_vector_t *preference,
                        igraph_bool_t directed) {

    long int no_of_nodes = nodes;
    igraph_psumtree_t sumtree;
    igraph_vector_t edges;
    long int i, j, k;
    long int *lastcit;
    long int *index;
    long int agebins = pagebins;
    long int binwidth = no_of_nodes / agebins + 1;

    if (agebins != igraph_vector_size(preference) - 1) {
        IGRAPH_ERROR("`preference' vector should be of length `agebins' plus one",
                     IGRAPH_EINVAL);
    }
    if (agebins <= 1 ) {
        IGRAPH_ERROR("at least two age bins are need for lastcit game",
                     IGRAPH_EINVAL);
    }
    if (VECTOR(*preference)[agebins] <= 0) {
        IGRAPH_ERROR("the last element of the `preference' vector needs to be positive",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    lastcit = igraph_Calloc(no_of_nodes, long int);
    if (!lastcit) {
        IGRAPH_ERROR("lastcit game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, lastcit);

    index = igraph_Calloc(no_of_nodes + 1, long int);
    if (!index) {
        IGRAPH_ERROR("lastcit game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, index);

    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes * edges_per_node));

    /* The first node */
    igraph_psumtree_update(&sumtree, 0, VECTOR(*preference)[agebins]);
    index[0] = 0;
    index[1] = 0;

    RNG_BEGIN();

    for (i = 1; i < no_of_nodes; i++) {

        /* Add new edges */
        for (j = 0; j < edges_per_node; j++) {
            long int to;
            igraph_real_t sum = igraph_psumtree_sum(&sumtree);
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, to);
            lastcit[to] = i + 1;
            igraph_psumtree_update(&sumtree, to, VECTOR(*preference)[0]);
        }

        /* Add the node itself */
        igraph_psumtree_update(&sumtree, i, VECTOR(*preference)[agebins]);
        index[i + 1] = index[i] + edges_per_node;

        /* Update the preference of some vertices if they got to another bin.
           We need to know the citations of some older vertices, this is in the index. */
        for (k = 1; i - binwidth * k >= 1; k++) {
            long int shnode = i - binwidth * k;
            long int m = index[shnode], n = index[shnode + 1];
            for (j = 2 * m; j < 2 * n; j += 2) {
                long int cnode = (long int) VECTOR(edges)[j + 1];
                if (lastcit[cnode] == shnode + 1) {
                    igraph_psumtree_update(&sumtree, cnode, VECTOR(*preference)[k]);
                }
            }
        }

    }

    RNG_END();

    igraph_psumtree_destroy(&sumtree);
    igraph_free(index);
    igraph_free(lastcit);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_cited_type_game
 * \brief Simulate a citation based on vertex types.
 *
 * Function to create a network based on some vertex categories. This
 * function creates a citation network, in each step a single vertex
 * and \p edges_per_step citating edges are added, nodes with
 * different categories (may) have different probabilities to get
 * cited, as given by the \p pref vector.
 *
 * </para><para>
 * Note that this function might generate networks with multiple edges
 * if \p edges_per_step is greater than one. You might want to call
 * \ref igraph_simplify() on the result to remove multiple edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the network.
 * \param types Numeric vector giving the categories of the vertices,
 *     so it should contain \p nodes non-negative integer
 *     numbers. Types are numbered from zero.
 * \param pref The attractivity of the different vertex categories in
 *     a vector. Its length should be the maximum element in \p types
 *     plus one (types are numbered from zero).
 * \param edges_per_step Integer constant, the number of edges to add
 *     in each time step.
 * \param directed Logical constant, whether to create a directed
 *     network.
 * \return Error code.
 *
 * \sa \ref igraph_citing_cited_type_game() for a bit more general
 * game.
 *
 * Time complexity: O((|V|+|E|)log|V|), |V| and |E| are number of
 * vertices and edges, respectively.
 */

int igraph_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
                           const igraph_vector_t *types,
                           const igraph_vector_t *pref,
                           igraph_integer_t edges_per_step,
                           igraph_bool_t directed) {

    igraph_vector_t edges;
    igraph_vector_t cumsum;
    igraph_real_t sum;
    long int i, j, nnval, type;

    if (igraph_vector_size(types) != nodes) {
        IGRAPH_ERROR("Invalid size of types", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    /* return an empty graph is nodes is zero */
    if (nodes == 0) {
        igraph_create(graph, &edges, nodes, directed);
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&cumsum, 2);
    IGRAPH_CHECK(igraph_vector_reserve(&cumsum, nodes + 1));
    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes * edges_per_step));

    /* first node */
    VECTOR(cumsum)[0] = 0;
    type = (long int) VECTOR(*types)[0];
    if (type >= igraph_vector_size(pref)) {
        IGRAPH_ERROR("pref is too short for the given types", IGRAPH_EINVAL);
    }
    nnval = VECTOR(*pref)[type];
    if (nnval < 0) {
        IGRAPH_ERROR("pref contains negative entries", IGRAPH_EINVAL);
    }
    sum = VECTOR(cumsum)[1] = nnval;

    RNG_BEGIN();

    for (i = 1; i < nodes; i++) {
        for (j = 0; j < edges_per_step; j++) {
            long int to;
            if (sum > 0) {
                igraph_vector_binsearch(&cumsum, RNG_UNIF(0, sum), &to);
            } else {
                to = i + 1;
            }
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, to - 1);
        }
        type = (long int) VECTOR(*types)[i];
        if (type >= igraph_vector_size(pref)) {
            IGRAPH_ERROR("pref is too short for the given types", IGRAPH_EINVAL);
        }
        nnval = VECTOR(*pref)[type];
        if (nnval < 0) {
            IGRAPH_ERROR("pref contains negative entries", IGRAPH_EINVAL);
        }
        sum += nnval;
        igraph_vector_push_back(&cumsum, sum);
    }

    RNG_END();

    igraph_vector_destroy(&cumsum);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static void igraph_i_citing_cited_type_game_free(igraph_i_citing_cited_type_game_struct_t *s) {
    long int i;
    if (!s->sumtrees) {
        return;
    }
    for (i = 0; i < s->no; i++) {
        igraph_psumtree_destroy(&s->sumtrees[i]);
    }
}

/**
 * \function igraph_citing_cited_type_game
 * \brief Simulate a citation network based on vertex types.
 *
 * This game is similar to \ref igraph_cited_type_game() but here the
 * category of the citing vertex is also considered.
 *
 * </para><para>
 * An evolving citation network is modeled here, a single vertex and
 * its \p edges_per_step citation are added in each time step. The
 * odds the a given vertex is cited by the new vertex depends on the
 * category of both the citing and the cited vertex and is given in
 * the \p pref matrix. The categories of the citing vertex correspond
 * to the rows, the categories of the cited vertex to the columns of
 * this matrix. Ie. the element in row \c i and column \c j gives the
 * probability that a \c j vertex is cited, if the category of the
 * citing vertex is \c i.
 *
 * </para><para>
 * Note that this function might generate networks with multiple edges
 * if \p edges_per_step is greater than one. You might want to call
 * \ref igraph_simplify() on the result to remove multiple edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the network.
 * \param types A numeric matrix of length \p nodes, containing the
 *    categories of the vertices. The categories are numbered from
 *    zero.
 * \param pref The preference matrix, a square matrix is required,
 *     both the number of rows and columns should be the maximum
 *     element in \p types plus one (types are numbered from zero).
 * \param directed Logical constant, whether to create a directed
 *     network.
 * \return Error code.
 *
 * Time complexity: O((|V|+|E|)log|V|), |V| and |E| are number of
 * vertices and edges, respectively.
 */

int igraph_citing_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
                                  const igraph_vector_t *types,
                                  const igraph_matrix_t *pref,
                                  igraph_integer_t edges_per_step,
                                  igraph_bool_t directed) {

    igraph_vector_t edges;
    igraph_i_citing_cited_type_game_struct_t str = { 0, 0 };
    igraph_psumtree_t *sumtrees;
    igraph_vector_t sums;
    long int nocats;
    long int i, j;

    if (igraph_vector_size(types) != nodes) {
        IGRAPH_ERROR("Invalid size of types", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    /* return an empty graph is nodes is zero */
    if (nodes == 0) {
        igraph_create(graph, &edges, nodes, directed);
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(2); /* str and edges */
        return 0;
    }

    nocats = igraph_matrix_ncol(pref);
    str.sumtrees = sumtrees = igraph_Calloc(nocats, igraph_psumtree_t);
    if (!sumtrees) {
        IGRAPH_ERROR("Citing-cited type game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_i_citing_cited_type_game_free, &str);

    for (i = 0; i < nocats; i++) {
        IGRAPH_CHECK(igraph_psumtree_init(&sumtrees[i], nodes));
        str.no++;
    }
    IGRAPH_VECTOR_INIT_FINALLY(&sums, nocats);

    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes * edges_per_step));

    /* First node */
    for (i = 0; i < nocats; i++) {
        long int type = (long int) VECTOR(*types)[0];
        if ( MATRIX(*pref, i, type) < 0) {
            IGRAPH_ERROR("pref contains negative entries", IGRAPH_EINVAL);
        }
        igraph_psumtree_update(&sumtrees[i], 0, MATRIX(*pref, i, type));
        VECTOR(sums)[i] = MATRIX(*pref, i, type);
    }

    RNG_BEGIN();

    for (i = 1; i < nodes; i++) {
        long int type = (long int) VECTOR(*types)[i];
        igraph_real_t sum = VECTOR(sums)[type];
        for (j = 0; j < edges_per_step; j++) {
            long int to;
            igraph_psumtree_search(&sumtrees[type], &to, RNG_UNIF(0, sum));
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, to);
        }

        /* add i */
        for (j = 0; j < nocats; j++) {
            if ( MATRIX(*pref, j, type) < 0) {
                IGRAPH_ERROR("pref contains negative entries", IGRAPH_EINVAL);
            }
            igraph_psumtree_update(&sumtrees[j], i, MATRIX(*pref, j,  type));
            VECTOR(sums)[j] += MATRIX(*pref, j, type);
        }
    }

    RNG_END();

    igraph_i_citing_cited_type_game_free(&str);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_create(graph, &edges, nodes, directed);
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}



/**
 * \ingroup generators
 * \function igraph_simple_interconnected_islands_game
 * \brief Generates a random graph made of several interconnected islands, each island being a random graph.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param islands_n The number of islands in the graph.
 * \param islands_size The size of islands in the graph.
 * \param islands_pin The probability to create each possible edge into each island .
 * \param n_inter The number of edges to create between two islands .

 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid parameter
 *         \c IGRAPH_ENOMEM: there is not enough
 *         memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 */
int igraph_simple_interconnected_islands_game(
        igraph_t *graph,
        igraph_integer_t islands_n,
        igraph_integer_t islands_size,
        igraph_real_t islands_pin,
        igraph_integer_t n_inter) {


    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int nbNodes;
    double maxpossibleedgesPerIsland;
    double maxedgesPerIsland;
    int nbEdgesInterIslands;
    double maxedges;
    int startIsland = 0;
    int endIsland = 0;
    int i, j, is;
    double myrand, last;

    if (islands_n < 0) {
        IGRAPH_ERROR("Invalid number of islands", IGRAPH_EINVAL);
    }
    if (islands_size < 0) {
        IGRAPH_ERROR("Invalid size for islands", IGRAPH_EINVAL);
    }
    if (islands_pin < 0 || islands_pin > 1) {
        IGRAPH_ERROR("Invalid probability for islands", IGRAPH_EINVAL);
    }
    if ( (n_inter < 0) || (n_inter > islands_size) ) {
        IGRAPH_ERROR("Invalid number of inter-islands links", IGRAPH_EINVAL);
    }

    /* how much memory ? */
    nbNodes = islands_n * islands_size;
    maxpossibleedgesPerIsland = ((double)islands_size * ((double)islands_size - (double)1)) / (double)2;
    maxedgesPerIsland = islands_pin * maxpossibleedgesPerIsland;
    nbEdgesInterIslands = n_inter * (islands_n * (islands_n - 1)) / 2;
    maxedges = maxedgesPerIsland * islands_n + nbEdgesInterIslands;    

    /* reserve enough space for all the edges */
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, (long int) maxedges));

    RNG_BEGIN();

    /* first create all the islands */
    for (is = 1; is <= islands_n; is++) { /* for each island */

        /* index for start and end of nodes in this island */
        startIsland = islands_size * (is - 1);
        endIsland = startIsland + islands_size - 1;

        /* create the random numbers to be used (into s) */
        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&s, (long int) maxedgesPerIsland));

        last = RNG_GEOM(islands_pin);
        while (last < maxpossibleedgesPerIsland) { /* maxedgesPerIsland */
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            myrand = RNG_GEOM(islands_pin);
            last += myrand; /* RNG_GEOM(islands_pin); */
            last += 1;
        }



        /* change this to edges ! */
        for (i = 0; i < igraph_vector_size(&s); i++) {
            long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
            long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2);
            to += startIsland;
            from += startIsland;

            igraph_vector_push_back(&edges, from);
            igraph_vector_push_back(&edges, to);
        }

        /* clear the memory used for random number for this island */
        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);


        /* create the links with other islands */
        for (i = is + 1; i <= islands_n; i++) { /* for each other island (not the previous ones) */

            for (j = 0; j < n_inter; j++) { /* for each link between islands */
                long int from = (long int) RNG_UNIF(startIsland, endIsland);
                long int to = (long int) RNG_UNIF((i - 1) * islands_size, i * islands_size);

                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }

        }
    }

    RNG_END();

    /* actually fill the graph object */
    IGRAPH_CHECK(igraph_create(graph, &edges, nbNodes, 0));

    /* clean remaining things */
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup generators
 * \function igraph_static_fitness_game
 * \brief Generates a non-growing random graph with edge probabilities
 *        proportional to node fitness scores.
 *
 * This game generates a directed or undirected random graph where the
 * probability of an edge between vertices i and j depends on the fitness
 * scores of the two vertices involved. For undirected graphs, each vertex
 * has a single fitness score. For directed graphs, each vertex has an out-
 * and an in-fitness, and the probability of an edge from i to j depends on
 * the out-fitness of vertex i and the in-fitness of vertex j.
 *
 * </para><para>
 * The generation process goes as follows. We start from N disconnected nodes
 * (where N is given by the length of the fitness vector). Then we randomly
 * select two vertices i and j, with probabilities proportional to their
 * fitnesses. (When the generated graph is directed, i is selected according to
 * the out-fitnesses and j is selected according to the in-fitnesses). If the
 * vertices are not connected yet (or if multiple edges are allowed), we
 * connect them; otherwise we select a new pair. This is repeated until the
 * desired number of links are created.
 *
 * </para><para>
 * It can be shown that the \em expected degree of each vertex will be
 * proportional to its fitness, although the actual, observed degree will not
 * be. If you need to generate a graph with an exact degree sequence, consider
 * \ref igraph_degree_sequence_game instead.
 *
 * </para><para>
 * This model is commonly used to generate static scale-free networks. To
 * achieve this, you have to draw the fitness scores from the desired power-law
 * distribution. Alternatively, you may use \ref igraph_static_power_law_game
 * which generates the fitnesses for you with a given exponent.
 *
 * </para><para>
 * Reference: Goh K-I, Kahng B, Kim D: Universal behaviour of load distribution
 * in scale-free networks. Phys Rev Lett 87(27):278701, 2001.
 *
 * \param graph        Pointer to an uninitialized graph object.
 * \param fitness_out  A numeric vector containing the fitness of each vertex.
 *                     For directed graphs, this specifies the out-fitness
 *                     of each vertex.
 * \param fitness_in   If \c NULL, the generated graph will be undirected.
 *                     If not \c NULL, this argument specifies the in-fitness
 *                     of each vertex.
 * \param no_of_edges  The number of edges in the generated graph.
 * \param loops        Whether to allow loop edges in the generated graph.
 * \param multiple     Whether to allow multiple edges in the generated graph.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid parameter
 *         \c IGRAPH_ENOMEM: there is not enough
 *         memory for the operation.
 *
 * Time complexity: O(|V| + |E| log |E|).
 */
int igraph_static_fitness_game(igraph_t *graph, igraph_integer_t no_of_edges,
                               igraph_vector_t* fitness_out, igraph_vector_t* fitness_in,
                               igraph_bool_t loops, igraph_bool_t multiple) {
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t no_of_nodes;
    igraph_integer_t outnodes, innodes, nodes;
    igraph_vector_t cum_fitness_in, cum_fitness_out;
    igraph_vector_t *p_cum_fitness_in, *p_cum_fitness_out;
    igraph_real_t x, max_in, max_out;
    igraph_real_t max_no_of_edges;
    igraph_bool_t is_directed = (fitness_in != 0);
    float num_steps;
    igraph_integer_t step_counter = 0;
    long int i, from, to, pos;

    if (fitness_out == 0) {
        IGRAPH_ERROR("fitness_out must not be null", IGRAPH_EINVAL);
    }

    if (no_of_edges < 0) {
        IGRAPH_ERROR("Invalid number of edges", IGRAPH_EINVAL);
    }

    no_of_nodes = (int) igraph_vector_size(fitness_out);
    if (no_of_nodes == 0) {
        IGRAPH_CHECK(igraph_empty(graph, 0, is_directed));
        return IGRAPH_SUCCESS;
    }

    if (is_directed && igraph_vector_size(fitness_in) != no_of_nodes) {
        IGRAPH_ERROR("fitness_in must have the same size as fitness_out", IGRAPH_EINVAL);
    }

    /* Sanity checks for the fitnesses */
    if (igraph_vector_min(fitness_out) < 0) {
        IGRAPH_ERROR("Fitness scores must be non-negative", IGRAPH_EINVAL);
    }
    if (fitness_in != 0 && igraph_vector_min(fitness_in) < 0) {
        IGRAPH_ERROR("Fitness scores must be non-negative", IGRAPH_EINVAL);
    }

    /* Avoid getting into an infinite loop when too many edges are requested */
    if (!multiple) {
        if (is_directed) {
            outnodes = innodes = nodes = 0;
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*fitness_out)[i] != 0) {
                    outnodes++;
                }
                if (VECTOR(*fitness_in)[i] != 0) {
                    innodes++;
                }
                if (VECTOR(*fitness_out)[i] != 0 && VECTOR(*fitness_in)[i] != 0) {
                    nodes++;
                }
            }
            max_no_of_edges = ((igraph_real_t) outnodes) * innodes - (loops ? 0 : nodes);
        } else {
            nodes = 0;
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*fitness_out)[i] != 0) {
                    nodes++;
                }
            }
            max_no_of_edges = loops
                              ? nodes * ((igraph_real_t)nodes + 1) / 2
                              : nodes * ((igraph_real_t)nodes - 1) / 2;
        }
        if (no_of_edges > max_no_of_edges) {
            IGRAPH_ERROR("Too many edges requested", IGRAPH_EINVAL);
        }
    }

    /* Calculate the cumulative fitness scores */
    IGRAPH_VECTOR_INIT_FINALLY(&cum_fitness_out, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_cumsum(&cum_fitness_out, fitness_out));
    max_out = igraph_vector_tail(&cum_fitness_out);
    p_cum_fitness_out = &cum_fitness_out;
    if (is_directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&cum_fitness_in, no_of_nodes);
        IGRAPH_CHECK(igraph_vector_cumsum(&cum_fitness_in, fitness_in));
        max_in = igraph_vector_tail(&cum_fitness_in);
        p_cum_fitness_in = &cum_fitness_in;
    } else {
        max_in = max_out;
        p_cum_fitness_in = &cum_fitness_out;
    }

    RNG_BEGIN();
    num_steps = no_of_edges;
    if (multiple) {
        /* Generating when multiple edges are allowed */

        IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&edges, 2 * no_of_edges));

        while (no_of_edges > 0) {
            /* Report progress after every 10000 edges */
            if ((step_counter++) % 10000 == 0) {
                IGRAPH_PROGRESS("Static fitness game", 100.0 * (1 - no_of_edges / num_steps), NULL);
                IGRAPH_ALLOW_INTERRUPTION();
            }

            x = RNG_UNIF(0, max_out);
            igraph_vector_binsearch(p_cum_fitness_out, x, &from);
            x = RNG_UNIF(0, max_in);
            igraph_vector_binsearch(p_cum_fitness_in, x, &to);

            /* Skip if loop edge and loops = false */
            if (!loops && from == to) {
                continue;
            }

            igraph_vector_push_back(&edges, from);
            igraph_vector_push_back(&edges, to);

            no_of_edges--;
        }

        /* Create the graph */
        IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, is_directed));

        /* Clear the edge list */
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Multiple edges are disallowed */
        igraph_adjlist_t al;
        igraph_vector_int_t* neis;

        IGRAPH_CHECK(igraph_adjlist_init_empty(&al, no_of_nodes));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
        while (no_of_edges > 0) {
            /* Report progress after every 10000 edges */
            if ((step_counter++) % 10000 == 0) {
                IGRAPH_PROGRESS("Static fitness game", 100.0 * (1 - no_of_edges / num_steps), NULL);
                IGRAPH_ALLOW_INTERRUPTION();
            }

            x = RNG_UNIF(0, max_out);
            igraph_vector_binsearch(p_cum_fitness_out, x, &from);
            x = RNG_UNIF(0, max_in);
            igraph_vector_binsearch(p_cum_fitness_in, x, &to);

            /* Skip if loop edge and loops = false */
            if (!loops && from == to) {
                continue;
            }

            /* For undirected graphs, ensure that from < to */
            if (!is_directed && from > to) {
                pos = from; from = to; to = pos;
            }

            /* Is there already an edge? If so, try again */
            neis = igraph_adjlist_get(&al, from);
            if (igraph_vector_int_binsearch(neis, to, &pos)) {
                continue;
            }

            /* Insert the edge */
            IGRAPH_CHECK(igraph_vector_int_insert(neis, pos, to));

            no_of_edges--;
        }

        /* Create the graph. We cannot use IGRAPH_ALL here for undirected graphs
         * because we did not add edges in both directions in the adjacency list.
         * We will use igraph_to_undirected in an extra step. */
        IGRAPH_CHECK(igraph_adjlist(graph, &al, IGRAPH_OUT, 1));
        if (!is_directed) {
            IGRAPH_CHECK(igraph_to_undirected(graph, IGRAPH_TO_UNDIRECTED_EACH, 0));
        }

        /* Clear the adjacency list */
        igraph_adjlist_destroy(&al);
        IGRAPH_FINALLY_CLEAN(1);
    }
    RNG_END();

    IGRAPH_PROGRESS("Static fitness game", 100.0, NULL);

    /* Cleanup before we create the graph */
    if (is_directed) {
        igraph_vector_destroy(&cum_fitness_in);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&cum_fitness_out);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup generators
 * \function igraph_static_power_law_game
 * \brief Generates a non-growing random graph with expected power-law degree distributions.
 *
 * This game generates a directed or undirected random graph where the
 * degrees of vertices follow power-law distributions with prescribed
 * exponents. For directed graphs, the exponents of the in- and out-degree
 * distributions may be specified separately.
 *
 * </para><para>
 * The game simply uses \ref igraph_static_fitness_game with appropriately
 * constructed fitness vectors. In particular, the fitness of vertex i
 * is i<superscript>-alpha</superscript>, where alpha = 1/(gamma-1)
 * and gamma is the exponent given in the arguments.
 *
 * </para><para>
 * To remove correlations between in- and out-degrees in case of directed
 * graphs, the in-fitness vector will be shuffled after it has been set up
 * and before \ref igraph_static_fitness_game is called.
 *
 * </para><para>
 * Note that significant finite size effects may be observed for exponents
 * smaller than 3 in the original formulation of the game. This function
 * provides an argument that lets you remove the finite size effects by
 * assuming that the fitness of vertex i is
 * (i+i0-1)<superscript>-alpha</superscript>,
 * where i0 is a constant chosen appropriately to ensure that the maximum
 * degree is less than the square root of the number of edges times the
 * average degree; see the paper of Chung and Lu, and Cho et al for more
 * details.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Goh K-I, Kahng B, Kim D: Universal behaviour of load distribution
 * in scale-free networks. Phys Rev Lett 87(27):278701, 2001.
 *
 * </para><para>
 * Chung F and Lu L: Connected components in a random graph with given
 * degree sequences. Annals of Combinatorics 6, 125-145, 2002.
 *
 * </para><para>
 * Cho YS, Kim JS, Park J, Kahng B, Kim D: Percolation transitions in
 * scale-free networks under the Achlioptas process. Phys Rev Lett
 * 103:135702, 2009.
 *
 * \param graph        Pointer to an uninitialized graph object.
 * \param no_of_nodes  The number of nodes in the generated graph.
 * \param no_of_edges  The number of edges in the generated graph.
 * \param exponent_out The power law exponent of the degree distribution.
 *                     For directed graphs, this specifies the exponent of the
 *                     out-degree distribution. It must be greater than or
 *                     equal to 2. If you pass \c IGRAPH_INFINITY here, you
 *                     will get back an Erdos-Renyi random network.
 * \param exponent_in  If negative, the generated graph will be undirected.
 *                     If greater than or equal to 2, this argument specifies
 *                     the exponent of the in-degree distribution. If
 *                     non-negative but less than 2, an error will be
 *                     generated.
 * \param loops        Whether to allow loop edges in the generated graph.
 * \param multiple     Whether to allow multiple edges in the generated graph.
 * \param finite_size_correction  Whether to use the proposed finite size
 *                     correction of Cho et al.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid parameter
 *         \c IGRAPH_ENOMEM: there is not enough
 *         memory for the operation.
 *
 * Time complexity: O(|V| + |E| log |E|).
 */
int igraph_static_power_law_game(igraph_t *graph,
                                 igraph_integer_t no_of_nodes, igraph_integer_t no_of_edges,
                                 igraph_real_t exponent_out, igraph_real_t exponent_in,
                                 igraph_bool_t loops, igraph_bool_t multiple,
                                 igraph_bool_t finite_size_correction) {

    igraph_vector_t fitness_out, fitness_in;
    igraph_real_t alpha_out = 0.0, alpha_in = 0.0;
    long int i;
    igraph_real_t j;

    if (no_of_nodes < 0) {
        IGRAPH_ERROR("Invalid number of nodes", IGRAPH_EINVAL);
    }

    /* Calculate alpha_out */
    if (exponent_out < 2) {
        IGRAPH_ERROR("out-degree exponent must be >= 2", IGRAPH_EINVAL);
    } else if (igraph_finite(exponent_out)) {
        alpha_out = -1.0 / (exponent_out - 1);
    } else {
        alpha_out = 0.0;
    }

    /* Construct the out-fitnesses */
    IGRAPH_VECTOR_INIT_FINALLY(&fitness_out, no_of_nodes);
    j = no_of_nodes;
    if (finite_size_correction && alpha_out < -0.5) {
        /* See the Cho et al paper, first page first column + footnote 7 */
        j += pow(no_of_nodes, 1 + 0.5 / alpha_out) *
             pow(10 * sqrt(2) * (1 + alpha_out), -1.0 / alpha_out) - 1;
    }
    if (j < no_of_nodes) {
        j = no_of_nodes;
    }
    for (i = 0; i < no_of_nodes; i++, j--) {
        VECTOR(fitness_out)[i] = pow(j, alpha_out);
    }

    if (exponent_in >= 0) {
        if (exponent_in < 2) {
            IGRAPH_ERROR("in-degree exponent must be >= 2; use negative numbers "
                         "for undirected graphs", IGRAPH_EINVAL);
        } else if (igraph_finite(exponent_in)) {
            alpha_in = -1.0 / (exponent_in - 1);
        } else {
            alpha_in = 0.0;
        }

        IGRAPH_VECTOR_INIT_FINALLY(&fitness_in, no_of_nodes);
        j = no_of_nodes;
        if (finite_size_correction && alpha_in < -0.5) {
            /* See the Cho et al paper, first page first column + footnote 7 */
            j += pow(no_of_nodes, 1 + 0.5 / alpha_in) *
                 pow(10 * sqrt(2) * (1 + alpha_in), -1.0 / alpha_in) - 1;
        }
        if (j < no_of_nodes) {
            j = no_of_nodes;
        }
        for (i = 0; i < no_of_nodes; i++, j--) {
            VECTOR(fitness_in)[i] = pow(j, alpha_in);
        }
        IGRAPH_CHECK(igraph_vector_shuffle(&fitness_in));

        IGRAPH_CHECK(igraph_static_fitness_game(graph, no_of_edges,
                                                &fitness_out, &fitness_in, loops, multiple));

        igraph_vector_destroy(&fitness_in);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        IGRAPH_CHECK(igraph_static_fitness_game(graph, no_of_edges,
                                                &fitness_out, 0, loops, multiple));
    }

    igraph_vector_destroy(&fitness_out);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup generators
 * \function igraph_k_regular_game
 * \brief Generates a random graph where each vertex has the same degree.
 *
 * This game generates a directed or undirected random graph where the
 * degrees of vertices are equal to a predefined constant k. For undirected
 * graphs, at least one of k and the number of vertices must be even.
 *
 * </para><para>
 * The game simply uses \ref igraph_degree_sequence_game with appropriately
 * constructed degree sequences.
 *
 * \param graph        Pointer to an uninitialized graph object.
 * \param no_of_nodes  The number of nodes in the generated graph.
 * \param k            The degree of each vertex in an undirected graph, or
 *                     the out-degree and in-degree of each vertex in a
 *                     directed graph.
 * \param directed     Whether the generated graph will be directed.
 * \param multiple     Whether to allow multiple edges in the generated graph.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid parameter; e.g., negative number of nodes,
 *                           or odd number of nodes and odd k for undirected
 *                           graphs.
 *         \c IGRAPH_ENOMEM: there is not enough memory for the operation.
 *
 * Time complexity: O(|V|+|E|) if \c multiple is true, otherwise not known.
 */
int igraph_k_regular_game(igraph_t *graph,
                          igraph_integer_t no_of_nodes, igraph_integer_t k,
                          igraph_bool_t directed, igraph_bool_t multiple) {
    igraph_vector_t degseq;
    igraph_degseq_t mode = multiple ? IGRAPH_DEGSEQ_SIMPLE : IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE;

    /* Note to self: we are not using IGRAPH_DEGSEQ_VL when multiple = false
     * because the VL method is not really good at generating k-regular graphs.
     * Actually, that's why we have added SIMPLE_NO_MULTIPLE. */

    if (no_of_nodes < 0) {
        IGRAPH_ERROR("number of nodes must be non-negative", IGRAPH_EINVAL);
    }
    if (k < 0) {
        IGRAPH_ERROR("degree must be non-negative", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&degseq, no_of_nodes);
    igraph_vector_fill(&degseq, k);
    IGRAPH_CHECK(igraph_degree_sequence_game(graph, &degseq, directed ? &degseq : 0, mode));

    igraph_vector_destroy(&degseq);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_correlated_game
 * Generate pairs of correlated random graphs
 *
 * Sample a new graph by perturbing the adjacency matrix of a
 * given graph and shuffling its vertices.
 *
 * \param old_graph The original graph.
 * \param new_graph The new graph will be stored here.
 * \param corr A scalar in the unit interval, the target Pearson
 *        correlation between the adjacency matrices of the original the
 *        generated graph (the adjacency matrix being used as a vector).
 * \param p A numeric scalar, the probability of an edge between two
 *        vertices, it must in the open (0,1) interval.
 * \param permutation A permutation to apply to the vertices of the
 *        generated graph. It can also be a null pointer, in which case
 *        the vertices will not be permuted.
 * \return Error code
 *
 * \sa \ref igraph_correlated_pair_game() for generating a pair
 * of correlated random graphs in one go.
 */

int igraph_correlated_game(const igraph_t *old_graph, igraph_t *new_graph,
                           igraph_real_t corr, igraph_real_t p,
                           const igraph_vector_t *permutation) {

    int no_of_nodes = igraph_vcount(old_graph);
    int no_of_edges = igraph_ecount(old_graph);
    igraph_bool_t directed = igraph_is_directed(old_graph);
    igraph_real_t no_of_all = directed ? no_of_nodes * (no_of_nodes - 1) :
                              no_of_nodes * (no_of_nodes - 1) / 2;
    igraph_real_t no_of_missing = no_of_all - no_of_edges;
    igraph_real_t q = p + corr * (1 - p);
    igraph_real_t p_del = 1 - q;
    igraph_real_t p_add = ((1 - q) * (p / (1 - p)));
    igraph_vector_t add, delete, edges, newedges;
    igraph_real_t last;
    int p_e = 0, p_a = 0, p_d = 0, no_add, no_del;
    igraph_real_t inf = IGRAPH_INFINITY;
    igraph_real_t next_e, next_a, next_d;
    int i;

    if (corr < -1 || corr > 1) {
        IGRAPH_ERROR("Correlation must be in [-1,1] in correlated "
                     "Erdos-Renyi game", IGRAPH_EINVAL);
    }
    if (p <= 0 || p >= 1) {
        IGRAPH_ERROR("Edge probability must be in (0,1) in correlated "
                     "Erdos-Renyi game", IGRAPH_EINVAL);
    }
    if (permutation) {
        if (igraph_vector_size(permutation) != no_of_nodes) {
            IGRAPH_ERROR("Invalid permutation length in correlated Erdos-Renyi game",
                         IGRAPH_EINVAL);
        }
    }

    /* Special cases */

    if (corr == 0) {
        return igraph_erdos_renyi_game(new_graph, IGRAPH_ERDOS_RENYI_GNP,
                                       no_of_nodes, p, directed,
                                       IGRAPH_NO_LOOPS);
    }
    if (corr == 1) {
        /* We don't copy, because we don't need the attributes.... */
        IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);
        IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));
        if (permutation) {
            int newec = igraph_vector_size(&edges);
            for (i = 0; i < newec; i++) {
                int tmp = VECTOR(edges)[i];
                VECTOR(edges)[i] = VECTOR(*permutation)[tmp];
            }
        }
        IGRAPH_CHECK(igraph_create(new_graph, &edges, no_of_nodes, directed));
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&newedges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&add, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&delete, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));

    RNG_BEGIN();

    if (p_del > 0) {
        last = RNG_GEOM(p_del);
        while (last < no_of_edges) {
            IGRAPH_CHECK(igraph_vector_push_back(&delete, last));
            last += RNG_GEOM(p_del);
            last += 1;
        }
    }
    no_del = igraph_vector_size(&delete);

    if (p_add > 0) {
        last = RNG_GEOM(p_add);
        while (last < no_of_missing) {
            IGRAPH_CHECK(igraph_vector_push_back(&add, last));
            last += RNG_GEOM(p_add);
            last += 1;
        }
    }
    no_add = igraph_vector_size(&add);

    RNG_END();

    IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));

    /* Now we are merging the original edges, the edges that are removed,
       and the new edges. We have the following pointers:
       - p_a: the next edge to add
       - p_d: the next edge to delete
       - p_e: the next original edge
       - next_e: the code of the next edge in 'edges'
       - next_a: the code of the next edge to add
       - next_d: the code of the next edge to delete */

#define D_CODE(f,t) (((t)==no_of_nodes-1 ? f : t) * no_of_nodes + (f))
#define U_CODE(f,t) ((t) * ((t)-1) / 2 + (f))
#define CODE(f,t) (directed ? D_CODE(f,t) : U_CODE(f,t))
#define CODEE() (CODE(VECTOR(edges)[2*p_e], VECTOR(edges)[2*p_e+1]))

    /* First we (re)code the edges to delete */

    for (i = 0; i < no_del; i++) {
        int td = VECTOR(delete)[i];
        int from = VECTOR(edges)[2 * td];
        int to = VECTOR(edges)[2 * td + 1];
        VECTOR(delete)[i] = CODE(from, to);
    }

    IGRAPH_CHECK(igraph_vector_reserve(&newedges,
                                       (no_of_edges - no_del + no_add) * 2));

    /* Now we can do the merge. Additional edges are tricky, because
       the code must be shifted by the edges in the original graph. */

#define UPD_E()                             \
    { if (p_e < no_of_edges) { next_e=CODEE(); } else { next_e = inf; } }
#define UPD_A()                             \
{ if (p_a < no_add) { \
            next_a = VECTOR(add)[p_a] + p_e; } else { next_a = inf; } }
#define UPD_D()                             \
{ if (p_d < no_del) { \
            next_d = VECTOR(delete)[p_d]; } else { next_d = inf; } }

    UPD_E(); UPD_A(); UPD_D();

    while (next_e != inf || next_a != inf || next_d != inf) {
        if (next_e <= next_a && next_e < next_d) {

            /* keep an edge */
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, VECTOR(edges)[2 * p_e]));
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, VECTOR(edges)[2 * p_e + 1]));
            p_e ++; UPD_E(); UPD_A()

        } else if (next_e <= next_a && next_e == next_d) {

            /* delete an edge */
            p_e ++; UPD_E(); UPD_A();
            p_d++; UPD_D();

        } else {

            /* add an edge */
            int to, from;
            if (directed) {
                to = (int) floor(next_a / no_of_nodes);
                from = (int) (next_a - ((igraph_real_t)to) * no_of_nodes);
                if (from == to) {
                    to = no_of_nodes - 1;
                }
            } else {
                to = (int) floor((sqrt(8 * next_a + 1) + 1) / 2);
                from = (int) (next_a - (((igraph_real_t)to) * (to - 1)) / 2);
            }
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, from));
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, to));
            p_a++; UPD_A();

        }
    }

    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&add);
    igraph_vector_destroy(&delete);
    IGRAPH_FINALLY_CLEAN(3);

    if (permutation) {
        int newec = igraph_vector_size(&newedges);
        for (i = 0; i < newec; i++) {
            int tmp = VECTOR(newedges)[i];
            VECTOR(newedges)[i] = VECTOR(*permutation)[tmp];
        }
    }

    IGRAPH_CHECK(igraph_create(new_graph, &newedges, no_of_nodes, directed));

    igraph_vector_destroy(&newedges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

#undef D_CODE
#undef U_CODE
#undef CODE
#undef CODEE
#undef UPD_E
#undef UPD_A
#undef UPD_D

/**
 * \function igraph_correlated_pair_game
 * Generate pairs of correlated random graphs
 *
 * Sample two random graphs, with given correlation.
 *
 * \param graph1 The first graph will be stored here.
 * \param graph2 The second graph will be stored here.
 * \param n The number of vertices in both graphs.
 * \param corr A scalar in the unit interval, the target Pearson
 *        correlation between the adjacency matrices of the original the
 *        generated graph (the adjacency matrix being used as a vector).
 * \param p A numeric scalar, the probability of an edge between two
 *        vertices, it must in the open (0,1) interval.
 * \param directed Whether to generate directed graphs.
 * \param permutation A permutation to apply to the vertices of the
 *        second graph. It can also be a null pointer, in which case
 *        the vertices will not be permuted.
 * \return Error code
 *
 * \sa \ref igraph_correlated_game() for generating a correlated pair
 * to a given graph.
 */

int igraph_correlated_pair_game(igraph_t *graph1, igraph_t *graph2,
                                int n, igraph_real_t corr, igraph_real_t p,
                                igraph_bool_t directed,
                                const igraph_vector_t *permutation) {

    IGRAPH_CHECK(igraph_erdos_renyi_game(graph1, IGRAPH_ERDOS_RENYI_GNP, n, p,
                                         directed, IGRAPH_NO_LOOPS));
    IGRAPH_CHECK(igraph_correlated_game(graph1, graph2, corr, p, permutation));
    return 0;
}


/* Uniform sampling of labelled trees (igraph_tree_game) */

/* The following implementation uniformly samples Prufer trees and converts
 * them to trees.
 */

static int igraph_i_tree_game_prufer(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {
    igraph_vector_int_t prufer;
    long i;

    if (directed) {
        IGRAPH_ERROR("The Prufer method for random tree generation does not support directed trees", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_init(&prufer, n - 2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &prufer);

    RNG_BEGIN();

    for (i = 0; i < n - 2; ++i) {
        VECTOR(prufer)[i] = RNG_INTEGER(0, n - 1);
    }

    RNG_END();

    IGRAPH_CHECK(igraph_from_prufer(graph, &prufer));

    igraph_vector_int_destroy(&prufer);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* The following implementation is based on loop-erased random walks and Wilson's algorithm
 * for uniformly sampling spanning trees. We effectively sample spanning trees of the complete
 * graph.
 */

/* swap two elements of a vector_int */
#define SWAP_INT_ELEM(vec, i, j) \
    { \
        igraph_integer_t temp; \
        temp = VECTOR(vec)[i]; \
        VECTOR(vec)[i] = VECTOR(vec)[j]; \
        VECTOR(vec)[j] = temp; \
    }

static int igraph_i_tree_game_loop_erased_random_walk(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {
    igraph_vector_t edges;
    igraph_vector_int_t vertices;
    igraph_vector_bool_t visited;
    long i, j, k;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * (n - 1));

    IGRAPH_CHECK(igraph_vector_bool_init(&visited, n));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &visited);

    /* The vertices vector contains visited vertices between 0..k-1, unvisited ones between k..n-1. */
    IGRAPH_CHECK(igraph_vector_int_init_seq(&vertices, 0, n - 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vertices);

    RNG_BEGIN();

    /* A simple implementation could be as below. This is for illustration only.
     * The actually implemented algorithm avoids unnecessary walking on the already visited
     * portion of the vertex set.
     */
    /*
    // pick starting point for the walk
    i = RNG_INTEGER(0, n-1);
    VECTOR(visited)[i] = 1;

    k=1;
    while (k < n) {
        // pick next vertex in the walk
        j = RNG_INTEGER(0, n-1);
        // if it has not been visited before, connect to the previous vertex in the sequence
        if (! VECTOR(visited)[j]) {
            VECTOR(edges)[2*k - 2] = i;
            VECTOR(edges)[2*k - 1] = j;
            VECTOR(visited)[j] = 1;
            k++;
        }
        i=j;
    }
    */

    i = RNG_INTEGER(0, n - 1);
    VECTOR(visited)[i] = 1;
    SWAP_INT_ELEM(vertices, 0, i);

    for (k = 1; k < n; ++k) {
        j = RNG_INTEGER(0, n - 1);
        if (VECTOR(visited)[VECTOR(vertices)[j]]) {
            i = VECTOR(vertices)[j];
            j = RNG_INTEGER(k, n - 1);
        }
        VECTOR(visited)[VECTOR(vertices)[j]] = 1;
        SWAP_INT_ELEM(vertices, k, j);
        VECTOR(edges)[2 * k - 2] = i;
        i = VECTOR(vertices)[k];
        VECTOR(edges)[2 * k - 1] = i;
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&vertices);
    igraph_vector_bool_destroy(&visited);
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

#undef SWAP_INT_ELEM

/**
 * \ingroup generators
 * \function igraph_tree_game
 * \brief Generates a random tree with the given number of nodes
 *
 * This function samples uniformly from the set of labelled trees,
 * i.e. it can generate each labelled tree with the same probability.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of nodes in the tree.
 * \param directed Whether to create a directed tree. The edges are oriented away from the root.
 * \param method The algorithm to use to generate the tree. Possible values:
 *        \clist
 *        \cli IGRAPH_RANDOM_TREE_PRUFER
 *          This algorithm samples Pr&uuml;fer sequences unformly, then converts them to trees.
 *          Directed trees are not currently supported.
 *        \cli IGRAPH_RANDOM_LERW
 *          This algorithm effectively performs a loop-erased random walk on the complete graph
 *          to uniformly sample its spanning trees (Wilson's algorithm).
 *        \endclist
 * \return Error code:
 *          \c IGRAPH_ENOMEM: there is not enough
 *           memory to perform the operation.
 *          \c IGRAPH_EINVAL: invalid tree size
 *
 * \sa \ref igraph_from_prufer()
 *
 */

int igraph_tree_game(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, igraph_random_tree_t method) {
    if (n < 2) {
        IGRAPH_CHECK(igraph_empty(graph, n, directed));
        return IGRAPH_SUCCESS;
    }

    switch (method) {
    case IGRAPH_RANDOM_TREE_PRUFER:
        return igraph_i_tree_game_prufer(graph, n, directed);
    case IGRAPH_RANDOM_TREE_LERW:
        return igraph_i_tree_game_loop_erased_random_walk(graph, n, directed);
    default:
        IGRAPH_ERROR("Invalid method for random tree construction", IGRAPH_EINVAL);
    }
}
