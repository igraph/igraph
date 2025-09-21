/*
   igraph library.
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

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_random.h"

#include "core/interruption.h"
#include "internal/utils.h"
#include "math/safe_intop.h"
#include "misc/graphicality.h"
#include "random/random_internal.h"

/**
 * \section about_erdos_renyi
 *
 * <para>
 * There are two classic random graph models referred to as the Erdős-Rényi
 * random graph, or sometimes simply \em the random graph. Both fix the vertex
 * count n, but while the G(n,m) model prescribes precisely m edges, the G(n,p)
 * model connects all vertex pairs independently with probability p. While
 * these models look superficially different, when n is large they behave in
 * a similar manner. G(n,m) graphs have a density of exactly
 * <code>p = m / m_max</code>, while G(n,p) graphs have <code>m = p m_max</code>
 * edges on \em average, where \c m_max is the number of vertex pairs. Indeed,
 * these two models turns out to be two sides of the same coin: both can be
 * understood as maximum entropy models with a constraint on the number of
 * edges. The G(n,m) is obtained from a sharp constraint, while G(n,p) from
 * an average constraint (soft constraint).
 * </para>
 *
 * <para>
 * The maximum entropy framework allows for rigorous generalizations of these
 * models to various scenarios, of which igraph supports many, such as models
 * defined over directed graphs, bipartite graphs, multigraphs, or even over
 * edge-labelled graphs. Constraining edge counts between various subsets of
 * vertices yields further families of related models, such as
 * \ref igraph_sbm_game() (given connection probabilities between categories)
 * or \ref igraph_degree_sequence_game() (given incident edge counts, i.e.
 * degrees, for each vertex).
 * </para>
 */


static igraph_error_t iea_game(
        igraph_t *graph,
        igraph_int_t n, igraph_int_t m,
        igraph_bool_t directed, igraph_bool_t loops) {

    igraph_vector_int_t edges;
    int iter = 0;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, m * 2));

    for (igraph_int_t i = 0; i < m; i++) {
        igraph_int_t from, to;
        from = RNG_INTEGER(0, n - 1);
        if (loops) {
            to = RNG_INTEGER(0, n - 1);
        } else {
            to = RNG_INTEGER(0, n - 2);
            if (from == to) {
                to = n - 1;
            }
        }
        igraph_vector_int_push_back(&edges, from); /* reserved */
        igraph_vector_int_push_back(&edges, to); /* reserved */
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Uniform sampling of multigraphs from G(n,m) */
static igraph_error_t gnm_multi(
        igraph_t *graph,
        igraph_int_t n, igraph_int_t m,
        igraph_bool_t directed, igraph_bool_t loops) {

    /* Conceptually, uniform multigraph sampling works as follows:
     *
     *  - Consider a list containing all vertex pairs. Use unordered pairs for
     *    undirected graphs, ordered pairs for directed ones, and include
     *    self-pairs if desired.
     *  - Pick a list element uniformly at random with replacement and append
     *    it to the list. Add the corresponding edge to the graph.
     *  - Continue picking elements form the list uniformly at random and
     *    appending the pick to the list until we have sampled the desired
     *    number of edges.
     *
     * Let's illustrate how this is implemented on the directed case with loops
     * allowed. Each element of the n*n adjacency matrix corresponds to an
     * ordered vertex pair. Uniformly sampling a matrix element is possible
     * by generating two random matrix indices. As an analog of appending to the
     * list of pairs, we extend the matrix with additional rows. The last row
     * will generally be incomplete, so we use rejection sampling to avoid
     * exceeding its length.
     *
     * In the undirected case, we still sample _ordered_ vertex pairs for
     * simplicity. To account for the resulting duplication, we append the
     * sampled pair _twice_. In the undirected case with loops, we must also
     * duplicate the matrix diagonal.
     */

    igraph_vector_int_t edges;
    igraph_int_t nrow, ncol;
    igraph_int_t last; /* column index of last element in last row */
    int iter = 0;

    /* Constraining n and m by IGRAPH_VCOUNT_MAX and IGRAPH_ECOUNT_MAX,
     * as done in the caller, is sufficient to prevent overflow,
     * except for the one special case below:
     */

    if (!directed && !loops &&
        n == IGRAPH_VCOUNT_MAX && m == IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Too many edges or vertices for G(n,m) multigraph.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 2*m);

    if (directed && loops) {
        nrow = ncol = n;
        last = ncol-1;
        for (igraph_int_t i=0; i < m; i++) {
            while (true) {
                igraph_int_t r = RNG_INTEGER(0, nrow-1);
                igraph_int_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_int_t j = (r - n) * ncol + c;
                    if (IGRAPH_UNLIKELY(j >= i)) continue; /* rejection sampling */
                    VECTOR(edges)[2*i]   = VECTOR(edges)[2*j];
                    VECTOR(edges)[2*i+1] = VECTOR(edges)[2*j+1];
                } else {
                    VECTOR(edges)[2*i]   = r;
                    VECTOR(edges)[2*i+1] = c;
                }

                last += 1;
                if (last >= ncol) {
                    last -= ncol;
                    nrow++;
                }

                break;
            }
            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }
    } else if (directed && !loops) {
        nrow = n;
        ncol = n-1;
        last = ncol-1;
        for (igraph_int_t i=0; i < m; i++) {
            while (true) {
                igraph_int_t r = RNG_INTEGER(0, nrow-1);
                igraph_int_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_int_t j = (r - n) * ncol + c;
                    if (IGRAPH_UNLIKELY(j >= i)) continue; /* rejection sampling */
                    VECTOR(edges)[2*i]   = VECTOR(edges)[2*j];
                    VECTOR(edges)[2*i+1] = VECTOR(edges)[2*j+1];
                } else {

                    /* Eliminate self-loops. */
                    if (c == r) {
                        c = n-1;
                    }

                    VECTOR(edges)[2*i]   = r;
                    VECTOR(edges)[2*i+1] = c;
                }

                last += 1;
                if (last >= ncol) {
                    last -= ncol;
                    nrow++;
                }

                break;
            }
            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }
    } else if (!directed && loops) {
        nrow = n;
        ncol = n+1;
        last = ncol-1;
        for (igraph_int_t i=0; i < m; i++) {
            while (true) {
                igraph_int_t r = RNG_INTEGER(0, nrow-1);
                igraph_int_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_int_t j = (r - n) * ncol + c;
                    if (IGRAPH_UNLIKELY(j >= 2*i)) continue; /* rejection sampling */
                    VECTOR(edges)[2*i]   = VECTOR(edges)[2*(j/2)];
                    VECTOR(edges)[2*i+1] = VECTOR(edges)[2*(j/2)+1];
                } else {

                    /* Two chances to sample from matrix diagonal,
                     * when c == r and when c == n. */
                    if (c == n) {
                        c = r;
                    }

                    VECTOR(edges)[2*i]   = r;
                    VECTOR(edges)[2*i+1] = c;
                }

                last += 2;
                while (last >= ncol) {
                    last -= ncol;
                    nrow++;
                }

                break;
            }
            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }
    } else  /* !directed && !loops */ {
        nrow = n;
        ncol = n-1;
        last = ncol-1;
        for (igraph_int_t i=0; i < m; i++) {
            while (true) {
                igraph_int_t r = RNG_INTEGER(0, nrow-1);
                igraph_int_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_int_t j = (r - n) * ncol + c;
                    if (IGRAPH_UNLIKELY(j >= 2*i)) continue; /* rejection sampling */
                    VECTOR(edges)[2*i]   = VECTOR(edges)[2*(j/2)];
                    VECTOR(edges)[2*i+1] = VECTOR(edges)[2*(j/2)+1];
                } else {

                    /* Eliminate self-loops. */
                    if (c == r) {
                        c = n-1;
                    }

                    VECTOR(edges)[2*i]   = r;
                    VECTOR(edges)[2*i+1] = c;
                }

                last += 2;
                while (last >= ncol) {
                    last -= ncol;
                    nrow++;
                }

                break;
            }
            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Uniform sampling of simple graphs (with loops) from G(n,m) */
static igraph_error_t gnm_simple(
        igraph_t *graph,
        igraph_int_t n, igraph_int_t m,
        igraph_bool_t directed, igraph_bool_t loops,
        igraph_bool_t edge_labeled) {

    /* This function uses doubles in its `s` vector, and for `maxedges` and `last`.
     * This is because on a system with 32-bit ints, maxedges will be larger than
     * IGRAPH_INTEGER_MAX and this will cause overflows when calculating `from` and `to`
     * for tests on large graphs. This is also why we need a 'real' version of random_sample.
    */
    igraph_real_t n_real = (igraph_real_t) n; /* for divisions below */
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int iter = 0;


    igraph_real_t maxedges = n;
    if (directed && loops) {
        maxedges *= n;
    } else if (directed && !loops) {
        maxedges *= (n - 1);
    } else if (!directed && loops) {
        maxedges *= (n + 1) / 2.0;
    } else {
        maxedges *= (n - 1) / 2.0;
    }

    if (m > maxedges) {
        IGRAPH_ERROR(
                "Too many edges requested compared to the number of vertices for G(n,m) model.",
                IGRAPH_EINVAL);
    }

    if (maxedges == m && ! edge_labeled) {
        /* TODO: Cannot use igraph_full() when edge_labeled as we must shuffle edges. */
        IGRAPH_CHECK(igraph_full(graph, n, directed, loops));
    } else {
        igraph_int_t slen;

        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_i_random_sample_real(&s, 0, maxedges - 1, m));

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, igraph_vector_size(&s) * 2));

        slen = igraph_vector_size(&s);
        if (directed && loops) {
            for (igraph_int_t i = 0; i < slen; i++) {
                igraph_int_t to = trunc(VECTOR(s)[i] / n_real);
                igraph_int_t from = VECTOR(s)[i] - to * n_real;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (directed && !loops) {
            for (igraph_int_t i = 0; i < slen; i++) {
                igraph_int_t from = trunc(VECTOR(s)[i] / (n_real - 1));
                igraph_int_t to = VECTOR(s)[i] - from * (n_real - 1);
                if (from == to) {
                    to = n - 1;
                }
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (!directed && loops) {
            for (igraph_int_t i = 0; i < slen; i++) {
                igraph_int_t to = trunc((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                igraph_int_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else { /* !directed && !loops */
            for (igraph_int_t i = 0; i < slen; i++) {
                igraph_int_t to = trunc((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                igraph_int_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        }

        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);

        if (edge_labeled) {
            IGRAPH_CHECK(igraph_i_vector_int_shuffle_pairs(&edges));
        }

        IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_erdos_renyi_game_gnm
 * \brief Generates a random (Erdős-Rényi) graph with a fixed number of edges.
 *
 * In the <code>G(n, m)</code> Erdős-Rényi model, a graph with \p n vertices
 * and \p m edges is generated uniformly at random.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param m The number of edges in the graph.
 * \param directed Whether to generate a directed graph.
 * \param allowed_edge_types Controls whether multi-edges and self-loops
 *     are generated. See \ref igraph_edge_type_sw_t.
 * \param edge_labeled If true, the sampling is done uniformly from the set
 *     of ordered edge lists. See \ref igraph_iea_game() for more information.
 *     Set this to \c false to select the classic Erdős-Rényi model.
 *     The constants \c IGRAPH_EDGE_UNLABELED and \c IGRAPH_EDGE_LABELED
 *     may be used instead of \c false and \c true for better readability.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid \p n or \p m parameter.
 *         \c IGRAPH_ENOMEM: there is not enough memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_erdos_renyi_game_gnp() to sample from the related
 * <code>G(n, p)</code> model, which constrains the \em expected edge count;
 * \ref igraph_iea_game() to generate multigraph by assigning edges to vertex
 * pairs uniformly and independently;
 * \ref igraph_degree_sequence_game() to constrain the degree sequence;
 * \ref igraph_bipartite_game_gnm() for the bipartite version of this model;
 * \ref igraph_barabasi_game() and \ref igraph_growing_random_game() for other
 * commonly used random graph models.
 *
 * \example examples/simple/igraph_erdos_renyi_game_gnm.c
 */
igraph_error_t igraph_erdos_renyi_game_gnm(
        igraph_t *graph,
        igraph_int_t n, igraph_int_t m,
        igraph_bool_t directed,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_bool_t edge_labeled) {

    igraph_bool_t loops, multiple;

    /* The multigraph implementation relies on the below checks to avoid overflow. */
    if (n < 0 || n > IGRAPH_VCOUNT_MAX) {
        IGRAPH_ERROR("Invalid number of vertices for G(n,m) model.", IGRAPH_EINVAL);
    }
    if (m < 0 || m > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Invalid number of edges for G(n,m) model.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_edge_type_to_loops_multiple(allowed_edge_types, &loops, &multiple));

    /* Special cases of "too many edges" that also apply to multigraphs:
     *  - The null graph cannot have edges.
     *  - The singleton graph cannot have edges unless loops are allowed.
     */
    if (m > 0 && ((n == 0) || (!loops && n == 1))) {
        IGRAPH_ERROR(
            "Too many edges requested compared to the number of vertices for G(n,m) model.",
             IGRAPH_EINVAL);
    }

    if (m == 0) {
        return igraph_empty(graph, n, directed);
    }

    if (edge_labeled) {
        if (multiple) {
            return iea_game(graph, n, m, directed, loops);
        } else {
            return gnm_simple(graph, n, m, directed, loops, /*edge_labeled=*/ true);
        }
    } else {
        if (multiple) {
            return gnm_multi(graph, n, m, directed, loops);
        } else {
            return gnm_simple(graph, n, m, directed, loops, /*edge_labeled=*/ false);
        }
    }
}


/**
 * \ingroup generators
 * \function igraph_iea_game
 * \brief Generates a random multigraph through independent edge assignment.
 *
 * \experimental
 *
 * This model generates random multigraphs on \p n vertices with \p m edges
 * through independent edge assignment (IEA). Each of the \p m edges is assigned
 * uniformly at random to an \em ordered vertex pair, independently of each
 * other.
 *
 * </para><para>
 * This model does not sample multigraphs uniformly. Undirected graphs are
 * generated with probability proportional to
 *
 * </para><para>
 * <code>(prod_(i&lt;j) A_ij ! prod_i A_ii !!)^(-1)</code>,
 *
 * </para><para>
 * where \c A denotes the adjacency matrix and <code>!!</code> denotes
 * the double factorial. Here \c A is assumed to have twice the number of
 * self-loops on its diagonal. The corresponding  expression for directed
 * graphs is
 *
 * </para><para>
 * <code>(prod_(i,j) A_ij !)^(-1)</code>.
 *
 * </para><para>
 * Thus the probability of all simple graphs (which only have 0s and 1s in
 * the adjacency matrix) is the same, while that of non-simple ones depends
 * on their edge and self-loop multiplicities.
 *
 * </para><para>
 * An alternative way to think of this model is that it performs uniform
 * sampling of \em edge-labeled graphs, i.e. graphs in which not only vertices,
 * but also edges carry unique identities and are distinguishable.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param m The number of edges in the graph.
 * \param directed Whether to generate a directed graph.
 * \param loops Whether to generate self-loops.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid \p n or \p m parameter.
 *         \c IGRAPH_ENOMEM: there is not enough
 *         memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_erdos_renyi_game_gnm() to uniformly sample graphs with
 * a given number of vertices and edges.
 */
igraph_error_t igraph_iea_game(
        igraph_t *graph,
        igraph_int_t n, igraph_int_t m,
        igraph_bool_t directed, igraph_bool_t loops) {

    igraph_edge_type_sw_t allowed_edge_types = IGRAPH_MULTI_SW;
    if (loops) {
        allowed_edge_types |= IGRAPH_LOOPS_SW;
    }
    return igraph_erdos_renyi_game_gnm(graph, n, m, directed, allowed_edge_types, true);
}

/* This G(n,p) implementation is used only with very large vertex counts, above
 * sqrt(MAX_EXACT_REAL) ~ 100 million, when the default implementation would
 * fail due to overflow. While this version avoids overflow and uses less memory,
 * it is also slower than the default implementation.
 *
 * This function expects that when multiple=true, the p parameter has already
 * been transformed by p = p / (1 + p). This is currently done by the caller.
 */
static igraph_error_t gnp_large(
    igraph_t *graph, igraph_int_t n, igraph_real_t p,
    igraph_bool_t directed, igraph_bool_t loops, igraph_bool_t multiple,
    igraph_int_t ecount_estimate
) {

    igraph_vector_int_t edges;
    int iter = 0;

    /* Necessitated by floating point arithmetic used in the implementation. */
    if (n >= IGRAPH_MAX_EXACT_REAL) {
        IGRAPH_ERROR("Number of vertices is too large.", IGRAPH_EOVERFLOW);
    }

    if (ecount_estimate > IGRAPH_ECOUNT_MAX) {
        ecount_estimate = IGRAPH_ECOUNT_MAX;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, 2*ecount_estimate));

    for (igraph_int_t i=0; i < n; i++) {
        igraph_int_t j = directed ? 0 : i;

        while (true) {
            igraph_real_t gap = RNG_GEOM(p);

            /* This formulation not only terminates the loop when necessary,
             * but also protects against overflow when 'p' is very small
             * and 'gap' becomes very large, perhaps larger than representable
             * in an igraph_int_t. */
            if (gap >= n - j) {
                break;
            }

            j += gap;

            if (loops || i != j) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
            }

            j += ! multiple; /* 1 for simple graph, 0 for multigraph */

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t gnp_edge_labeled(
        igraph_t *graph,
        igraph_int_t n, igraph_real_t p,
        igraph_bool_t directed,
        igraph_bool_t loops, igraph_bool_t multiple) {

    if (multiple) {
        igraph_real_t maxedges = n;

        if (directed && loops) {
            maxedges *= n;
        } else if (directed && !loops) {
            maxedges *= (n - 1);
        } else if (!directed && loops) {
            maxedges *= (n + 1) / 2.0;
        } else {
            maxedges *= (n - 1) / 2.0;
        }

        igraph_real_t m;
        do {
            m = RNG_GEOM( 1.0 / (1.0 + maxedges * p) );
        } while (m > (igraph_real_t) IGRAPH_INTEGER_MAX);

        return iea_game(graph, n, m, directed, loops);
    } else {
        IGRAPH_ERROR("The edge-labeled G(n,p) model is not yet implemented for graphs without multi-edges.",
                     IGRAPH_UNIMPLEMENTED);
    }
}

/**
 * \ingroup generators
 * \function igraph_erdos_renyi_game_gnp
 * \brief Generates a random (Erdős-Rényi) graph with fixed edge probabilities.
 *
 * In the <code>G(n, p)</code> Erdős-Rényi model, also known as the Gilbert model,
 * or Bernoulli random graph, a graph with \p n vertices is generated such that
 * every possible edge is included in the graph independently with probability
 * \p p. This is equivalent to a maximum entropy random graph model model with
 * a constraint on the \em expected edge count. The maximum entropy view allows
 * for extending the model to multigraphs, as discussed by Park and Newman (2004),
 * section III.D. In this case, \p p is interpreted as the expected number of
 * edges between any vertex pair.
 *
 * </para><para>
 * Setting <code>p = 1/2</code> and <code>multiple = false</code> generates all
 * graphs without multi-edges on \p n vertices with the same probability.
 *
 * </para><para>
 * For both simple and multigraphs, the expected mean degree of the graph is
 * approximately <code>p n</code>; set <code>p = k/n</code> when a mean degree
 * of approximately \c k is desired. More precisely, the expected mean degree is
 * <code>p(n-1)</code> in (undirected or directed) graphs without self-loops,
 * <code>p(n+1)</code> in undirected graphs with self-loops, and
 * <code>p n</code> in directed graphs with self-loops.
 *
 * </para><para>
 * When generating multigraphs, the distribution of the edge multiplicities is
 * geometric, i.e. the probability of finding \c m edges between two vertices
 * is <code>q (1-q)^m</code>, where <code>q = 1 / (1+p)</code>.
 *
 * </para><para>
 * This function uses the sequential geometric sampling technique described in
 * Batagelj and Brandes (2005), with a modification to handle multigraphs.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * J. Park and M. E. J. Newman: "Statistical Mechanics of Networks".
 * Phys. Rev. E 70, 066117 (2004).
 * https://doi.org/10.1103/PhysRevE.70.066117
 *
 * </para><para>
 * V. Batagelj and U. Brandes: "Efficient Generation of Large Random Networks".
 * Phys. Rev. E 71, 036113 (2005).
 * https://doi.org/10.1103/PhysRevE.71.036113
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param p The expected number of edges between any vertex pair.
 *    When multi-edges are disallowed, this is equivalent to the probability
 *    of having a connection between any two vertices.
 * \param directed Whether to generate a directed graph.
 * \param allowed_edge_types Controls whether multi-edges and self-loops
 *     are generated. See \ref igraph_edge_type_sw_t.
 * \param edge_labeled If true, the model is defined over the set of ordered
 *     edge lists, i.e. over the set of edge-labeled graphs. Set it to
 *     \c false to select the classic Erdős-Rényi model.
 *     The constants \c IGRAPH_EDGE_UNLABELED and \c IGRAPH_EDGE_LABELED
 *     may be used instead of \c false and \c true for better readability.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid \p n or \p p parameter.
 *         \c IGRAPH_ENOMEM: there is not enough memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_erdos_renyi_game_gnm() to generate random graphs with
 * a sharply fixed edge count; \ref igraph_chung_lu_game() and
 * \ref igraph_static_fitness_game() to generate random graphs with a
 * fixed expected degree sequence; \ref igraph_bipartite_game_gnm() for the
 * bipartite version of this model; \ref igraph_barabasi_game() and
 * \ref igraph_growing_random_game() for other commonly used random graph models.
 *
 * \example examples/simple/igraph_erdos_renyi_game_gnp.c
 */
igraph_error_t igraph_erdos_renyi_game_gnp(
        igraph_t *graph,
        igraph_int_t n, igraph_real_t p,
        igraph_bool_t directed,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_bool_t edge_labeled) {

    /* This function uses doubles in its `s` vector, and for `maxedges` and `last`.
     * This is because on a system with 32-bit ints, maxedges will be larger than
     * IGRAPH_INTEGER_MAX and this will cause overflows when calculating `from` and `to`
     * for tests on large graphs.
    */
    igraph_int_t no_of_nodes = n;
    igraph_real_t no_of_nodes_real = (igraph_real_t) no_of_nodes;   /* for divisions below */
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    igraph_bool_t loops, multiple;
    int iter = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices for G(n,p) model.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_edge_type_to_loops_multiple(allowed_edge_types, &loops, &multiple));

    if (multiple) {
        if (p < 0.0) {
            IGRAPH_ERROR("Invalid expected edge multiplicity given for G(n,p) multigraph model.", IGRAPH_EINVAL);
        }
    } else {
        if (p < 0.0 || p > 1.0) {
            IGRAPH_ERROR("Invalid probability given for G(n,p) model.", IGRAPH_EINVAL);
        }
    }

    if (edge_labeled) {
        return gnp_edge_labeled(graph, n, p, directed, loops, multiple);
    }

    if (multiple) {
        /* Convert the expected edge count to the appropriate probability parameter
         * of the geometric distribution when sampling lengths of runs of 0s in the
         * adjacency matrix. */
        p = p / (1 + p);
    }

    if (p == 0.0 || no_of_nodes == 0) {
        IGRAPH_CHECK(igraph_empty(graph, n, directed));
    } else if (! multiple && p == 1.0) {
        IGRAPH_CHECK(igraph_full(graph, n, directed, loops));
    } else {
        igraph_real_t maxedges = n, last;
        igraph_int_t ecount_estimate, ecount;

        if (directed && loops) {
            maxedges *= n;
        } else if (directed && !loops) {
            maxedges *= (n - 1);
        } else if (!directed && loops) {
            maxedges *= (n + 1) / 2.0;
        } else {
            maxedges *= (n - 1) / 2.0;
        }

        IGRAPH_CHECK(igraph_i_safe_floor(maxedges * p * 1.1, &ecount_estimate));

        if (maxedges > IGRAPH_MAX_EXACT_REAL) {
            /* Use a slightly slower, but overflow-free implementation. */
            return gnp_large(graph, n, p, directed, loops, multiple, ecount_estimate);
        }

        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&s, ecount_estimate));

        last = RNG_GEOM(p);
        while (last < maxedges) {
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            last += RNG_GEOM(p);
            last += ! multiple; /* 1 for simple graph, 0 for multigraph */
            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }

        ecount = igraph_vector_size(&s);
        if (ecount > IGRAPH_ECOUNT_MAX) {
            IGRAPH_ERROR("Overflow in number of edges.", IGRAPH_EOVERFLOW);
        }

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, 2*ecount));

        iter = 0;
        if (directed && loops) {
            for (igraph_int_t i = 0; i < ecount; i++) {
                igraph_int_t to = trunc(VECTOR(s)[i] / no_of_nodes_real);
                igraph_int_t from = VECTOR(s)[i] - to * no_of_nodes_real;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (directed && !loops) {
            for (igraph_int_t i = 0; i < ecount; i++) {
                igraph_int_t to = trunc(VECTOR(s)[i] / no_of_nodes_real);
                igraph_int_t from = VECTOR(s)[i] - to * no_of_nodes_real;
                if (from == to) {
                    to = no_of_nodes - 1;
                }
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (!directed && loops) {
            for (igraph_int_t i = 0; i < ecount; i++) {
                igraph_int_t to = trunc((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                igraph_int_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else { /* !directed && !loops */
            for (igraph_int_t i = 0; i < ecount; i++) {
                igraph_int_t to = trunc((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                igraph_int_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        }

        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
