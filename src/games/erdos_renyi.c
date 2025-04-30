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

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_random.h"

#include "core/interruption.h"
#include "math/safe_intop.h"
#include "random/random_internal.h"

/**
 * \section about_games
 *
 * <para>Games are random graph generators, i.e. they generate a different
 * graph every time they are called. igraph includes many such generators.
 * Some implement stochastic graph construction processes inspired by real-world
 * mechanics, such as preferential attachment, while others are designed to
 * produce graphs with certain used properties (e.g. fixed number of edges,
 * fixed degrees, etc.)</para>
 */

/* This implementation is used only with very large vertex counts, above
 * sqrt(MAX_EXACT_REAL) ~ 100 million, when the default implementation would
 * fail due to overflow. While this version avoids overflow and uses less memory,
 * it is also slower than the default implementation. */
static igraph_error_t gnp_large(
    igraph_t *graph, igraph_integer_t n, igraph_real_t p,
    igraph_bool_t directed, igraph_bool_t loops, igraph_integer_t ecount_estimate
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

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < n; i++) {
        igraph_integer_t j = directed ? 0 : i;

        while (true) {
            igraph_real_t gap = RNG_GEOM(p);

            /* This formulation not only terminates the loop when necessary,
             * but also protects against overflow when 'p' is very small
             * and 'gap' becomes very large, perhaps larger than representable
             * in an igraph_integer_t. */
            if (gap >= n - j) {
                break;
            }

            j += gap;

            if (loops || i != j) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
            }

            j++;

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }
    }
    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
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
 * a constraint on the \em expected edge count. Setting <code>p = 1/2</code>
 * generates all graphs on \p n vertices with the same probability.
 *
 * </para><para>
 * The expected mean degree of the graph is approximately <code>p n</code>;
 * set <code>p = k/n</code> when a mean degree of approximately \c k is
 * desired. More precisely, the expected mean degree is <code>p(n-1)</code>
 * in (undirected or directed) graphs without self-loops,
 * <code>p(n+1)</code> in undirected graphs with self-loops, and
 * <code>p n</code> in directed graphs with self-loops.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param p The probability of the existence of an edge in the graph.
 * \param directed Whether to generate a directed graph.
 * \param loops Whether to generate self-loops.
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
    igraph_t *graph, igraph_integer_t n, igraph_real_t p,
    igraph_bool_t directed, igraph_bool_t loops
) {
    /* This function uses doubles in its `s` vector, and for `maxedges` and `last`.
     * This is because on a system with 32-bit ints, maxedges will be larger than
     * IGRAPH_INTEGER_MAX and this will cause overflows when calculating `from` and `to`
     * for tests on large graphs.
    */
    igraph_integer_t no_of_nodes = n;
    igraph_real_t no_of_nodes_real = (igraph_real_t) no_of_nodes;   /* for divisions below */
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int iter = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices for G(n,p) model.", IGRAPH_EINVAL);
    }
    if (p < 0.0 || p > 1.0) {
        IGRAPH_ERROR("Invalid probability given for G(n,p) model.", IGRAPH_EINVAL);
    }

    if (p == 0.0 || no_of_nodes == 0) {
        IGRAPH_CHECK(igraph_empty(graph, n, directed));
    } else if (p == 1.0) {
        IGRAPH_CHECK(igraph_full(graph, n, directed, loops));
    } else {
        igraph_real_t maxedges = n, last;
        igraph_integer_t ecount_estimate, ecount;

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
            return gnp_large(graph, n, p, directed, loops, ecount_estimate);
        }

        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&s, ecount_estimate));

        RNG_BEGIN();

        last = RNG_GEOM(p);
        while (last < maxedges) {
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            last += RNG_GEOM(p);
            last += 1;
            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }

        RNG_END();

        ecount = igraph_vector_size(&s);
        if (ecount > IGRAPH_ECOUNT_MAX) {
            IGRAPH_ERROR("Overflow in number of edges.", IGRAPH_EOVERFLOW);
        }

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, 2*ecount));

        iter = 0;
        if (directed && loops) {
            for (igraph_integer_t i = 0; i < ecount; i++) {
                igraph_integer_t to = floor(VECTOR(s)[i] / no_of_nodes_real);
                igraph_integer_t from = VECTOR(s)[i] - to * no_of_nodes_real;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (directed && !loops) {
            for (igraph_integer_t i = 0; i < ecount; i++) {
                igraph_integer_t to = floor(VECTOR(s)[i] / no_of_nodes_real);
                igraph_integer_t from = VECTOR(s)[i] - to * no_of_nodes_real;
                if (from == to) {
                    to = no_of_nodes - 1;
                }
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (!directed && loops) {
            for (igraph_integer_t i = 0; i < ecount; i++) {
                igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                igraph_integer_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else { /* !directed && !loops */
            for (igraph_integer_t i = 0; i < ecount; i++) {
                igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                igraph_integer_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2;
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
 * Thus the probability of all simple graphs (which only have 0s and 1s
 * in the adjacency matrix) is the same, while that of
 * non-simple ones depends on their edge and self-loop multiplicities.
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
        igraph_integer_t n, igraph_integer_t m,
        igraph_bool_t directed, igraph_bool_t loops) {

    igraph_vector_int_t edges;
    int iter = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices for IEA model.", IGRAPH_EINVAL);
    }
    if (m < 0 || m > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Invalid number of edges for IEA model.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, m * 2));

    RNG_BEGIN();
    for (igraph_integer_t i = 0; i < m; i++) {
        igraph_integer_t from, to;
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
    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Uniform sampling of multigraphs from G(n,m) */
static igraph_error_t gnm_multi(
        igraph_t *graph,
        igraph_integer_t n, igraph_integer_t m,
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
    igraph_integer_t nrow, ncol;
    igraph_integer_t last; /* column index of last element in last row */
    int iter = 0;

    /* Constraining n and m by IGRAPH_VCOUNT_MAX and IGRAPH_ECOUNT_MAX,
     * as done in the caller, is sufficient to prevent overflow,
     * except for the one special case below:
     */

    if (!directed && !loops &&
        n == IGRAPH_VCOUNT_MAX && m == IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Too many edges or vertices.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 2*m);

    RNG_BEGIN();
    if (directed && loops) {
        nrow = ncol = n;
        last = ncol-1;
        for (igraph_integer_t i=0; i < m; i++) {
            while (true) {
                igraph_integer_t r = RNG_INTEGER(0, nrow-1);
                igraph_integer_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_integer_t j = (r - n) * ncol + c;
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
        for (igraph_integer_t i=0; i < m; i++) {
            while (true) {
                igraph_integer_t r = RNG_INTEGER(0, nrow-1);
                igraph_integer_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_integer_t j = (r - n) * ncol + c;
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
        for (igraph_integer_t i=0; i < m; i++) {
            while (true) {
                igraph_integer_t r = RNG_INTEGER(0, nrow-1);
                igraph_integer_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_integer_t j = (r - n) * ncol + c;
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
        for (igraph_integer_t i=0; i < m; i++) {
            while (true) {
                igraph_integer_t r = RNG_INTEGER(0, nrow-1);
                igraph_integer_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n) {
                    igraph_integer_t j = (r - n) * ncol + c;
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
    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

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
 * \param loops Whether to generate self-loops.
 * \param multiple Whether it is allowed to generate more than one edge between
 *    the same pair of vertices.
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
    igraph_t *graph, igraph_integer_t n, igraph_integer_t m,
    igraph_bool_t directed, igraph_bool_t loops, igraph_bool_t multiple
) {

    /* This function uses doubles in its `s` vector, and for `maxedges` and `last`.
     * This is because on a system with 32-bit ints, maxedges will be larger than
     * IGRAPH_INTEGER_MAX and this will cause overflows when calculating `from` and `to`
     * for tests on large graphs. This is also why we need a 'real' version of random_sample.
    */
    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_edges = m;
    igraph_real_t no_of_nodes_real = (igraph_real_t) no_of_nodes; /* for divisions below */
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int iter = 0;

    /* The multigraph implementation relies on the below checks to avoid overflow. */
    if (n < 0 || n > IGRAPH_VCOUNT_MAX) {
        IGRAPH_ERROR("Invalid number of vertices for G(n,m) model.", IGRAPH_EINVAL);
    }
    if (m < 0 || m > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Invalid number of edges for G(n,m) model..", IGRAPH_EINVAL);
    }

    if (no_of_edges == 0 || no_of_nodes == 0) {
        IGRAPH_CHECK(igraph_empty(graph, n, directed));
        return IGRAPH_SUCCESS;
    }

    if (multiple) {
        return gnm_multi(graph, n, m, directed, loops);
    }

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

    if (no_of_edges > maxedges) {
        IGRAPH_ERROR("Too many edges requested compared to the number of vertices.", IGRAPH_EINVAL);
    }

    if (maxedges == no_of_edges) {
        IGRAPH_CHECK(igraph_full(graph, n, directed, loops));
    } else {
        igraph_integer_t slen;

        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_random_sample_real(&s, 0, maxedges - 1, no_of_edges));

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, igraph_vector_size(&s) * 2));

        slen = igraph_vector_size(&s);
        if (directed && loops) {
            for (igraph_integer_t i = 0; i < slen; i++) {
                igraph_integer_t to = floor(VECTOR(s)[i] / no_of_nodes_real);
                igraph_integer_t from = VECTOR(s)[i] - to * no_of_nodes_real;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (directed && !loops) {
            for (igraph_integer_t i = 0; i < slen; i++) {
                igraph_integer_t from = floor(VECTOR(s)[i] / (no_of_nodes_real - 1));
                igraph_integer_t to = VECTOR(s)[i] - from * (no_of_nodes_real - 1);
                if (from == to) {
                    to = no_of_nodes - 1;
                }
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else if (!directed && loops) {
            for (igraph_integer_t i = 0; i < slen; i++) {
                igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                igraph_integer_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2;
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
                IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
            }
        } else { /* !directed && !loops */
            for (igraph_integer_t i = 0; i < slen; i++) {
                igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                igraph_integer_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2;
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
