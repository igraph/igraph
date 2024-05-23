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
 * \param directed Logical, whether to generate a directed graph.
 * \param loops Logical, whether to generate self-loops.
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
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    }
    if (p < 0.0 || p > 1.0) {
        IGRAPH_ERROR("Invalid probability given.", IGRAPH_EINVAL);
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
 * \function igraph_erdos_renyi_game_gnm
 * \brief Generates a random (Erdős-Rényi) graph with a fixed number of edges.
 *
 * In the <code>G(n, m)</code> Erdős-Rényi model, a graph with \p n vertices
 * and \p m edges is generated uniformly at random.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param m The number of edges in the graph.
 * \param directed Logical, whether to generate a directed graph.
 * \param loops Logical, whether to generate self-loops.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid \p n or \p m parameter.
 *         \c IGRAPH_ENOMEM: there is not enough memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_erdos_renyi_game_gnp() to sample from the related
 * <code>G(n, p)</code> model, which constrains the \em expected edge count;
 * \ref igraph_degree_sequence_game() to constrain the degree sequence;
 * \ref igraph_bipartite_game_gnm() for the bipartite version of this model;
 * \ref igraph_barabasi_game() and \ref igraph_growing_random_game() for other
 * commonly used random graph models.
 *
 * \example examples/simple/igraph_erdos_renyi_game_gnm.c
 */
igraph_error_t igraph_erdos_renyi_game_gnm(
    igraph_t *graph, igraph_integer_t n, igraph_integer_t m,
    igraph_bool_t directed, igraph_bool_t loops
) {

    /* This function uses doubles in its `s` vector, and for `maxedges` and `last`.
     * This is because on a system with 32-bit ints, maxedges will be larger than
     * IGRAPH_INTEGER_MAX and this will cause overflows when calculating `from` and `to`
     * for tests on large graphs. This is also why we need a 'real' version of random_sample.
    */
    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_edges = m;
    igraph_real_t no_of_nodes_real = (igraph_real_t) no_of_nodes;   /* for divisions below */
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int iter = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    }
    if (m < 0 || m > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Invalid number of edges.", IGRAPH_EINVAL);
    }

    if (m == 0.0 || no_of_nodes == 0) {
        IGRAPH_CHECK(igraph_empty(graph, n, directed));
    } else {

        igraph_integer_t i;
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
                for (i = 0; i < slen; i++) {
                    igraph_integer_t to = floor(VECTOR(s)[i] / no_of_nodes_real);
                    igraph_integer_t from = VECTOR(s)[i] - to * no_of_nodes_real;
                    igraph_vector_int_push_back(&edges, from);
                    igraph_vector_int_push_back(&edges, to);
                    IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
                }
            } else if (directed && !loops) {
                for (i = 0; i < slen; i++) {
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
                for (i = 0; i < slen; i++) {
                    igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                    igraph_integer_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2;
                    igraph_vector_int_push_back(&edges, from);
                    igraph_vector_int_push_back(&edges, to);
                    IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
                }
            } else { /* !directed && !loops */
                for (i = 0; i < slen; i++) {
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
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_erdos_renyi_game
 * \brief Generates a random (Erdős-Rényi) graph.
 *
 * This function is deprecated; use \ref igraph_erdos_renyi_game_gnm() or
 * \ref igraph_erdos_renyi_game_gnp() instead.
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
 * \sa \ref igraph_barabasi_game(), \ref igraph_growing_random_game(),
 * \ref igraph_erdos_renyi_game_gnm(), \ref igraph_erdos_renyi_game_gnp()
 */
igraph_error_t igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
                            igraph_integer_t n, igraph_real_t p_or_m,
                            igraph_bool_t directed, igraph_bool_t loops) {

    if (type == IGRAPH_ERDOS_RENYI_GNP) {
        return igraph_erdos_renyi_game_gnp(graph, n, p_or_m, directed, loops);
    } else if (type == IGRAPH_ERDOS_RENYI_GNM) {
        return igraph_erdos_renyi_game_gnm(graph, n, (igraph_integer_t) p_or_m, directed, loops);
    } else {
        IGRAPH_ERROR("Invalid type", IGRAPH_EINVAL);
    }
}
