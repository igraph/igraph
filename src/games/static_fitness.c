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

#include "igraph_adjlist.h"
#include "igraph_conversion.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_random.h"

#include "core/interruption.h"

/**
 * \ingroup generators
 * \function igraph_static_fitness_game
 * \brief Non-growing random graph with edge probabilities proportional to node fitness scores.
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
                               const igraph_vector_t *fitness_out, const igraph_vector_t *fitness_in,
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
        IGRAPH_ERROR("fitness_out must not be null.", IGRAPH_EINVAL);
    }

    if (no_of_edges < 0) {
        IGRAPH_ERRORF("Number of edges cannot be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, no_of_edges);
    }

    no_of_nodes = (igraph_integer_t) igraph_vector_size(fitness_out);
    if (no_of_nodes == 0) {
        IGRAPH_CHECK(igraph_empty(graph, 0, is_directed));
        return IGRAPH_SUCCESS;
    }

    if (is_directed && igraph_vector_size(fitness_in) != no_of_nodes) {
        IGRAPH_ERROR("fitness_in must have the same size as fitness_out.", IGRAPH_EINVAL);
    }

    /* Sanity checks for the fitnesses */
    if (igraph_vector_min(fitness_out) < 0) {
        IGRAPH_ERROR("Fitness scores must be non-negative.", IGRAPH_EINVAL);
    }
    if (fitness_in != 0 && igraph_vector_min(fitness_in) < 0) {
        IGRAPH_ERROR("Fitness scores must be non-negative.", IGRAPH_EINVAL);
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
            IGRAPH_ERROR("Too many edges requested.", IGRAPH_EINVAL);
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
        IGRAPH_ERRORF("Number of nodes cannot be negative, got %" IGRAPH_PRId".", IGRAPH_EINVAL, no_of_nodes);
    }

    /* Calculate alpha_out */
    if (exponent_out < 2) {
        IGRAPH_ERRORF("Out-degree exponent must be >= 2, got %g.", IGRAPH_EINVAL, exponent_out);
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
            IGRAPH_ERRORF("For directed graphs the in-degree exponent must be >= 2, got %g.",
                          IGRAPH_EINVAL, exponent_in);
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
