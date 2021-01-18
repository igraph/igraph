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
 * Currently, this game simply uses \ref igraph_degree_sequence_game with 
 * the \c SIMPLE_NO_MULTIPLE method and appropriately constructed degree sequences.
 * Thefore, it does not sample uniformly: while it can generate all k-regular graphs
 * with the given number of vertices, it does not generate each one with the same
 * probability.
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
 * i.e. it generates each labelled tree with the same probability.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of nodes in the tree.
 * \param directed Whether to create a directed tree. The edges are oriented away from the root.
 * \param method The algorithm to use to generate the tree. Possible values:
 *        \clist
 *        \cli IGRAPH_RANDOM_TREE_PRUFER
 *          This algorithm samples Pr&uuml;fer sequences uniformly, then converts them to trees.
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
