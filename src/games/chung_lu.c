/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_games.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_random.h"
#include "igraph_vector.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

#include <math.h> /* isfinite() */

/* This implementation follows the ideas in:
 *
 * Joel Miller and Aric Hagberg,
 * Efficient generation of networks with given expected degrees
 * (2011)
 *
 * It is analogous to the method used in igraph_erdos_renyi_game_gnp()
 * and has linear complexity in the number of edges.
 */

static igraph_error_t check_expected_degrees(const igraph_vector_t *weights) {
    igraph_real_t minw, maxw;

    igraph_vector_minmax(weights, &minw, &maxw);

    if (minw < 0) {
        IGRAPH_ERRORF("Vertex weights must not be negative in Chung-Lu model, got %g.", IGRAPH_EINVAL, minw);
    }

    /* Catches both NaN and +Inf. */
    if (! isfinite(maxw)) {
        IGRAPH_ERRORF("Vertex weights must be finite, got %g.", IGRAPH_EINVAL, maxw);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_chung_lu_game
 * \brief Samples graphs from the Chung-Lu model.
 *
 * \experimental
 *
 * The Chung-Lu model is useful for generating random graphs with fixed
 * expected degrees. This function implements both the original model of Chung
 * and Lu, as well as some additional variants with useful properties.
 *
 * </para><para>
 * In the original Chung-Lu model, each pair of vertices \c i and \c j is
 * connected with independent probability <code>p_ij = w_i w_j / S</code>,
 * where \c w_i is a weight associated with vertex \c i and
 * <code>S = sum_k w_k</code> is the sum of weights. In the directed variant,
 * vertices have both out-weights, <code>w^out</code>, and in-weights,
 * <code>w^in</code>, with equal sums,
 * <code>S = sum_k w^out_k = sum_k w^in_k</code>.
 * The connection probability between \c i and \c j is
 * <code>p_ij = w^out_i w^in_j / S</code>.
 *
 * </para><para>
 * This model is commonly used to create random graphs with a fixed \em expected
 * degree sequence. The expected degree of vertex \c i is approximately equal
 * to the weight \c w_i. Specifically, if the graph is directed and self-loops
 * are allowed, then the expected out- and in-degrees are precisely
 * <code>w^out</code> and <code>w^in</code>. If self-loops are disallowed,
 * then the expected out- and in-degrees are <code>w^out (S - w^in) / S</code>
 * and <code>w^in (S - w^out) / S</code>, respectively. If the graph is
 * undirected, then the expected degrees with and without self-loops are
 * <code>w (S + w) / S</code> and <code>w (S - w) / S</code>, respectively.
 *
 * </para><para>
 * A limitation of the original Chung-Lu model is that when some of the
 * weights are large, the formula for \c p_ij yields values larger than 1.
 * Chung and Lu's original paper excludes the use of such weights. When
 * <code>p_ij > 1</code>, this function simply issues a warning and creates
 * a connection between \c i and \c j. However, in this case the expected degrees
 * will no longer relate to the weights in the manner stated above. Thus the
 * original Chung-Lu model cannot produce certain (large) expected degrees.
 *
 * </para><para>
 * The overcome this limitation, this function implements additional variants of
 * the model, with modified expressions for the connection probability \c p_ij
 * between vertices \c i and \c j. Let <code>q_ij = w_i w_j / S</code>, or
 * <code>q_ij = w^out_i w^in_j / S</code> in the directed case. All model
 * variants become equivalent in the limit of sparse graphs where \c q_ij
 * approaches zero. In the original Chung-Lu model, selectable by setting
 * \p variant to \c IGRAPH_CHUNG_LU_ORIGINAL, <code>p_ij = min(q_ij, 1)</code>.
 * The \c IGRAPH_CHUNG_LU_MAXENT variant, sometiems referred to a the generalized
 * random graph, uses <code>p_ij = q_ij / (1 + q_ij)</code>, and is equivalent
 * to a maximum entropy model (i.e. exponential random graph model) with
 * a constraint on expected degrees; see Park and Newman (2004), Section B,
 * setting <code>exp(-Theta_ij) = w_i w_j / S</code>. This model is also
 * discussed by Britton, Deijfen and Martin-Löf (2006). By virtue of being
 * a degree-constrained maximum entropy model, it produces graphs with the
 * same degree sequence with the same probability.
 * A third variant can be requested with \c IGRAPH_CHUNG_LU_NR, and uses
 * <code>p_ij = 1 - exp(-q_ij)</code>. This is the underlying simple graph
 * of a multigraph model introduced by Norros and Reittu (2006).
 * For a discussion of these three model variants, see Section 16.4 of
 * Bollobás, Janson, Riordan (2007), as well as Van Der Hofstad (2013).
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Chung F and Lu L: Connected components in a random graph with given
 * degree sequences. Annals of Combinatorics 6, 125-145 (2002).
 * https://doi.org/10.1007/PL00012580
 *
 * </para><para>
 * Miller JC and Hagberg A:
 * Efficient Generation of Networks with Given Expected Degrees (2011).
 * https://doi.org/10.1007/978-3-642-21286-4_10
 *
 * </para><para>
 * Park J and Newman MEJ: Statistical mechanics of networks.
 * Physical Review E 70, 066117 (2004).
 * https://doi.org/10.1103/PhysRevE.70.066117
 *
 * </para><para>
 * Britton T, Deijfen M, Martin-Löf A:
 * Generating Simple Random Graphs with Prescribed Degree Distribution.
 * J Stat Phys 124, 1377–1397 (2006).
 * https://doi.org/10.1007/s10955-006-9168-x
 *
 * </para><para>
 * Norros I and Reittu H: On a conditionally Poissonian graph process.
 * Advances in Applied Probability 38, 59–75 (2006).
 * https://doi.org/10.1239/aap/1143936140
 *
 * </para><para>
 * Bollobás B, Janson S, Riordan O:
 * The phase transition in inhomogeneous random graphs.
 * Random Struct Algorithms 31, 3–122 (2007).
 * https://doi.org/10.1002/rsa.20168
 *
 * </para><para>
 * Van Der Hofstad R: Critical behavior in inhomogeneous random graphs.
 * Random Struct Algorithms 42, 480–508 (2013).
 * https://doi.org/10.1002/rsa.20450
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param out_weights A vector of non-negative vertex weights (or out-weights).
 *    In sparse graphs these will be approximately equal to the expected
 *    (out-)degrees.
 * \param in_weights A vector of non-negative in-weights, approximately equal
 *    to the expected in-degrees in sparse graphs. May be set to \c NULL,
 *    in which case undirected graphs are generated.
 * \param loops Whether to allow the creation of self-loops. Since vertex
 *    pairs are connected independently, setting this to false is equivalent
 *    to simply discarding self-loops from an existing loopy Chung-Lu graph.
 * \param variant The model variant to sample from, with different definitions
 *    of the connection probability between vertices \c i and \c j. Given
 *    <code>q_ij = w_i w_j / S</code>, the following formulations are available:
 *    \clist
 *    \cli IGRAPH_CHUNG_LU_ORIGINAL
 *         the original Chung-Lu model, <code>p_ij = min(q_ij, 1)</code>.
 *    \cli IGRAPH_CHUNG_LU_MAXENT
 *         maximum entropy model with fixed expected degrees,
 *         <code>p_ij = q_ij / (1 + q_ij)</code>.
 *    \cli IGRAPH_CHUNG_LU_NR
 *         Norros and Reittu's model, <code>p_ij = 1 - exp(-q_ij)</code>.
 *    \endclist
 * \return Error code.
 *
 * \sa \ref igraph_static_fitness_game() implements a similar model with
 * a sharp constraint on the number of edges;
 * \ref igraph_degree_sequence_game() samples random graphs with sharply
 * specified degrees; \ref igraph_erdos_renyi_game_gnp() creates random
 * graphs with a fixed connection probability \c p between all vertex pairs.
 *
 * Time complexity: O(|E| + |V|), linear in the number of edges.
 */
igraph_error_t igraph_chung_lu_game(igraph_t *graph,
                                    const igraph_vector_t *out_weights,
                                    const igraph_vector_t *in_weights,
                                    igraph_bool_t loops,
                                    igraph_chung_lu_t variant) {

    const igraph_int_t no_of_nodes = igraph_vector_size(out_weights);
    const igraph_bool_t directed = in_weights != NULL;
    igraph_vector_int_t edges, idx;
    igraph_real_t wsum = igraph_vector_sum(out_weights);
    igraph_bool_t warned = false;
    int iter = 0;

    /* Necessitated by floating point arithmetic used in the implementation. */
    if (no_of_nodes >= IGRAPH_MAX_EXACT_REAL) {
        IGRAPH_ERROR("Number of vertices is too large.", IGRAPH_EOVERFLOW);
    }

    if (! directed) {
        in_weights = out_weights;
    } else if (igraph_vector_size(in_weights) != no_of_nodes) {
        IGRAPH_ERROR("Vertex out- and in-weight vectors must have the same length.", IGRAPH_EINVAL);
    }

    if (no_of_nodes == 0) {
        return igraph_empty(graph, 0, directed);
    }

    IGRAPH_CHECK(check_expected_degrees(out_weights));

    if (directed) {
        IGRAPH_CHECK(check_expected_degrees(in_weights));

        if (igraph_vector_sum(in_weights) != wsum) {
            IGRAPH_ERRORF("Sum of out- and in-weights must be the same, got %g and %g, respectively.",
                          IGRAPH_EINVAL, wsum, igraph_vector_sum(in_weights));
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&idx, 0);

    IGRAPH_CHECK(igraph_vector_sort_ind(in_weights, &idx, IGRAPH_DESCENDING));

    for (igraph_int_t i=0; i < no_of_nodes; i++) {
        igraph_int_t vi, vj;
        igraph_real_t wi, wj;
        igraph_real_t p, q;

        igraph_int_t j = directed ? 0 : i;

        vi = VECTOR(idx)[i];
        wi = VECTOR(*out_weights)[vi];

        if (wi == 0) {
            if (directed) continue; else break;
        }

        p = 1;

        while (true) {
            igraph_real_t gap = RNG_GEOM(p);

            /* This formulation not only terminates the loop when necessary,
             * but also protects against overflow when 'p' is very small
             * and 'gap' becomes very large, perhaps larger than representable
             * in an igraph_int_t. */
            if (gap >= no_of_nodes-j) {
                break;
            }

            j += gap;

            vj = VECTOR(idx)[j];
            wj = VECTOR(*in_weights)[vj];

            q = wi * wj / wsum;
            switch (variant) {
            case IGRAPH_CHUNG_LU_ORIGINAL:
                if (q > 1) {
                    q = 1;
                    if (! warned && (loops || vi != vj)) {
                        IGRAPH_WARNINGF(
                            "Expected degrees %g and %g lead to a calculated connection probability "
                            "larger than 1 in Chung-Lu model. The degrees of the resulting graph will "
                            "not be consistent with the given input.", wi, wj);
                        warned = true;
                    }
                }
                break;
            case IGRAPH_CHUNG_LU_MAXENT:
                q = q / (1 + q);
                break;
            case IGRAPH_CHUNG_LU_NR:
                q = 1 - exp(-q);
                break;
            default:
                IGRAPH_ERROR("Invalid Chung-Lu variant.", IGRAPH_EINVAL);
            }

            /* A probability of zero must not be passed to RNG_GEOM(),
             * so we catch this case here. */
            if (q == 0) {
                break;
            }

            if (RNG_UNIF01() < q/p && (loops || vi != vj)) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, vi));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, vj));
            }

            p = q;

            j++;

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 16);
        }
    }

    igraph_vector_int_destroy(&idx);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
