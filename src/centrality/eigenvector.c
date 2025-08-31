/*
   igraph library.
   Copyright (C) 2007-2021  The igraph development team <igraph@igraph.org>

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

#include "igraph_centrality.h"

#include "igraph_adjlist.h"
#include "igraph_cycles.h"
#include "igraph_interface.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "centrality/centrality_internal.h"

#include <float.h>
#include <limits.h>
#include <math.h>

/* Multiplies vector 'from' by the unweighted adjacency matrix and stores the result in 'to'. */
static igraph_error_t adjmat_mul_unweighted(igraph_real_t *to, const igraph_real_t *from,
                                           int n, void *extra) {
    igraph_adjlist_t *adjlist = extra;
    igraph_vector_int_t *neis;
    igraph_int_t i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(adjlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_int_t nei = VECTOR(*neis)[j];
            to[i] += from[nei];
        }
    }

    return IGRAPH_SUCCESS;
}

typedef struct igraph_i_eigenvector_centrality_t {
    const igraph_t *graph;
    const igraph_inclist_t *inclist;
    const igraph_vector_t *weights;
} igraph_i_eigenvector_centrality_t;

/* Multiplies vector 'from' by the weighted adjacency matrix and stores the result in 'to'. */
static igraph_error_t adjmat_mul_weighted(igraph_real_t *to, const igraph_real_t *from,
                                            int n, void *extra) {

    igraph_i_eigenvector_centrality_t *data = extra;
    const igraph_t *graph = data->graph;
    const igraph_inclist_t *inclist = data->inclist;
    const igraph_vector_t *weights = data->weights;
    igraph_vector_int_t *edges;
    igraph_int_t i, j, nlen;

    for (i = 0; i < n; i++) {
        edges = igraph_inclist_get(inclist, i);
        nlen = igraph_vector_int_size(edges);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_int_t edge = VECTOR(*edges)[j];
            igraph_int_t nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            to[i] += w * from[nei];
        }
    }

    return IGRAPH_SUCCESS;
}

/* Checks if any eigenvector centrality values are zero. If they are, it indicates that the
 * graph is not (strongly) connected. Eigenvector centrality is not meaningful for such graphs.
 * To account for numerical inaccuracies, a threshold of 'eps' is used when testing for zero.
 * This function is intended to be used with eigenvector centrality values scaled such that
 * the maximum is 1. 'eps' is chosen accordinly.
 */
static void warn_zero_entries(const igraph_vector_t *evcent) {
    /* This is a conservative tolerance that will still catch most values
     * which should be zero without being triggered by small yet truly
     * nonzero values.
     * See https://github.com/igraph/igraph/pull/2592 */
    const igraph_real_t tol = 10 * DBL_EPSILON;
    const igraph_int_t n = igraph_vector_size(evcent);

    for (igraph_int_t i=0; i < n; i++) {
        igraph_real_t x = VECTOR(*evcent)[i];
        if (-tol < x && x < tol) {
            IGRAPH_WARNING("Some eigenvector centralities are nearly zero, indicating that the graph may not be (strongly) connected. "
                           "Eigenvector centrality is not meaningful for disconnected graphs.");
            return;
        }
    }
}

static igraph_error_t igraph_i_eigenvector_centrality_undirected(const igraph_t *graph, igraph_vector_t *vector,
                                                      igraph_real_t *value,
                                                      const igraph_vector_t *weights,
                                                      igraph_arpack_options_t *options) {

    igraph_vector_t values;
    igraph_matrix_t vectors;
    igraph_vector_t degree;
    igraph_int_t i;
    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t negative_weights = false;

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph has too many vertices for ARPACK.", IGRAPH_EOVERFLOW);
    }

    if (no_of_edges == 0) {
        /* special case: empty graph */
        if (value) {
            *value = 0;
        }
        if (vector) {
            IGRAPH_CHECK(igraph_vector_resize(vector, igraph_vcount(graph)));
            igraph_vector_fill(vector, 1);
        }
        if (no_of_nodes > 1) {
            IGRAPH_WARNING("The graph has no edges and is disconnected. "
                           "Eigenvector centrality is not meaningful for disconnected graphs.");
        }
        return IGRAPH_SUCCESS;
    }

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERRORF("Weights vector length (%" IGRAPH_PRId ") not equal to "
                    "number of edges (%" IGRAPH_PRId ").", IGRAPH_EINVAL,
                    igraph_vector_size(weights), no_of_edges);
        }
        /* Safe to call minmax, ecount == 0 case was caught earlier */
        igraph_vector_minmax(weights, &min, &max);
        if (min == 0 && max == 0) {
            /* special case: all weights are zeros */
            if (value) {
                *value = 0;
            }
            if (vector) {
                IGRAPH_CHECK(igraph_vector_resize(vector, igraph_vcount(graph)));
                igraph_vector_fill(vector, 1);
            }
            if (no_of_nodes > 1) {
                IGRAPH_WARNING("All edge weights are zero, making the graph effectively disconnected. "
                               "Eigenvector centrality is not meaningful for disconnected graphs.");
            }
            return IGRAPH_SUCCESS;
        }

        if (min < 0) {
            /* When there are negative weights, the eigenvalue and the eigenvector are no
             * longer guaranteed to be non-negative. */
            negative_weights = true;
            IGRAPH_WARNING("Negative weight in graph. The largest eigenvalue "
                           "will be selected, but it may not be the largest in magnitude. "
                           "Some eigenvector centralities may be negative.");
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, no_of_nodes, 1);

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(),
                                 IGRAPH_ALL, IGRAPH_LOOPS, weights));

    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(degree)[i]) {
            /* Note: Keep random perturbation non-negative. */
            MATRIX(vectors, i, 0) = VECTOR(degree)[i] + RNG_UNIF(0, 1e-4);
        } else if (! negative_weights) {
            /* The eigenvector centrality of zero degree vertices is also zero. */
            MATRIX(vectors, i, 0) = 0.0;
        } else {
            /* When negative weights are present, a zero strength may occur even
             * if the degree is not zero, and some edges have non-zero weight. */
            igraph_int_t deg;
            IGRAPH_CHECK(igraph_degree_1(graph, &deg, i, IGRAPH_ALL, IGRAPH_LOOPS));
            MATRIX(vectors, i, 0) = deg == 0 ? 0.0 : 1.0;
        }

    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    options->n = (int) no_of_nodes;
    options->nev = 1;
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->which[0] = 'L'; options->which[1] = 'A';
    options->start = 1;   /* no random start vector */

    if (!weights) {

        igraph_adjlist_t adjlist;

        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_CHECK(igraph_arpack_rssolve(adjmat_mul_unweighted,
                                           &adjlist, options, 0, &values, &vectors));

        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);

    } else {

        igraph_inclist_t inclist;
        igraph_i_eigenvector_centrality_t data;

        data.graph = graph;
        data.inclist = &inclist;
        data.weights = weights;

        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

        IGRAPH_CHECK(igraph_arpack_rssolve(adjmat_mul_weighted,
                                           &data, options, 0, &values, &vectors));

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (vector) {
        IGRAPH_CHECK(igraph_vector_resize(vector, no_of_nodes));

        /* Note: With non-negative weights, a zero eigenvalue should only occur
         * when there are no edges, or all edge have zero weight. This case is
         * caught earlier. Thus we do not handle the case of zero eigenvalues here .*/

        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
        }

        /* Scale result so that the largest value is 1.0. */
        igraph_i_vector_scale_by_max_abs(vector);

        /* Correction for numeric inaccuracies (eliminating -0.0) */
        if (! negative_weights) {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*vector)[i] < 0) {
                    VECTOR(*vector)[i] = 0;
                }
            }

            warn_zero_entries(vector);
        }
    }

    if (value) {
        *value = VECTOR(values)[0];
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine.");
    }

    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_eigenvector_centrality_directed(const igraph_t *graph, igraph_vector_t *vector,
                                                    igraph_real_t *value,
                                                    igraph_neimode_t mode,
                                                    const igraph_vector_t *weights,
                                                    igraph_arpack_options_t *options) {

    igraph_matrix_t values;
    igraph_matrix_t vectors;
    igraph_vector_t indegree;
    igraph_bool_t dag;
    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t negative_weights = false;

    mode = IGRAPH_REVERSE_MODE(mode);

    if (no_of_edges == 0) {
        /* special case: empty graph */
        if (value) {
            *value = 0;
        }
        if (vector) {
            IGRAPH_CHECK(igraph_vector_resize(vector, igraph_vcount(graph)));
            igraph_vector_fill(vector, 1);
        }
        if (no_of_nodes > 1) {
            IGRAPH_WARNING("The graph has no edges and is disconnected. "
                           "Eigenvector centrality is not meaningful for disconnected graphs.");
        }
        return IGRAPH_SUCCESS;
    }

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERRORF("Weights vector length (%" IGRAPH_PRId ") not equal to "
                    "number of edges (%" IGRAPH_PRId ").", IGRAPH_EINVAL,
                    igraph_vector_size(weights), no_of_edges);
        }

        /* Safe to call minmax, ecount == 0 case was caught earlier */
        igraph_vector_minmax(weights, &min, &max);

        if (min < 0.0) {
            /* When there are negative weights, the eigenvalue and the eigenvector are no
             * longer guaranteed to be non-negative, or even real-valued. */
            negative_weights = true;
            IGRAPH_WARNING("Negative weights in directed graph, eigenpair may be complex.");
        }
        if (min == 0.0 && max == 0.0) {
            /* special case: all weights are zeros */
            if (value) {
                *value = 0;
            }
            if (vector) {
                IGRAPH_CHECK(igraph_vector_resize(vector, igraph_vcount(graph)));
                igraph_vector_fill(vector, 1);
            }
            if (no_of_nodes > 1) {
                IGRAPH_WARNING("All edge weights are zero, making the graph effectively disconnected. "
                               "Eigenvector centrality is not meaningful for disconnected graphs.");
            }
            return IGRAPH_SUCCESS;
        }
    }

    /* If the graph is a DAG, the eigenvalue is zero, and degenerate. Vectors that
     * are zero in all vertices that have some out-edges, as well as their linear
     * combinations, are all valid eigenvectors.
     *
     * ARPACK may converge to one of these, or it may converge to an all-zero
     * vector, depending on chance. To eliminate this uncertainty, we return
     * an eigenvector that has 1s in sinks (vertices with no out-edges) and 0s
     * everywhere else. This ensures that the in-star, which is treated as the
     * "most centralized network" for eigenvector centralization, has a non-zero
     * centrality score in its centre. See https://github.com/igraph/igraph/issues/2679
     * for a discussion of the topic.
     *
     * For simplicity, we identify sinks by checking for zero out-strength. This
     * also identifies effective sinks where out-edges exist, but they all have zero
     * weight. This test is valid only when there are no negative weights, otherwise
     * non-zero weights may still add up to a zero strength. Since negative weights
     * are problematic anyway, we leave that case entirely up to ARPACK, without
     * special treatment.
     */
    if (! negative_weights) {
        IGRAPH_CHECK(igraph_is_dag(graph, &dag));
        if (dag) {
            /* special case: graph is a DAG */
            IGRAPH_WARNING("Graph is directed and acyclic; "
                           "returning eigenvector centralities of 1 in sink vertices, "
                           "and 0 everywhere else.");
            if (value) {
                *value = 0;
            }
            if (vector) {
                IGRAPH_CHECK(igraph_strength(
                    graph, vector, igraph_vss_all(), IGRAPH_REVERSE_MODE(mode),
                    IGRAPH_LOOPS, weights
                ));
                for (igraph_int_t i=0; i < no_of_nodes; i++) {
                    if (VECTOR(*vector)[i] == 0) {
                        VECTOR(*vector)[i] = 1;
                    } else {
                        VECTOR(*vector)[i] = 0;
                    }
                }
            }
            return IGRAPH_SUCCESS;
        }
    }

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph has too many vertices for ARPACK.", IGRAPH_EOVERFLOW);
    }

    options->n = (int) no_of_nodes;
    options->start = 1;
    options->nev = 1;
    /* Use higher NCV than usual as long directed cycles cause convergence issues.
     * According to numerical experiments, a cycle graph C_n tends to need about
     * NCV = n/4 at minimum. */
    options->ncv = no_of_nodes > 30 ? 30 : options->n;
    /* LM mode is not OK here because +1 and -1 can be eigenvalues at the
     * same time, e.g.: a -> b -> a, c -> a */
    options->which[0] = 'L' ; options->which[1] = 'R';

    IGRAPH_MATRIX_INIT_FINALLY(&values, 0, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, no_of_nodes, 1);

    IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &indegree, igraph_vss_all(),
                                 mode, IGRAPH_LOOPS, weights));

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        if (VECTOR(indegree)[i]) {
            /* Note: Keep random perturbation non-negative. */
            MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(0, 1e-4);
        } else if (! negative_weights) {
            /* The eigenvector centrality of zero in-degree vertices is also zero. */
            MATRIX(vectors, i, 0) = 0.0;
        } else {
            /* When negative weights are present, a zero in-strength may occur even
             * if the in-degree is not zero, and some in-edges have non-zero weight. */
            igraph_int_t deg;
            IGRAPH_CHECK(igraph_degree_1(graph, &deg, i, mode, IGRAPH_LOOPS));
            MATRIX(vectors, i, 0) = deg == 0 ? 0.0 : 1.0;
        }
    }

    igraph_vector_destroy(&indegree);
    IGRAPH_FINALLY_CLEAN(1);

    if (!weights) {
        igraph_adjlist_t adjlist;

        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_CHECK(igraph_arpack_rnsolve(adjmat_mul_unweighted,
                                           &adjlist, options, NULL, &values,
                                           &vectors));

        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_inclist_t inclist;
        igraph_i_eigenvector_centrality_t data;

        data.graph = graph;
        data.inclist = &inclist;
        data.weights = weights;

        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

        IGRAPH_CHECK(igraph_arpack_rnsolve(adjmat_mul_weighted,
                                           &data, options, NULL, &values, &vectors));

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (vector) {
        IGRAPH_CHECK(igraph_vector_resize(vector, options->n));

        if (!negative_weights && MATRIX(values, 0, 0) <= 0) {
            /* Pathological case: largest eigenvalue is zero, therefore all the
             * scores can also be zeros, this will be a valid eigenvector.
             * This usually happens with graphs that have lots of sinks and
             * sources only. */
            igraph_vector_fill(vector, 0);
            MATRIX(values, 0, 0) = 0;
        } else {
            for (igraph_int_t i = 0; i < no_of_nodes; i++) {
                VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
            }

            /* Scale result so that the largest value is 1.0. */
            igraph_i_vector_scale_by_max_abs(vector);
        }

        /* Correction for numeric inaccuracies (eliminating -0.0) */
        if (! negative_weights) {
            for (igraph_int_t i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*vector)[i] < 0) {
                    VECTOR(*vector)[i] = 0;
                }
            }

            warn_zero_entries(vector);
        }
    }

    if (value) {
        *value = MATRIX(values, 0, 0);
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine.");
    }

    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_eigenvector_centrality
 * \brief Eigenvector centrality of the vertices.
 *
 * Eigenvector centrality is a measure of the importance of a node in a
 * network. It assigns relative scores to all nodes in the network based
 * on the principle that connections from high-scoring nodes contribute
 * more to the score of the node in question than equal connections from
 * low-scoring nodes. Specifically, the eigenvector centrality of each
 * vertex is proportional to the sum of eigenvector centralities of its
 * neighbors. In practice, the centralities are determined by calculating the
 * eigenvector corresponding to the largest positive eigenvalue of the
 * adjacency matrix. This is motivated by the fact that the principal
 * eigenvector is guaranteed to be non-negative, assuming that edge weights
 * are also non-negative. In fact, in connected undirected graphs, this is
 * the \em only non-negative eigenvector.
 *
 * </para><para>
 * In the undirected case, this function considers
 * the diagonal entries of the adjacency matrix to be \em twice the number of
 * self-loops on the corresponding vertex.
 *
 * </para><para>
 * In the weighted case, the eigenvector centrality of a vertex is proportional
 * to the weighted sum of centralities of its neighbours, i.e.
 * <code>c_j = sum_i w_ij c_i</code>, where <code>w_ij</code> is the weight
 * of the edge connecting vertex \c i to \c j. The weights of parallel edges
 * are added up.
 *
 * </para><para>
 * The centrality scores returned by igraph are normalized such that the largest
 * eigenvector centrality score is 1, unless all scores are zeros.
 *
 * </para><para>
 * Eigenvector centrality is meaningful only for (strongly) connected graphs.
 * Undirected graphs that are not connected should be decomposed into connected
 * components, and the eigenvector centrality calculated for each separately.
 * This function does not directly verify that the graph is connected. If it is
 * not, in the undirected case the scores of all but one component will typically
 * be zeros. When zeros are detected, a warning is issued.
 *
 * </para><para>
 * Also note that the adjacency matrix of a directed acyclic graph or the
 * adjacency matrix of an empty graph does not possess positive eigenvalues,
 * therefore the eigenvector centrality is not meaningful for these graphs.
 * igraph will return an eigenvalue of zero in such cases. The returned
 * eigenvector centralities will all be equal for vertices with zero out-degree
 * or zero in-degrees (depending on whether \p mode is \c IGRAPH_OUT or \c IGRAPH_IN)
 * and zeros for other vertices. Such pathological cases can be detected
 * by asking igraph to calculate the eigenvalue as well (using the \p value
 * parameter, see below) and checking whether the eigenvalue is very close
 * to zero.
 *
 * </para><para>
 * Eigenvector centrality was developed for networks with non-negative edge
 * weights. While igraph does not refuse to carry out the calculation with
 * negative weights, it will issue a warning.
 *
 * </para><para>
 * When working with directed graphs, consider using hub and authority
 * scores instead, see \ref igraph_hub_and_authority_scores().
 *
 * \param graph The input graph. It may be directed.
 * \param vector Pointer to an initialized vector, it will be resized
 *     as needed. The result of the computation is stored here. It can
 *     be a null pointer, then it is ignored.
 * \param value If not a null pointer, then the eigenvalue
 *     corresponding to the found eigenvector is stored here.
 * \param mode How to consider edge directions in directed graphs.
 *     It is ignored for undirected graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the left eigenvector of the adjacency matrix is calculated,
 *          i.e. the centrality of a vertex is proportional to the sum
 *          of centralities of vertices pointing to it. This is the standard
 *          eigenvector centrality.
 *        \cli IGRAPH_IN
 *          the right eigenvector of the adjacency matrix is calculated,
 *          i.e. the centrality of a vertex is proportional to the sum
 *          of centralities of vertices it points to.
 *        \cli IGRAPH_ALL
 *          edge directions are ignored, and the unweighted eigenvector
 *          centrality is calculated.
 *        \endclist
 * \param weights A null pointer (indicating no edge weights), or a vector
 *     giving the weights of the edges. Weights should be positive to guarantee
 *     a meaningful result. The algorithm might produce complex numbers when some
 *     weights are negative and the graph is directed. In this case only
 *     the real part is reported.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Supply \c NULL here to use the defaults. Note that the
 *    function overwrites the <code>n</code> (number of vertices) parameter and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 *
 * Time complexity: depends on the input graph, usually it is O(|V|+|E|).
 *
 * \sa \ref igraph_pagerank and \ref igraph_personalized_pagerank for
 *   modifications of eigenvector centrality.
 * \ref igraph_hub_and_authority_scores() for a similar pair of measures
 * intended for directed graphs.
 *
 * \example examples/simple/eigenvector_centrality.c
 */

igraph_error_t igraph_eigenvector_centrality(const igraph_t *graph,
                                  igraph_vector_t *vector,
                                  igraph_real_t *value,
                                  igraph_neimode_t mode,
                                  const igraph_vector_t *weights,
                                  igraph_arpack_options_t *options) {

    if (!options) {
        options = igraph_arpack_options_get_default();
    }

    if (mode != IGRAPH_ALL && mode != IGRAPH_OUT && mode != IGRAPH_IN) {
        IGRAPH_ERROR("Invalid mode for eigenvector centrality.", IGRAPH_EINVAL);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (mode == IGRAPH_ALL) {
        return igraph_i_eigenvector_centrality_undirected(graph, vector, value,
                                                          weights, options);
    } else {
        return igraph_i_eigenvector_centrality_directed(graph, vector, value,
                                                        mode, weights, options);
    }
}
