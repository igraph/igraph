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
#include "igraph_interface.h"
#include "igraph_structural.h"
#include "igraph_random.h"

#include "centrality/centrality_internal.h"

#include <float.h>
#include <limits.h>

/* struct for the unweighted variant of the HITS algorithm */
typedef struct igraph_i_kleinberg_data_t {
    igraph_adjlist_t *in;
    igraph_adjlist_t *out;
    igraph_vector_t *tmp;
} igraph_i_kleinberg_data_t;

/* struct for the weighted variant of the HITS algorithm */
typedef struct igraph_i_kleinberg_data2_t {
    const igraph_t *graph;
    igraph_inclist_t *in;
    igraph_inclist_t *out;
    igraph_vector_t *tmp;
    const igraph_vector_t *weights;
} igraph_i_kleinberg_data2_t;

/* Checks if at least a certain fraction of HITS centrality scores are zero.
 * Any zero value indicates that the graphs corresponding to A A^T and A^T A are
 * not connected, and therefore the solution is not unique. However, this situation
 * is fairly common, and difficult to control. Thefore we only warn if the number
 * of zero values exceeds a certain fraction.
 *
 * To account for numerical inaccuracies, a threshold of 'eps' is used when testing for zero.
 * This function is intended to be used with centrality values scaled such that
 * the maximum is 1. 'eps' is chosen accordinly.
 *
 * See the analogous function used in igraph_eigenvector_centrality() for details
 * on the choice of 'eps'.
 */
static void warn_zero_entries(const igraph_vector_t *cent) {
    const igraph_real_t tol = 10 * DBL_EPSILON;
    const igraph_real_t frac = 0.3; /* warn if at least this fraction of centralities is zero */
    const igraph_int_t n = igraph_vector_size(cent);

    /* Skip check for small graphs */
    if (n < 10) {
        return;
    }

    const igraph_int_t max_zero_cnt = (igraph_int_t) frac*n;
    igraph_int_t zero_cnt = 0;

    for (igraph_int_t i=0; i < n; i++) {
        igraph_real_t x = VECTOR(*cent)[i];
        if (-tol < x && x < tol) {
            if (++zero_cnt > max_zero_cnt) {
                IGRAPH_WARNINGF(
                    "More than %d%% of hub or authority scores are zeros. The presence of zero values "
                    "indicates that the solution is not unique, thus the returned result may not be meaningful.",
                    (int) (frac * 100)
                );
                return;
            }
        }
    }
}

static igraph_error_t igraph_i_kleinberg_unweighted_hub_to_auth(
        igraph_int_t n, igraph_vector_t *to, const igraph_real_t *from,
        igraph_adjlist_t *in) {
    igraph_vector_int_t *neis;
    igraph_int_t i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(in, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*to)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_int_t nei = VECTOR(*neis)[j];
            VECTOR(*to)[i] += from[nei];
        }
    }
    return IGRAPH_SUCCESS;
}

/* ARPACK auxiliary routine for the unweighted HITS algorithm */
static igraph_error_t igraph_i_kleinberg_unweighted(igraph_real_t *to,
                                         const igraph_real_t *from,
                                         int n, void *extra) {
    igraph_i_kleinberg_data_t *data = (igraph_i_kleinberg_data_t*)extra;
    igraph_adjlist_t *out = data->out;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    igraph_int_t i, j, nlen;

    igraph_i_kleinberg_unweighted_hub_to_auth(n, tmp, from, data->in);

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(out, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_int_t nei = VECTOR(*neis)[j];
            to[i] += VECTOR(*tmp)[nei];
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_kleinberg_weighted_hub_to_auth(igraph_int_t n,
        igraph_vector_t *to, const igraph_real_t *from, igraph_inclist_t *in,
        const igraph_t *g, const igraph_vector_t *weights) {
    igraph_vector_int_t *neis;
    igraph_int_t nlen, i, j;
    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(in, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*to)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_int_t nei_edge = VECTOR(*neis)[j];
            igraph_int_t nei = IGRAPH_OTHER(g, nei_edge, i);
            VECTOR(*to)[i] += from[nei] * VECTOR(*weights)[nei_edge];
        }
    }
    return IGRAPH_SUCCESS;
}

/* ARPACK auxiliary routine for the weighted HITS algorithm */
static igraph_error_t igraph_i_kleinberg_weighted(igraph_real_t *to,
                                       const igraph_real_t *from,
                                       int n, void *extra) {

    igraph_i_kleinberg_data2_t *data = (igraph_i_kleinberg_data2_t*)extra;
    igraph_inclist_t *out = data->out;
    igraph_vector_t *tmp = data->tmp;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *g = data->graph;
    igraph_vector_int_t *neis;
    igraph_int_t i, j, nlen;

    igraph_i_kleinberg_weighted_hub_to_auth(n, tmp, from, data->in, g, weights);

    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(out, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_int_t nei_edge = VECTOR(*neis)[j];
            igraph_int_t nei = IGRAPH_OTHER(g, nei_edge, i);
            to[i] += VECTOR(*tmp)[nei] * VECTOR(*weights)[nei_edge];
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hub_and_authority_scores
 * \brief Kleinberg's hub and authority scores (HITS).
 *
 * Hub and authority scores are a generalization of the ideas behind
 * eigenvector centrality to directed graphs. The authority score of
 * a vertex is proportional to the sum of the hub scores of vertices
 * that point to it. Conversely, the hub score of a vertex is proportional
 * to the sum of authority scores of vertices that it points to. These
 * concepts are also known under the name Hyperlink-Induced Topic Search (HITS).
 *
 * </para><para>
 * The hub and authority scores of the vertices are defined as the principal
 * eigenvectors of <code>A A^T</code> and <code>A^T A</code>, respectively,
 * where <code>A</code> is the adjacency matrix of the graph and <code>A^T</code>
 * is its transpose. The motivation for choosing the principal eigenvector
 * is that it is guaranteed to be non-negative when edge weights are also
 * non-negative.
 *
 * </para><para>
 * If vectors \c h and \c a contain hub and authority scores, then the two
 * scores are related by <code>h = A a</code> and <code>a = A^T h</code>.
 * When the principal eigenvalue of <code>A A^T</code> is degenerate, there
 * is no unique solution to the hub- and authority-score problem.
 * igraph guarantees that the scores that are returned are matching, i.e. are
 * related by these formulas, even in this situation.
 *
 * </para><para>
 * Note that hub and authority scores are not well behaved in extremely sparse
 * graphs where no single connected component dominates the undirected graphs
 * corresponding to <code>A A^T</code> and <code>A^T A</code>. In these cases,
 * there are many different non-negative eigenvectors, all reasonable solutions
 * to the HITS equations. The symptom of such a situation is that a large
 * fraction of the scores are zeros. igraph issues a warning when this is
 * detected.
 *
 * </para><para>
 * Results are scaled so that the largest hub and authority scores are both 1.
 *
 * </para><para>
 * The concept of hub and authority scores were developed for \em directed graphs.
 * In undirected graphs, both the hub and authority scores are equal to the
 * eigenvector centrality, which can be computed using
 * \ref igraph_eigenvector_centrality().
 *
 * </para><para>
 * HITS scores were developed for networks with non-negative edge weights.
 * While igraph does not refuse to carry out the calculation with negative
 * weights, it will issue a warning.
 *
 * </para><para>
 * See the following reference on the meaning of this score:
 * J. Kleinberg. Authoritative sources in a hyperlinked
 * environment. \emb Proc. 9th ACM-SIAM Symposium on Discrete
 * Algorithms, \eme 1998. Extended version in \emb Journal of the
 * ACM \eme 46 (1999).
 * https://doi.org/10.1145/324133.324140
 * Also appears as IBM Research Report RJ 10076, May 1997.
 *
 * \param graph The input graph. Can be directed and undirected.
 * \param hub_vector Pointer to an initialized vector, the hub scores are
 *    stored here. If a null pointer then it is ignored.
 * \param authority_vector Pointer to an initialized vector, the authority
 *    scores are stored here. If a null pointer then it is ignored.
 * \param value If not a null pointer then the eigenvalue
 *    corresponding to the calculated eigenvectors is stored here.
 * \param weights A null pointer (meaning no edge weights), or a vector
 *     giving the weights of the edges.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Supply \c NULL here to use the defaults. Note that the
 *    function overwrites the <code>n</code> (number of vertices) parameter
 *    and it always starts the calculation from a vector calculated based on
 *    the degree of the vertices.
 * \return Error code.
 *
 * Time complexity: depends on the input graph, usually it is O(|V|),
 * the number of vertices.
 *
 * \sa \ref igraph_pagerank(), \ref igraph_personalized_pagerank();
 * \ref igraph_eigenvector_centrality() for a similar measure intended
 * for undirected graphs.
 */
igraph_error_t igraph_hub_and_authority_scores(const igraph_t *graph,
        igraph_vector_t *hub_vector, igraph_vector_t *authority_vector,
        igraph_real_t *value,
        const igraph_vector_t *weights, igraph_arpack_options_t *options) {

    /* The current implementation computes hub scores, i.e the principal
     * eigenvector of A A^T, and transforms these to authority scores as
     * authority = A^T hub. */

    igraph_adjlist_t inadjlist, outadjlist;
    igraph_inclist_t ininclist, outinclist;
    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_t tmp;
    igraph_vector_t values;
    igraph_matrix_t vectors;
    igraph_i_kleinberg_data_t extra;
    igraph_i_kleinberg_data2_t extra2;
    igraph_vector_t *my_hub_vector_p;
    igraph_vector_t my_hub_vector;
    igraph_bool_t negative_weights = false;

    if (! igraph_is_directed(graph)) {
        /* In undirected graphs, hub and authority scores are the same as eigenvector
         * centralities. We issue a warning to avoid user confusion.
         * If Ax = lambda x then A^2 x = lambda^2 x. Therefore the principal
         * eigenvector of A is also a principal eigenvector of A^2. However,
         * if both lambda and -lambda are eigenvalues of A, then the lambda^2
         * eigenvalue of A^2 will be degenerate. This happens for example in
         * an even cycle graph where 2 and -2 are the largest eigenvalues in magnitude,
         * therefore 4 is a degenerate eigenvalue of A^2, with eigenspace spanned by
         * (0, 1, 0, 1, ...) and (1, 0, 1, ...). The vector (1, 1, ...), which
         * would be expected by users based on symmetry considerations, may not be
         * returned. We avoid such issues by falling back to igraph_eigenvector_centrality() */

        IGRAPH_WARNING("Hub and authority scores requested for undirected graph. "
                       "These are the same as eigenvector centralities.");
        if (! hub_vector) {
            IGRAPH_VECTOR_INIT_FINALLY(&my_hub_vector, no_of_nodes);
            my_hub_vector_p = &my_hub_vector;
        } else {
            my_hub_vector_p = hub_vector;
        }
        IGRAPH_CHECK(igraph_eigenvector_centrality(graph, my_hub_vector_p, value, IGRAPH_ALL, weights, options));
        *value = (*value) * (*value); /* adjust the eigenvalue, see comment at top */
        if (authority_vector) {
            IGRAPH_CHECK(igraph_vector_update(authority_vector, my_hub_vector_p));
        }
        if (! hub_vector) {
            igraph_vector_destroy(&my_hub_vector);
            IGRAPH_FINALLY_CLEAN(1);
        }
        return IGRAPH_SUCCESS;
    }

    if (igraph_ecount(graph) == 0) {
        /* special case: empty graph */
        if (value) {
            *value = 0;
        }
        if (hub_vector) {
            IGRAPH_CHECK(igraph_vector_resize(hub_vector, no_of_nodes));
            igraph_vector_fill(hub_vector, 1.0);
        }
        if (authority_vector) {
            IGRAPH_CHECK(igraph_vector_resize(authority_vector, no_of_nodes));
            igraph_vector_fill(authority_vector, 1.0);
        }
        if (no_of_nodes > 1) {
            IGRAPH_WARNING("The graph has no edges. Hub and authority scores are not meaningful.");
        }
        return IGRAPH_SUCCESS;
    }

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != igraph_ecount(graph)) {
            IGRAPH_ERRORF(
                    "Weights vector length (%" IGRAPH_PRId ") should match number of "
                    "edges (%" IGRAPH_PRId ") when calculating "
                    "hub or authority scores.",
                    IGRAPH_EINVAL,
                    igraph_vector_size(weights),
                    igraph_ecount(graph));
        }

        /* Safe to call minmax, ecount == 0 case was caught earlier */
        igraph_vector_minmax(weights, &min, &max);

        if (min < 0.0) {
            /* When there are negative weights, the principal eigenvalue and the eigenvector
             * are no longer guaranteed to be non-negative. */
            negative_weights = true;
            IGRAPH_WARNING("Negative weight in graph. The largest eigenvalue "
                           "will be selected, but it may not be the largest in magnitude. "
                           "Some hub and authority scores may be negative.");
        }

        if (min == 0 && max == 0) {
            /* special case: all weights are zeros */
            if (value) {
                *value = 0;
            }
            if (hub_vector) {
                IGRAPH_CHECK(igraph_vector_resize(hub_vector, no_of_nodes));
                igraph_vector_fill(hub_vector, 1);
            }
            if (authority_vector) {
                IGRAPH_CHECK(igraph_vector_resize(authority_vector, no_of_nodes));
                igraph_vector_fill(authority_vector, 1);
            }
            IGRAPH_WARNING("All edge weights are zero. Hub and authority scores are not meaningful.");
            return IGRAPH_SUCCESS;
        }
    }

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph has too many vertices for ARPACK.", IGRAPH_EOVERFLOW);
    }

    if (!options) {
        options = igraph_arpack_options_get_default();
    }

    options->n = (int) no_of_nodes;
    options->start = 1;   /* no random start vector */

    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);

    if (weights == NULL) {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &inadjlist, IGRAPH_IN, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &inadjlist);
        IGRAPH_CHECK(igraph_adjlist_init(graph, &outadjlist, IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &outadjlist);
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &ininclist, IGRAPH_IN, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &ininclist);
        IGRAPH_CHECK(igraph_inclist_init(graph, &outinclist, IGRAPH_OUT, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &outinclist);
    }

    /* We calculate hub scores, which correlate with out-degrees / out-strengths.
     * Thus we use out-strengths as starting values. */
    IGRAPH_CHECK(igraph_strength(graph, &tmp, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS, weights));

    for (igraph_int_t i = 0; i < options->n; i++) {
        if (VECTOR(tmp)[i] != 0) {
            /* Note: Keep random perturbation non-negative. */
            MATRIX(vectors, i, 0) = VECTOR(tmp)[i] + RNG_UNIF(0, 1e-4);
        } else if (! negative_weights) {
            /* The hub score of zero out-degree vertices is also zero. */
            MATRIX(vectors, i, 0) = 0.0;
        } else {
            /* When negative weights are present, a zero out-strength may occur even
             * if the out-degree is not zero, and some out-edges have non-zero weight. */
            igraph_int_t deg;
            IGRAPH_CHECK(igraph_degree_1(graph, &deg, i, IGRAPH_OUT, /* loops */ true));
            MATRIX(vectors, i, 0) = deg == 0 ? 0.0 : 1.0;
        }
    }

    extra.in = &inadjlist; extra.out = &outadjlist; extra.tmp = &tmp;
    extra2.in = &ininclist; extra2.out = &outinclist; extra2.tmp = &tmp;
    extra2.graph = graph; extra2.weights = weights;

    options->nev = 1;
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->which[0] = 'L'; options->which[1] = 'A';

    if (weights == NULL) {
        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_kleinberg_unweighted, &extra,
                                           options, 0, &values, &vectors));
    } else {
        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_kleinberg_weighted, &extra2,
                                           options, 0, &values, &vectors));
    }

    if (value) {
        *value = VECTOR(values)[0];
    }

    if (hub_vector || authority_vector) {
        if (!hub_vector) {
            IGRAPH_VECTOR_INIT_FINALLY(&my_hub_vector, options->n);
            my_hub_vector_p = &my_hub_vector;
        } else {
            my_hub_vector_p = hub_vector;
        }

        IGRAPH_CHECK(igraph_vector_resize(my_hub_vector_p, options->n));
        for (igraph_int_t i = 0; i < options->n; i++) {
            VECTOR(*my_hub_vector_p)[i] = MATRIX(vectors, i, 0);
        }

        igraph_i_vector_scale_by_max_abs(my_hub_vector_p);

        /* Correction for numeric inaccuracies (eliminating -0.0) */
        if (! negative_weights) {
            for (igraph_int_t i = 0; i < options->n; i++) {
                if (VECTOR(*my_hub_vector_p)[i] < 0) {
                    VECTOR(*my_hub_vector_p)[i] = 0;
                }
            }
        }

        warn_zero_entries(my_hub_vector_p);
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
    }
    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    if (authority_vector) {
        IGRAPH_CHECK(igraph_vector_resize(authority_vector, no_of_nodes));
        igraph_vector_null(authority_vector);
        if (weights == NULL) {
            igraph_i_kleinberg_unweighted_hub_to_auth(no_of_nodes, authority_vector, &VECTOR(*my_hub_vector_p)[0], &inadjlist);
        } else {
            igraph_i_kleinberg_weighted_hub_to_auth(no_of_nodes, authority_vector, &VECTOR(*my_hub_vector_p)[0], &ininclist, graph, weights);
        }
        igraph_i_vector_scale_by_max_abs(authority_vector);
    }

    if (!hub_vector && authority_vector) {
        igraph_vector_destroy(&my_hub_vector);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (weights == NULL) {
        igraph_adjlist_destroy(&outadjlist);
        igraph_adjlist_destroy(&inadjlist);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        igraph_inclist_destroy(&outinclist);
        igraph_inclist_destroy(&ininclist);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
