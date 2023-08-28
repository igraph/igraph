/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
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
#include "igraph_random.h"
#include "igraph_structural.h"

#include "centrality/prpack_internal.h"

#include <limits.h>

static igraph_error_t igraph_i_personalized_pagerank_arpack(const igraph_t *graph,
                                                 igraph_vector_t *vector,
                                                 igraph_real_t *value, const igraph_vs_t vids,
                                                 igraph_bool_t directed, igraph_real_t damping,
                                                 const igraph_vector_t *reset,
                                                 const igraph_vector_t *weights,
                                                 igraph_arpack_options_t *options);

typedef struct {
    const igraph_t *graph;
    igraph_adjlist_t *adjlist;
    igraph_real_t damping;
    igraph_vector_t *outdegree;
    igraph_vector_t *tmp;
    igraph_vector_t *reset;
} pagerank_data_t;

typedef struct {
    const igraph_t *graph;
    igraph_inclist_t *inclist;
    const igraph_vector_t *weights;
    igraph_real_t damping;
    igraph_vector_t *outdegree;
    igraph_vector_t *tmp;
    igraph_vector_t *reset;
} pagerank_data_weighted_t;

/* The two pagerank_operator functions below update the probabilities of a random walker
 * being in each of the vertices after one step of the walk. */

static igraph_error_t pagerank_operator_unweighted(igraph_real_t *to, const igraph_real_t *from,
                             int n, void *extra) {

    pagerank_data_t *data = extra;
    igraph_adjlist_t *adjlist = data->adjlist;
    igraph_vector_t *outdegree = data->outdegree;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_t *reset = data->reset;
    igraph_vector_int_t *neis;
    igraph_integer_t i, j, nlen;
    igraph_real_t sumfrom = 0.0;
    igraph_real_t fact = 1 - data->damping;

    /* Calculate p(x) / outdegree(x) in advance for all the vertices.
     * Note that we may divide by zero here; this is intentional since
     * we won't use those values and we save a comparison this way.
     * At the same time, we calculate the global probability of a
     * random jump in `sumfrom`. For vertices with no outgoing edges,
     * we will surely jump from there if we are there, hence those
     * vertices contribute p(x) to the teleportation probability.
     * For vertices with some outgoing edges, we jump from there with
     * probability `fact` if we are there, hence they contribute
     * p(x)*fact */
    for (i = 0; i < n; i++) {
        sumfrom += VECTOR(*outdegree)[i] != 0 ? from[i] * fact : from[i];
        VECTOR(*tmp)[i] = from[i] / VECTOR(*outdegree)[i];
    }

    /* Here we calculate the part of the `to` vector that results from
     * moving along links (and not from teleportation) */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(adjlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_integer_t nei = VECTOR(*neis)[j];
            to[i] += VECTOR(*tmp)[nei];
        }
        to[i] *= data->damping;
    }

    /* Now we add the contribution from random jumps. `reset` is a vector
     * that defines the probability of ending up in vertex i after a jump.
     * `sumfrom` is the global probability of jumping as mentioned above. */
    /* printf("sumfrom = %.6f\n", (float)sumfrom); */

    if (reset) {
        /* Running personalized PageRank */
        for (i = 0; i < n; i++) {
            to[i] += sumfrom * VECTOR(*reset)[i];
        }
    } else {
        /* Traditional PageRank with uniform reset vector */
        sumfrom /= n;
        for (i = 0; i < n; i++) {
            to[i] += sumfrom;
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t pagerank_operator_weighted(igraph_real_t *to, const igraph_real_t *from,
                              int n, void *extra) {

    pagerank_data_weighted_t *data = extra;
    const igraph_t *graph = data->graph;
    igraph_inclist_t *inclist = data->inclist;
    const igraph_vector_t *weights = data->weights;
    igraph_vector_t *outdegree = data->outdegree;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_t *reset = data->reset;
    igraph_integer_t i, j, nlen;
    igraph_real_t sumfrom = 0.0;
    igraph_vector_int_t *neis;
    igraph_real_t fact = 1 - data->damping;

    /*
    printf("PageRank weighted: multiplying vector: ");
    for (i=0; i<n; i++) { printf(" %.4f", from[i]); }
    printf("\n");
    */

    for (i = 0; i < n; i++) {
        if (VECTOR(*outdegree)[i] > 0) {
            sumfrom += from[i] * fact;
            VECTOR(*tmp)[i] = from[i] / VECTOR(*outdegree)[i];
        } else {
            sumfrom += from[i];
            /* The following value is used only when all outgoing edges have
             * weight zero (as opposed to there being no outgoing edges at all).
             * We set it to zero to avoid a 0.0*inf situation when computing
             * to[i] below. */
            VECTOR(*tmp)[i] = 0;
        }
    }

    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(inclist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t nei = IGRAPH_OTHER(graph, edge, i);
            to[i] += VECTOR(*weights)[edge] * VECTOR(*tmp)[nei];
        }
        to[i] *= data->damping;
    }

    /* printf("sumfrom = %.6f\n", (float)sumfrom); */

    if (reset) {
        /* Running personalized PageRank */
        for (i = 0; i < n; i++) {
            to[i] += sumfrom * VECTOR(*reset)[i];
        }
    } else {
        /* Traditional PageRank with uniform reset vector */
        sumfrom /= n;
        for (i = 0; i < n; i++) {
            to[i] += sumfrom;
        }
    }

    /*
    printf("PageRank weighted: multiplied vector: ");
    for (i=0; i<n; i++) { printf(" %.4f", to[i]); }
    printf("\n");
    */

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_pagerank
 * \brief Calculates the Google PageRank for the specified vertices.
 *
 * The PageRank centrality of a vertex is the fraction of time a
 * random walker traversing the graph would spend on that vertex.
 * The walker follows the out-edges with probabilities proportional
 * to their weights. Additionally, in each step, it restarts the walk
 * from a random vertex with probability <code>1 - damping</code>.
 * If the random walker gets stuck in a sink vertex, it will also restart
 * from a random vertex.
 *
 * </para><para>
 * The PageRank centrality is mainly useful for directed graphs. In undirected
 * graphs it converges to trivial values proportional to degrees as the damping
 * factor approaches 1.
 *
 * </para><para>
 * Starting from version 0.9, igraph has two PageRank implementations,
 * and the user can choose between them. The first implementation is
 * \c IGRAPH_PAGERANK_ALGO_ARPACK, which phrases the PageRank calculation
 * as an eigenvalue problem, which is then solved using the ARPACK library.
 * This was the default before igraph version 0.7. The second and recommended
 * implementation is \c IGRAPH_PAGERANK_ALGO_PRPACK. This is using the
 * PRPACK package, see https://github.com/dgleich/prpack. PRPACK uses an
 * algebraic method, i.e. solves a linear system to obtain the PageRank
 * scores.
 *
 * </para><para>
 * Note that the PageRank of a given vertex depends on the PageRank
 * of all other vertices, so even if you want to calculate the PageRank for
 * only some of the vertices, all of them must be calculated. Requesting
 * the PageRank for only some of the vertices does not result in any
 * performance increase at all.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual
 * Web Search Engine. Proceedings of the 7th World-Wide Web Conference,
 * Brisbane, Australia, April 1998.
 * https://doi.org/10.1016/S0169-7552(98)00110-X
 *
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_ARPACK, \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable. When using \c IGRAPH_PAGERANK_ALGO_ARPACK,
 *    the eigenvalue corresponding to the PageRank vector is stored here. It is
 *    expected to be exactly one. Checking this value can be used to diagnose cases
 *    when ARPACK failed to converge to the leading eigenvector.
 *    When using \c IGRAPH_PAGERANK_ALGO_PRPACK, this is always set to 1.0.
 * \param vids The vertex IDs for which the PageRank is returned. This parameter
 *    is only for convenience. Computing PageRank for fewer than all vertices will
 *    not speed up the calculation.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper).
 *    Must be a probability in the range [0, 1]. A commonly used value is 0.85.
 * \param weights Optional edge weights. May be a \c NULL pointer,
 *    meaning unweighted edges, or a vector of non-negative values
 *    of the same length as the number of edges.
 * \param options Options for the ARPACK method. See \ref igraph_arpack_options_t
 *    for details. Supply \c NULL here to use the defaults. Note that the function
 *    overwrites the <code>n</code> (number of vertices), <code>nev</code> (1),
 *    <code>ncv</code> (3) and <code>which</code> (LM) parameters and it always
 *    starts the calculation from a non-random vector calculated based on the
 *    degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for temporary data.
 *         \c IGRAPH_EINVVID, invalid vertex ID in \p vids.
 *
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 *
 * \sa \ref igraph_personalized_pagerank() and \ref igraph_personalized_pagerank_vs()
 * for the personalized PageRank measure. See \ref igraph_arpack_rssolve() and
 * \ref igraph_arpack_rnsolve() for the underlying machinery used by
 * \c IGRAPH_PAGERANK_ALGO_ARPACK.
 *
 * \example examples/simple/igraph_pagerank.c
 */

igraph_error_t igraph_pagerank(const igraph_t *graph, igraph_pagerank_algo_t algo,
                    igraph_vector_t *vector,
                    igraph_real_t *value, const igraph_vs_t vids,
                    igraph_bool_t directed, igraph_real_t damping,
                    const igraph_vector_t *weights, igraph_arpack_options_t *options) {
    return igraph_personalized_pagerank(graph, algo, vector, value, vids,
                                        directed, damping, NULL, weights,
                                        options);
}

/**
 * \function igraph_personalized_pagerank_vs
 * \brief Calculates the personalized Google PageRank for the specified vertices.
 *
 * The personalized PageRank is similar to the original PageRank measure, but
 * when the random walk is restarted, a new starting vertex is chosen according to
 * a specified distribution.
 * This distribution is used both when restarting randomly with probability
 * <code>1 - damping</code>, and when the walker is forced to restart due to being
 * stuck in a sink vertex (a vertex with no outgoing edges).
 *
 * </para><para>
 * This simplified interface takes a vertex sequence and resets the random walk to
 * one of the vertices in the specified vertex sequence, chosen uniformly. A typical
 * application of personalized PageRank is when the random walk is reset to the same
 * vertex every time - this can easily be achieved using \ref igraph_vss_1() which
 * generates a vertex sequence containing only a single vertex.
 *
 * </para><para>
 * Note that the personalized PageRank of a given vertex depends on the
 * personalized PageRank of all other vertices, so even if you want to calculate
 * the personalized PageRank for only some of the vertices, all of them must be
 * calculated. Requesting the personalized PageRank for only some of the vertices
 * does not result in any performance increase at all.
 *
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_ARPACK, \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable. When using \c IGRAPH_PAGERANK_ALGO_ARPACK,
 *    the eigenvalue corresponding to the PageRank vector is stored here. It is
 *    expected to be exactly one. Checking this value can be used to diagnose cases
 *    when ARPACK failed to converge to the leading eigenvector.
 *    When using \c IGRAPH_PAGERANK_ALGO_PRPACK, this is always set to 1.0.
 * \param vids The vertex IDs for which the PageRank is returned. This parameter
 *    is only for convenience. Computing PageRank for fewer than all vertices will
 *    not speed up the calculation.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper).
 *    Must be a probability in the range [0, 1]. A commonly used value is 0.85.
 * \param reset_vids IDs of the vertices used when resetting the random walk.
 * \param weights Optional edge weights, it is either a null pointer,
 *    then the edges are not weighted, or a vector of the same length
 *    as the number of edges.
 * \param options Options for the ARPACK method. See \ref igraph_arpack_options_t
 *    for details. Supply \c NULL here to use the defaults. Note that the function
 *    overwrites the <code>n</code> (number of vertices), <code>nev</code> (1),
 *    <code>ncv</code> (3) and <code>which</code> (LM) parameters and it always
 *    starts the calculation from a non-random vector calculated based on the
 *    degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *         \c IGRAPH_EINVVID, invalid vertex ID in
 *         \p vids or an empty reset vertex sequence in
 *         \p vids_reset.
 *
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 *
 * \sa \ref igraph_pagerank() for the non-personalized implementation.
 */

igraph_error_t igraph_personalized_pagerank_vs(const igraph_t *graph,
                                    igraph_pagerank_algo_t algo, igraph_vector_t *vector,
                                    igraph_real_t *value, const igraph_vs_t vids,
                                    igraph_bool_t directed, igraph_real_t damping,
                                    igraph_vs_t reset_vids,
                                    const igraph_vector_t *weights,
                                    igraph_arpack_options_t *options) {
    igraph_vector_t reset;
    igraph_vit_t vit;

    IGRAPH_VECTOR_INIT_FINALLY(&reset, igraph_vcount(graph));
    IGRAPH_CHECK(igraph_vit_create(graph, reset_vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    while (!IGRAPH_VIT_END(vit)) {
        VECTOR(reset)[IGRAPH_VIT_GET(vit)]++;
        IGRAPH_VIT_NEXT(vit);
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_personalized_pagerank(graph, algo, vector,
                 value, vids, directed,
                 damping, &reset, weights,
                 options));

    igraph_vector_destroy(&reset);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_personalized_pagerank
 * \brief Calculates the personalized Google PageRank for the specified vertices.
 *
 * The personalized PageRank is similar to the original PageRank measure, but
 * when the random walk is restarted, a new starting vertex is chosen non-uniformly,
 * according to the distribution specified in \p reset
 * (instead of the uniform distribution in the original PageRank measure).
 * The \p reset distribution is used both when restarting randomly with probability
 * <code>1 - damping</code>, and when the walker is forced to restart due to being
 * stuck in a sink vertex (a vertex with no outgoing edges).
 *
 * </para><para>
 * Note that the personalized PageRank of a given vertex depends on the
 * personalized PageRank of all other vertices, so even if you want to calculate
 * the personalized PageRank for only some of the vertices, all of them must be
 * calculated. Requesting the personalized PageRank for only some of the vertices
 * does not result in any performance increase at all.
 *
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_ARPACK, \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable. When using \c IGRAPH_PAGERANK_ALGO_ARPACK,
 *    the eigenvalue corresponding to the PageRank vector is stored here. It is
 *    expected to be exactly one. Checking this value can be used to diagnose cases
 *    when ARPACK failed to converge to the leading eigenvector.
 *    When using \c IGRAPH_PAGERANK_ALGO_PRPACK, this is always set to 1.0.
 * \param vids The vertex IDs for which the PageRank is returned. This parameter
 *    is only for convenience. Computing PageRank for fewer than all vertices will
 *    not speed up the calculation.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper).
 *    Must be a probability in the range [0, 1]. A commonly used value is 0.85.
 * \param reset The probability distribution over the vertices used when
 *    resetting the random walk. It is either a \c NULL pointer (denoting
 *    a uniform choice that results in the original PageRank measure)
 *    or a vector of the same length as the number of vertices.
 * \param weights Optional edge weights. May be a \c NULL pointer,
 *    meaning unweighted edges, or a vector of non-negative values
 *    of the same length as the number of edges.
 * \param options Options for the ARPACK method. See \ref igraph_arpack_options_t
 *    for details. Supply \c NULL here to use the defaults. Note that the function
 *    overwrites the <code>n</code> (number of vertices), <code>nev</code> (1),
 *    <code>ncv</code> (3) and <code>which</code> (LM) parameters and it always
 *    starts the calculation from a non-random vector calculated based on the
 *    degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *         \c IGRAPH_EINVVID, invalid vertex ID in
 *         \p vids or an invalid reset vector in \p reset.
 *
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 *
 * \sa \ref igraph_pagerank() for the non-personalized implementation,
 * \ref igraph_personalized_pagerank_vs() for a personalized implementation
 * with resetting to specific vertices.
 */
igraph_error_t igraph_personalized_pagerank(const igraph_t *graph,
                                 igraph_pagerank_algo_t algo, igraph_vector_t *vector,
                                 igraph_real_t *value, const igraph_vs_t vids,
                                 igraph_bool_t directed, igraph_real_t damping,
                                 const igraph_vector_t *reset,
                                 const igraph_vector_t *weights,
                                 igraph_arpack_options_t *options) {

    if (damping < 0.0 || damping > 1.0) {
        IGRAPH_ERROR("The PageRank damping factor must be in the range [0,1].", IGRAPH_EINVAL);
    }

    if (algo == IGRAPH_PAGERANK_ALGO_ARPACK) {
        return igraph_i_personalized_pagerank_arpack(graph, vector, value, vids,
                directed, damping, reset,
                weights, options ? options : igraph_arpack_options_get_default()
        );
    } else if (algo == IGRAPH_PAGERANK_ALGO_PRPACK) {
        return igraph_i_personalized_pagerank_prpack(graph, vector, value, vids,
                directed, damping, reset,
                weights);
    }

    IGRAPH_ERROR("Unknown PageRank algorithm", IGRAPH_EINVAL);
}

/*
 * ARPACK-based implementation of \c igraph_personalized_pagerank.
 *
 * See \c igraph_personalized_pagerank for the documentation of the parameters.
 */
static igraph_error_t igraph_i_personalized_pagerank_arpack(const igraph_t *graph, igraph_vector_t *vector,
                                                 igraph_real_t *value, const igraph_vs_t vids,
                                                 igraph_bool_t directed, igraph_real_t damping,
                                                 const igraph_vector_t *reset,
                                                 const igraph_vector_t *weights,
                                                 igraph_arpack_options_t *options) {
    igraph_matrix_t values;
    igraph_matrix_t vectors;
    igraph_neimode_t dirmode;
    igraph_vector_t outdegree;
    igraph_vector_t indegree;
    igraph_vector_t tmp;
    igraph_vector_t normalized_reset;

    igraph_integer_t i;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);

    igraph_real_t reset_sum; /* used only when reset != NULL */

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph has too many vertices for ARPACK.", IGRAPH_EOVERFLOW);
    }

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid length of weights vector when calculating PageRank scores.", IGRAPH_EINVAL);
    }

    if (reset && igraph_vector_size(reset) != no_of_nodes) {
        IGRAPH_ERROR("Invalid length of reset vector when calculating personalized PageRank scores.", IGRAPH_EINVAL);
    }

    if (reset) {
        reset_sum = igraph_vector_sum(reset);
        if (no_of_nodes > 0 && reset_sum == 0) {
            IGRAPH_ERROR("The sum of the elements in the reset vector must not be zero.", IGRAPH_EINVAL);
        }

        igraph_real_t reset_min = igraph_vector_min(reset);
        if (reset_min < 0) {
            IGRAPH_ERROR("The reset vector must not contain negative elements.", IGRAPH_EINVAL);
        }
        if (isnan(reset_min)) {
            IGRAPH_ERROR("The reset vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    if (no_of_edges == 0) {
        /* Special case: graph with no edges. Result is the same as the personalization vector. */
        if (value) {
            *value = 1.0;
        }
        if (vector) {
            if (reset && no_of_nodes > 0) {
                IGRAPH_CHECK(igraph_vector_update(vector, reset));
                igraph_vector_scale(vector, 1.0 / reset_sum);
            } else {
                IGRAPH_CHECK(igraph_vector_resize(vector, no_of_nodes));
                igraph_vector_fill(vector, 1.0 / no_of_nodes);
            }
        }
        return IGRAPH_SUCCESS;
    }

    options->n = (int) no_of_nodes;
    options->nev = 1;
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rnsolve */
    options->which[0] = 'L'; options->which[1] = 'R';
    options->start = 1; /* no random start vector */

    directed = directed && igraph_is_directed(graph);

    if (weights) {
        igraph_real_t min, max;

        /* Safe to call minmax, ecount == 0 case was caught earlier */
        igraph_vector_minmax(weights, &min, &max);
        if (min < 0) {
            IGRAPH_ERROR("Edge weights must not be negative.", IGRAPH_EINVAL);
        }
        if (isnan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
        if (min == 0 && max == 0) {
            /* Special case: all weights are zeros. Result is the same as the personalization vector. */
            if (value) {
                *value = 1.0;
            }
            if (vector) {
                IGRAPH_CHECK(igraph_vector_resize(vector, no_of_nodes));
                if (reset) {
                    for (i=0; i < no_of_nodes; ++i) {
                        VECTOR(*vector)[i] = VECTOR(*reset)[i];
                    }
                    igraph_vector_scale(vector, 1.0 / igraph_vector_sum(vector));
                } else {
                    igraph_vector_fill(vector, 1.0 / no_of_nodes);
                }
            }
            return IGRAPH_SUCCESS;
        }
    }

    IGRAPH_MATRIX_INIT_FINALLY(&values, 0, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);

    if (directed) {
        dirmode = IGRAPH_IN;
    } else {
        dirmode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&indegree, options->n);
    IGRAPH_VECTOR_INIT_FINALLY(&outdegree, options->n);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);

    RNG_BEGIN();

    if (reset) {
        /* Normalize reset vector so the sum is 1 */
        IGRAPH_CHECK(igraph_vector_init_copy(&normalized_reset, reset));
        IGRAPH_FINALLY(igraph_vector_destroy, &normalized_reset);

        igraph_vector_scale(&normalized_reset, 1.0 / reset_sum);
    }

    IGRAPH_CHECK(igraph_strength(graph, &outdegree, igraph_vss_all(),
                                 directed ? IGRAPH_OUT : IGRAPH_ALL, IGRAPH_LOOPS, weights));
    IGRAPH_CHECK(igraph_strength(graph, &indegree, igraph_vss_all(),
                                 directed ? IGRAPH_IN : IGRAPH_ALL, IGRAPH_LOOPS, weights));

    /* Set up an appropriate starting vector. We start from the (possibly weight) in-degrees
     * plus some small random noise to avoid convergence problems. */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(indegree)[i] > 0) {
            MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(-1e-4, 1e-4);
        } else {
            MATRIX(vectors, i, 0) = 1;
        }
    }

    if (!weights) {

        igraph_adjlist_t adjlist;
        pagerank_data_t data;

        data.graph = graph;
        data.adjlist = &adjlist;
        data.damping = damping;
        data.outdegree = &outdegree;
        data.tmp = &tmp;
        data.reset = reset ? &normalized_reset : NULL;

        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, dirmode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_CHECK(igraph_arpack_rnsolve(pagerank_operator_unweighted,
                                           &data, options, NULL, &values, &vectors));

        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);

    } else {

        igraph_inclist_t inclist;
        pagerank_data_weighted_t data;

        data.graph = graph;
        data.inclist = &inclist;
        data.weights = weights;
        data.damping = damping;
        data.outdegree = &outdegree;
        data.tmp = &tmp;
        data.reset = reset ? &normalized_reset : NULL;

        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, dirmode, IGRAPH_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

        IGRAPH_CHECK(igraph_arpack_rnsolve(pagerank_operator_weighted,
                                           &data, options, NULL, &values, &vectors));

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    RNG_END();

    if (reset) {
        igraph_vector_destroy(&normalized_reset);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&outdegree);
    igraph_vector_destroy(&indegree);
    IGRAPH_FINALLY_CLEAN(3);

    if (value) {
        *value = MATRIX(values, 0, 0);
    }

    if (vector) {
        igraph_vit_t vit;
        igraph_integer_t nodes_to_calc;
        igraph_real_t sum = 0;

        for (i = 0; i < no_of_nodes; i++) {
            sum += MATRIX(vectors, i, 0);
        }

        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        nodes_to_calc = IGRAPH_VIT_SIZE(vit);

        IGRAPH_CHECK(igraph_vector_resize(vector, nodes_to_calc));
        for (IGRAPH_VIT_RESET(vit), i = 0; !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), i++) {
            VECTOR(*vector)[i] = MATRIX(vectors, IGRAPH_VIT_GET(vit), 0);
            VECTOR(*vector)[i] /= sum;
        }

        igraph_vit_destroy(&vit);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
    }

    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
