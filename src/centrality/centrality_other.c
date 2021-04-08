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

#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_structural.h"
#include "igraph_topology.h"
#include "igraph_stack.h"
#include "igraph_dqueue.h"

#include "centrality/prpack_internal.h"
#include "core/indheap.h"
#include "core/interruption.h"
#include "core/math.h"

#include "config.h"

#include <math.h>
#include <string.h>    /* memset */

static int igraph_i_personalized_pagerank_arpack(const igraph_t *graph,
                                                 igraph_vector_t *vector,
                                                 igraph_real_t *value, const igraph_vs_t vids,
                                                 igraph_bool_t directed, igraph_real_t damping,
                                                 const igraph_vector_t *reset,
                                                 const igraph_vector_t *weights,
                                                 igraph_arpack_options_t *options);

static igraph_bool_t igraph_i_vector_mostly_negative(const igraph_vector_t *vector) {
    /* Many of the centrality measures correspond to the eigenvector of some
     * matrix. When v is an eigenvector, c*v is also an eigenvector, therefore
     * it may happen that all the scores in the eigenvector are negative, in which
     * case we want to negate them since the centrality scores should be positive.
     * However, since ARPACK is not always stable, sometimes it happens that
     * *some* of the centrality scores are small negative numbers. This function
     * helps distinguish between the two cases; it should return true if most of
     * the values are relatively large negative numbers, in which case we should
     * negate the eigenvector.
     */
    long int n = igraph_vector_size(vector);
    igraph_real_t mi, ma;

    if (n == 0) {
        return 0;
    }

    igraph_vector_minmax(vector, &mi, &ma);

    if (mi >= 0) {
        return 0;
    }
    if (ma <= 0) {
        return 1;
    }

    /* is the most negative value larger in magnitude than the most positive? */
    return (-mi/ma > 1);
}

static int igraph_i_eigenvector_centrality(igraph_real_t *to, const igraph_real_t *from,
                                           int n, void *extra) {
    igraph_adjlist_t *adjlist = extra;
    igraph_vector_int_t *neis;
    long int i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(adjlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            to[i] += from[nei];
        }
    }


    return 0;
}

typedef struct igraph_i_eigenvector_centrality_t {
    const igraph_t *graph;
    const igraph_inclist_t *inclist;
    const igraph_vector_t *weights;
} igraph_i_eigenvector_centrality_t;

static int igraph_i_eigenvector_centrality2(igraph_real_t *to, const igraph_real_t *from,
                                            int n, void *extra) {

    igraph_i_eigenvector_centrality_t *data = extra;
    const igraph_t *graph = data->graph;
    const igraph_inclist_t *inclist = data->inclist;
    const igraph_vector_t *weights = data->weights;
    igraph_vector_int_t *edges;
    long int i, j, nlen;

    for (i = 0; i < n; i++) {
        edges = igraph_inclist_get(inclist, i);
        nlen = igraph_vector_int_size(edges);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*edges)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            to[i] += w * from[nei];
        }
    }

    return IGRAPH_SUCCESS;
}

static int igraph_i_eigenvector_centrality_undirected(const igraph_t *graph, igraph_vector_t *vector,
                                                      igraph_real_t *value, igraph_bool_t scale,
                                                      const igraph_vector_t *weights,
                                                      igraph_arpack_options_t *options) {

    igraph_vector_t values;
    igraph_matrix_t vectors;
    igraph_vector_t degree;
    long int i;

    options->n = igraph_vcount(graph);
    options->start = 1;   /* no random start vector */

    if (igraph_ecount(graph) == 0) {
        /* special case: empty graph */
        if (value) {
            *value = 0;
        }
        if (vector) {
            igraph_vector_resize(vector, igraph_vcount(graph));
            igraph_vector_fill(vector, 1);
        }
        return IGRAPH_SUCCESS;
    }

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != igraph_ecount(graph)) {
            IGRAPH_ERROR("Invalid length of weights vector when calculating "
                         "eigenvector centrality", IGRAPH_EINVAL);
        }
        /* Safe to call minmax, ecount == 0 case was caught earlier */
        IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));
        if (min == 0 && max == 0) {
            /* special case: all weights are zeros */
            if (value) {
                *value = 0;
            }
            if (vector) {
                igraph_vector_resize(vector, igraph_vcount(graph));
                igraph_vector_fill(vector, 1);
            }
            return IGRAPH_SUCCESS;
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);

    IGRAPH_VECTOR_INIT_FINALLY(&degree, options->n);
    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                               IGRAPH_ALL, /*loops=*/ 0));
    RNG_BEGIN();
    for (i = 0; i < options->n; i++) {
        if (VECTOR(degree)[i]) {
            MATRIX(vectors, i, 0) = VECTOR(degree)[i] + RNG_UNIF(-1e-4, 1e-4);
        } else {
            MATRIX(vectors, i, 0) = 1.0;
        }
    }
    RNG_END();
    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    options->n = igraph_vcount(graph);
    options->nev = 1;
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->which[0] = 'L'; options->which[1] = 'A';
    options->start = 1;   /* no random start vector */

    if (!weights) {

        igraph_adjlist_t adjlist;

        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigenvector_centrality,
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

        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigenvector_centrality2,
                                           &data, options, 0, &values, &vectors));

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (value) {
        *value = VECTOR(values)[0];
    }

    if (vector) {
        igraph_real_t amax = 0;
        long int which = 0;
        long int i;
        IGRAPH_CHECK(igraph_vector_resize(vector, options->n));

        if (VECTOR(values)[0] <= 0) {
            /* Pathological case: largest eigenvalue is zero, therefore all the
             * scores can also be zeros, this will be a valid eigenvector.
             * This usually happens with graphs that have lots of sinks and
             * sources only. */
            igraph_vector_fill(vector, 0);
        } else {
            for (i = 0; i < options->n; i++) {
                igraph_real_t tmp;
                VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
                tmp = fabs(VECTOR(*vector)[i]);
                if (tmp > amax) {
                    amax = tmp;
                    which = i;
                }
            }
            if (scale && amax != 0) {
                igraph_vector_scale(vector, 1 / VECTOR(*vector)[which]);
            } else if (igraph_i_vector_mostly_negative(vector)) {
                igraph_vector_scale(vector, -1.0);
            }

            /* Correction for numeric inaccuracies (eliminating -0.0) */
            for (i = 0; i < options->n; i++) {
                if (VECTOR(*vector)[i] < 0) {
                    VECTOR(*vector)[i] = 0;
                }
            }
        }
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
    }

    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/* int igraph_i_evcent_dir(igraph_real_t *to, const igraph_real_t *from, */
/*          long int n, void *extra) { */
/*   /\* TODO *\/ */
/*   return 0; */
/* } */

/* int igraph_i_evcent_dir2(igraph_real_t *to, const igraph_real_t *from, */
/*           long int n, void *extra) { */
/*   /\* TODO *\/ */
/*   return 0; */
/* } */

static int igraph_i_eigenvector_centrality_directed(const igraph_t *graph, igraph_vector_t *vector,
                                                    igraph_real_t *value, igraph_bool_t scale,
                                                    const igraph_vector_t *weights,
                                                    igraph_arpack_options_t *options) {

    igraph_matrix_t values;
    igraph_matrix_t vectors;
    igraph_vector_t indegree;
    igraph_bool_t dag;
    long int i;

    if (igraph_ecount(graph) == 0) {
        /* special case: empty graph */
        if (value) {
            *value = 0;
        }
        if (vector) {
            igraph_vector_resize(vector, igraph_vcount(graph));
            igraph_vector_fill(vector, 1);
        }
        return IGRAPH_SUCCESS;
    }

    /* Quick check: if the graph is a DAG, all the eigenvector centralities are
     * zeros, and so is the eigenvalue */
    IGRAPH_CHECK(igraph_is_dag(graph, &dag));
    if (dag) {
        /* special case: graph is a DAG */
        IGRAPH_WARNING("graph is directed and acyclic; eigenvector centralities "
                       "will be zeros");
        if (value) {
            *value = 0;
        }
        if (vector) {
            igraph_vector_resize(vector, igraph_vcount(graph));
            igraph_vector_fill(vector, 0);
        }
        return IGRAPH_SUCCESS;
    }

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != igraph_ecount(graph)) {
            IGRAPH_ERROR("Invalid length of weights vector when calculating "
                         "eigenvector centrality", IGRAPH_EINVAL);
        }
        if (igraph_is_directed(graph)) {
            IGRAPH_WARNING("Weighted directed graph in eigenvector centrality");
        }

        /* Safe to call minmax, ecount == 0 case was caught earlier */
        IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));

        if (min < 0.0) {
            IGRAPH_WARNING("Negative weights, eigenpair might be complex");
        }
        if (min == 0.0 && max == 0.0) {
            /* special case: all weights are zeros */
            if (value) {
                *value = 0;
            }
            if (vector) {
                igraph_vector_resize(vector, igraph_vcount(graph));
                igraph_vector_fill(vector, 1);
            }
            return IGRAPH_SUCCESS;
        }
    }

    options->n = igraph_vcount(graph);
    options->start = 1;
    options->nev = 1;
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rnsolve */
    /* LM mode is not OK here because +1 and -1 can be eigenvalues at the
     * same time, e.g.: a -> b -> a, c -> a */
    options->which[0] = 'L' ; options->which[1] = 'R';

    IGRAPH_MATRIX_INIT_FINALLY(&values, 0, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);

    IGRAPH_VECTOR_INIT_FINALLY(&indegree, options->n);
    IGRAPH_CHECK(igraph_strength(graph, &indegree, igraph_vss_all(),
                                 IGRAPH_IN, /*loops=*/ 1, weights));
    RNG_BEGIN();
    for (i = 0; i < options->n; i++) {
        if (VECTOR(indegree)[i]) {
            MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(-1e-4, 1e-4);
        } else {
            MATRIX(vectors, i, 0) = 1.0;
        }
    }
    RNG_END();
    igraph_vector_destroy(&indegree);
    IGRAPH_FINALLY_CLEAN(1);

    if (!weights) {
        igraph_adjlist_t adjlist;

        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_eigenvector_centrality,
                                           &adjlist, options, 0, &values,
                                           &vectors));

        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_inclist_t inclist;
        igraph_i_eigenvector_centrality_t data;

        data.graph = graph;
        data.inclist = &inclist;
        data.weights = weights;

        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_IN, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

        IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_eigenvector_centrality2,
                                           &data, options, 0, &values, &vectors));

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (value) {
        *value = MATRIX(values, 0, 0);
    }

    if (vector) {
        igraph_real_t amax = 0;
        long int which = 0;
        long int i;
        IGRAPH_CHECK(igraph_vector_resize(vector, options->n));

        if (MATRIX(values, 0, 0) <= 0) {
            /* Pathological case: largest eigenvalue is zero, therefore all the
             * scores can also be zeros, this will be a valid eigenvector.
             * This usually happens with graphs that have lots of sinks and
             * sources only. */
            igraph_vector_fill(vector, 0);
            MATRIX(values, 0, 0) = 0;
        } else {
            for (i = 0; i < options->n; i++) {
                igraph_real_t tmp;
                VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
                tmp = fabs(VECTOR(*vector)[i]);
                if (tmp > amax) {
                    amax = tmp;
                    which = i;
                }
            }
            if (scale && amax != 0) {
                igraph_vector_scale(vector, 1 / VECTOR(*vector)[which]);
            } else if (igraph_i_vector_mostly_negative(vector)) {
                igraph_vector_scale(vector, -1.0);
            }
        }

        /* Correction for numeric inaccuracies (eliminating -0.0) */
        for (i = 0; i < options->n; i++) {
            if (VECTOR(*vector)[i] < 0) {
                VECTOR(*vector)[i] = 0;
            }
        }
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
    }

    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_eigenvector_centrality
 * Eigenvector centrality of the vertices
 *
 * Eigenvector centrality is a measure of the importance of a node in a
 * network. It assigns relative scores to all nodes in the network based
 * on the principle that connections from high-scoring nodes contribute
 * more to the score of the node in question than equal connections from
 * low-scoring nodes. Specifically, the eigenvector centrality of each
 * vertex is proportional to the sum of eigenvector centralities of its
 * neighbors. In practice, the centralities are determined by calculating the
 * eigenvector corresponding to the largest positive eigenvalue of the
 * adjacency matrix. In the undirected case, this function considers
 * the diagonal entries of the adjacency matrix to be \em twice the number of
 * self-loops on the corresponding vertex.
 *
 * </para><para>
 * The centrality scores returned by igraph can be normalized
 * (using the \p scale parameter) such that the largest eigenvector centrality
 * score is 1 (with one exception, see below).
 *
 * </para><para>
 * In the directed case, the left eigenvector of the adjacency matrix is
 * calculated. In other words, the centrality of a vertex is proportional
 * to the sum of centralities of vertices pointing to it.
 *
 * </para><para>
 * Eigenvector centrality is meaningful only for connected graphs.
 * Graphs that are not connected should be decomposed into connected
 * components, and the eigenvector centrality calculated for each separately.
 * This function does not verify that the graph is connected. If it is not,
 * in the undirected case the scores of all but one component will be zeros.
 *
 * </para><para>
 * Also note that the adjacency matrix of a directed acyclic graph or the
 * adjacency matrix of an empty graph does not possess positive eigenvalues,
 * therefore the eigenvector centrality is not defined for these graphs.
 * igraph will return an eigenvalue of zero in such cases. The eigenvector
 * centralities will all be equal for an empty graph and will all be zeros
 * for a directed acyclic graph. Such pathological cases can be detected
 * by asking igraph to calculate the eigenvalue as well (using the \p value
 * parameter, see below) and checking whether the eigenvalue is very close
 * to zero.
 *
 * \param graph The input graph. It may be directed.
 * \param vector Pointer to an initialized vector, it will be resized
 *     as needed. The result of the computation is stored here. It can
 *     be a null pointer, then it is ignored.
 * \param value If not a null pointer, then the eigenvalue
 *     corresponding to the found eigenvector is stored here.
 * \param directed Boolean scalar, whether to consider edge directions
 *     in a directed graph. It is ignored for undirected graphs.
 * \param scale If not zero then the result will be scaled such that
 *     the absolute value of the maximum centrality is one.
 * \param weights A null pointer (=no edge weights), or a vector
 *     giving the weights of the edges. The algorithm might result
 *     complex numbers is some weights are negative. In this case only
 *     the real part is reported.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 *
 * Time complexity: depends on the input graph, usually it is O(|V|+|E|).
 *
 * \sa \ref igraph_pagerank and \ref igraph_personalized_pagerank for
 *   modifications of eigenvector centrality.
 *
 * \example examples/simple/eigenvector_centrality.c
 */

int igraph_eigenvector_centrality(const igraph_t *graph,
                                  igraph_vector_t *vector,
                                  igraph_real_t *value,
                                  igraph_bool_t directed, igraph_bool_t scale,
                                  const igraph_vector_t *weights,
                                  igraph_arpack_options_t *options) {

    if (directed && igraph_is_directed(graph)) {
        return igraph_i_eigenvector_centrality_directed(graph, vector, value,
                scale, weights, options);
    } else {
        return igraph_i_eigenvector_centrality_undirected(graph, vector, value,
                scale, weights, options);
    }
}

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

/* ARPACK auxiliary routine for the unweighted HITS algorithm */
static int igraph_i_kleinberg_unweighted(igraph_real_t *to,
                                         const igraph_real_t *from,
                                         int n, void *extra) {
    igraph_i_kleinberg_data_t *data = (igraph_i_kleinberg_data_t*)extra;
    igraph_adjlist_t *in = data->in;
    igraph_adjlist_t *out = data->out;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    long int i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(in, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            VECTOR(*tmp)[i] += from[nei];
        }
    }

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(out, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            to[i] += VECTOR(*tmp)[nei];
        }
    }

    return 0;
}

/* ARPACK auxiliary routine for the weighted HITS algorithm */
static int igraph_i_kleinberg_weighted(igraph_real_t *to,
                                       const igraph_real_t *from,
                                       int n, void *extra) {

    igraph_i_kleinberg_data2_t *data = (igraph_i_kleinberg_data2_t*)extra;
    igraph_inclist_t *in = data->in;
    igraph_inclist_t *out = data->out;
    igraph_vector_t *tmp = data->tmp;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *g = data->graph;
    igraph_vector_int_t *neis;
    long int i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(in, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei_edge = (long int) VECTOR(*neis)[j];
            long int nei = IGRAPH_OTHER(g, nei_edge, i);
            VECTOR(*tmp)[i] += from[nei] * VECTOR(*weights)[nei_edge];
        }
    }

    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(out, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei_edge = (long int) VECTOR(*neis)[j];
            long int nei = IGRAPH_OTHER(g, nei_edge, i);
            to[i] += VECTOR(*tmp)[nei] * VECTOR(*weights)[nei_edge];
        }
    }

    return 0;
}

static int igraph_i_kleinberg(const igraph_t *graph, igraph_vector_t *vector,
                              igraph_real_t *value, igraph_bool_t scale,
                              const igraph_vector_t *weights,
                              igraph_arpack_options_t *options, int inout) {

    igraph_adjlist_t myinadjlist, myoutadjlist;
    igraph_inclist_t myininclist, myoutinclist;
    igraph_adjlist_t *inadjlist, *outadjlist;
    igraph_inclist_t *ininclist, *outinclist;
    igraph_vector_t tmp;
    igraph_vector_t values;
    igraph_matrix_t vectors;
    igraph_i_kleinberg_data_t extra;
    igraph_i_kleinberg_data2_t extra2;
    long int i;

    if (igraph_ecount(graph) == 0 || igraph_vcount(graph) == 1) {
        /* special case: empty graph or single vertex */
        if (value) {
            *value = igraph_ecount(graph) ? 1.0 : IGRAPH_NAN;
        }
        if (vector) {
            igraph_vector_resize(vector, igraph_vcount(graph));
            igraph_vector_fill(vector, 1);
        }
        return IGRAPH_SUCCESS;
    }

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != igraph_ecount(graph)) {
            IGRAPH_ERROR("Invalid length of weights vector when calculating "
                         "hub or authority scores", IGRAPH_EINVAL);
        }
        /* Safe to call minmax, ecount == 0 case was caught earlier */
        IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));
        if (min == 0 && max == 0) {
            /* special case: all weights are zeros */
            if (value) {
                *value = IGRAPH_NAN;
            }
            if (vector) {
                igraph_vector_resize(vector, igraph_vcount(graph));
                igraph_vector_fill(vector, 1);
            }
            return IGRAPH_SUCCESS;
        }
    }

    options->n = igraph_vcount(graph);
    options->start = 1;   /* no random start vector */

    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);

    if (inout == 0) {
        inadjlist = &myinadjlist;
        outadjlist = &myoutadjlist;
        ininclist = &myininclist;
        outinclist = &myoutinclist;
    } else if (inout == 1) {
        inadjlist = &myoutadjlist;
        outadjlist = &myinadjlist;
        ininclist = &myoutinclist;
        outinclist = &myininclist;
    } else {
        /* This should not happen */
        IGRAPH_ERROR("Invalid 'inout' argument, please do not call "
                     "this function directly", IGRAPH_FAILURE);
    }

    if (weights == 0) {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &myinadjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &myinadjlist);
        IGRAPH_CHECK(igraph_adjlist_init(graph, &myoutadjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &myoutadjlist);
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &myininclist, IGRAPH_IN, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &myininclist);
        IGRAPH_CHECK(igraph_inclist_init(graph, &myoutinclist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &myoutinclist);
    }

    IGRAPH_CHECK(igraph_degree(graph, &tmp, igraph_vss_all(), IGRAPH_ALL, 0));
    for (i = 0; i < options->n; i++) {
        if (VECTOR(tmp)[i] != 0) {
            MATRIX(vectors, i, 0) = VECTOR(tmp)[i];
        } else {
            MATRIX(vectors, i, 0) = 1.0;
        }
    }

    extra.in = inadjlist; extra.out = outadjlist; extra.tmp = &tmp;
    extra2.in = ininclist; extra2.out = outinclist; extra2.tmp = &tmp;
    extra2.graph = graph; extra2.weights = weights;

    options->nev = 1;
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->which[0] = 'L'; options->which[1] = 'M';

    if (weights == 0) {
        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_kleinberg_unweighted, &extra,
                                           options, 0, &values, &vectors));
        igraph_adjlist_destroy(&myoutadjlist);
        igraph_adjlist_destroy(&myinadjlist);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_kleinberg_weighted, &extra2,
                                           options, 0, &values, &vectors));
        igraph_inclist_destroy(&myoutinclist);
        igraph_inclist_destroy(&myininclist);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    if (value) {
        *value = VECTOR(values)[0];
    }

    if (vector) {
        igraph_real_t amax = 0;
        long int which = 0;
        long int i;
        IGRAPH_CHECK(igraph_vector_resize(vector, options->n));
        for (i = 0; i < options->n; i++) {
            igraph_real_t tmp;
            VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
            tmp = fabs(VECTOR(*vector)[i]);
            if (tmp > amax) {
                amax = tmp;
                which = i;
            }
        }
        if (scale && amax != 0) {
            igraph_vector_scale(vector, 1 / VECTOR(*vector)[which]);
        } else if (igraph_i_vector_mostly_negative(vector)) {
            igraph_vector_scale(vector, -1.0);
        }

        /* Correction for numeric inaccuracies (eliminating -0.0) */
        for (i = 0; i < options->n; i++) {
            if (VECTOR(*vector)[i] < 0) {
                VECTOR(*vector)[i] = 0;
            }
        }
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
    }
    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_hub_score
 * \brief Kleinberg's hub scores.
 *
 * The hub scores of the vertices are defined as the principal
 * eigenvector of <code>A*A^T</code>, where <code>A</code> is the adjacency
 * matrix of the graph, <code>A^T</code> is its transposed.
 * </para><para>
 * See the following reference on the meaning of this score:
 * J. Kleinberg. Authoritative sources in a hyperlinked
 * environment. \emb Proc. 9th ACM-SIAM Symposium on Discrete
 * Algorithms, \eme 1998. Extended version in \emb Journal of the
 * ACM \eme 46(1999). Also appears as IBM Research Report RJ 10076, May
 * 1997.
 * \param graph The input graph. Can be directed and undirected.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. If a null pointer then it is ignored.
 * \param value If not a null pointer then the eigenvalue
 *    corresponding to the calculated eigenvector is stored here.
 * \param scale If not zero then the result will be scaled such that
 *     the absolute value of the maximum centrality is one.
 * \param weights A null pointer (=no edge weights), or a vector
 *     giving the weights of the edges.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 *
 * Time complexity: depends on the input graph, usually it is O(|V|),
 * the number of vertices.
 *
 * \sa \ref igraph_authority_score() for the companion measure,
 * \ref igraph_pagerank(), \ref igraph_personalized_pagerank(),
 * \ref igraph_eigenvector_centrality() for similar measures.
 */

int igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
                     igraph_real_t *value, igraph_bool_t scale,
                     const igraph_vector_t *weights,
                     igraph_arpack_options_t *options) {

    return igraph_i_kleinberg(graph, vector, value, scale, weights, options, 0);
}

/**
 * \function igraph_authority_score
 * \brief Kleinerg's authority scores.
 *
 * The authority scores of the vertices are defined as the principal
 * eigenvector of <code>A^T*A</code>, where <code>A</code> is the adjacency
 * matrix of the graph, <code>A^T</code> is its transposed.
 * </para><para>
 * See the following reference on the meaning of this score:
 * J. Kleinberg. Authoritative sources in a hyperlinked
 * environment. \emb Proc. 9th ACM-SIAM Symposium on Discrete
 * Algorithms, \eme 1998. Extended version in \emb Journal of the
 * ACM \eme 46(1999). Also appears as IBM Research Report RJ 10076, May
 * 1997.
 * \param graph The input graph. Can be directed and undirected.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. If a null pointer then it is ignored.
 * \param value If not a null pointer then the eigenvalue
 *    corresponding to the calculated eigenvector is stored here.
 * \param scale If not zero then the result will be scaled such that
 *     the absolute value of the maximum centrality is one.
 * \param weights A null pointer (=no edge weights), or a vector
 *     giving the weights of the edges.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 *
 * Time complexity: depends on the input graph, usually it is O(|V|),
 * the number of vertices.
 *
 * \sa \ref igraph_hub_score() for the companion measure,
 * \ref igraph_pagerank(), \ref igraph_personalized_pagerank(),
 * \ref igraph_eigenvector_centrality() for similar measures.
 */

int igraph_authority_score(const igraph_t *graph, igraph_vector_t *vector,
                           igraph_real_t *value, igraph_bool_t scale,
                           const igraph_vector_t *weights,
                           igraph_arpack_options_t *options) {

    return igraph_i_kleinberg(graph, vector, value, scale, weights, options, 1);
}

typedef struct igraph_i_pagerank_data_t {
    const igraph_t *graph;
    igraph_adjlist_t *adjlist;
    igraph_real_t damping;
    igraph_vector_t *outdegree;
    igraph_vector_t *tmp;
    igraph_vector_t *reset;
} igraph_i_pagerank_data_t;

typedef struct igraph_i_pagerank_data2_t {
    const igraph_t *graph;
    igraph_inclist_t *inclist;
    const igraph_vector_t *weights;
    igraph_real_t damping;
    igraph_vector_t *outdegree;
    igraph_vector_t *tmp;
    igraph_vector_t *reset;
} igraph_i_pagerank_data2_t;

static int igraph_i_pagerank(igraph_real_t *to, const igraph_real_t *from,
                             int n, void *extra) {

    igraph_i_pagerank_data_t *data = extra;
    igraph_adjlist_t *adjlist = data->adjlist;
    igraph_vector_t *outdegree = data->outdegree;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_t *reset = data->reset;
    igraph_vector_int_t *neis;
    long int i, j, nlen;
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
            long int nei = (long int) VECTOR(*neis)[j];
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

    return 0;
}

static int igraph_i_pagerank2(igraph_real_t *to, const igraph_real_t *from,
                              int n, void *extra) {

    igraph_i_pagerank_data2_t *data = extra;
    const igraph_t *graph = data->graph;
    igraph_inclist_t *inclist = data->inclist;
    const igraph_vector_t *weights = data->weights;
    igraph_vector_t *outdegree = data->outdegree;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_t *reset = data->reset;
    long int i, j, nlen;
    igraph_real_t sumfrom = 0.0;
    igraph_vector_int_t *neis;
    igraph_real_t fact = 1 - data->damping;

    /*
    printf("PageRank weighted: multiplying vector: ");
    for (i=0; i<n; i++) { printf(" %.4f", from[i]); }
    printf("\n");
    */

    for (i = 0; i < n; i++) {
        sumfrom += VECTOR(*outdegree)[i] != 0 ? from[i] * fact : from[i];
        VECTOR(*tmp)[i] = from[i] / VECTOR(*outdegree)[i];
    }

    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(inclist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = (long int) VECTOR(*neis)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
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

    return 0;
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
 * \c IGRAPH_PAGERANK_ALGO_ARPACK, based on the ARPACK library. This
 * was the default before igraph version 0.7. The second and recommended
 * implementation is \c IGRAPH_PAGERANK_ALGO_PRPACK. This is using the
 * PRPACK package, see https://github.com/dgleich/prpack .
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
 *
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_ARPACK, \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable, the eigenvalue
 *    corresponding to the PageRank vector is stored here. It should
 *    be always exactly one.
 * \param vids The vertex ids for which the PageRank is returned.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper).
 *    Must be a probability in the range [0, 1]. A commonly used value is 0.85.
 * \param weights Optional edge weights. May be a \c NULL pointer,
 *    meaning unweighted edges, or a vector of non-negative values
 *    of the same length as the number of edges.
 * \param options Options for the ARPACK method. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the <code>n</code> (number
 *    of vertices), <code>nev</code> (1), <code>ncv</code> (3) and <code>which</code>
 *    (LM) parameters and it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for temporary data.
 *         \c IGRAPH_EINVVID, invalid vertex id in \p vids.
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

int igraph_pagerank(const igraph_t *graph, igraph_pagerank_algo_t algo,
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
 * </para>
 *
 * <para>
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_ARPACK, \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable, the eigenvalue
 *    corresponding to the PageRank vector is stored here. It should
 *    be always exactly one.
 * \param vids The vertex ids for which the PageRank is returned.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper).
 *    Must be a probability in the range [0, 1]. A commonly used value is 0.85.
 * \param reset_vids IDs of the vertices used when resetting the random walk.
 * \param weights Optional edge weights, it is either a null pointer,
 *    then the edges are not weighted, or a vector of the same length
 *    as the number of edges.
 * \param options Options for the ARPACK method. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the <code>n</code> (number
 *    of vertices), <code>nev</code> (1), <code>ncv</code> (3) and <code>which</code>
 *    (LM) parameters and it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids or an empty reset vertex sequence in
 *         \p vids_reset.
 *
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 *
 * \sa \ref igraph_pagerank() for the non-personalized implementation.
 */

int igraph_personalized_pagerank_vs(const igraph_t *graph,
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
        VECTOR(reset)[(long int)IGRAPH_VIT_GET(vit)]++;
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

    return 0;
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
 * </para>
 *
 * <para>
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_ARPACK, \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable, the eigenvalue
 *    corresponding to the PageRank vector is stored here. It should
 *    be always exactly one.
 * \param vids The vertex ids for which the PageRank is returned.
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
 *    for details. Note that the function overwrites the <code>n</code> (number
 *    of vertices), <code>nev</code> (1), <code>ncv</code> (3) and <code>which</code>
 *    (LM) parameters and it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids or an invalid reset vector in \p reset.
 *
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 *
 * \sa \ref igraph_pagerank() for the non-personalized implementation,
 * \ref igraph_personalized_pagerank_vs() for a personalized implementation
 * with resetting to specific vertices.
 */
int igraph_personalized_pagerank(const igraph_t *graph,
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
                weights, options);
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
static int igraph_i_personalized_pagerank_arpack(const igraph_t *graph, igraph_vector_t *vector,
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

    long int i;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);

    if (reset && igraph_vector_size(reset) != no_of_nodes) {
        IGRAPH_ERROR("Invalid length of reset vector when calculating personalized PageRank scores.", IGRAPH_EINVAL);
    }

    if (no_of_edges == 0) {
        /* Special case: graph with no edges. Result is the same as the personalization vector. */
        if (value) {
            *value = 1.0;
        }
        if (vector) {
            IGRAPH_CHECK(igraph_vector_resize(vector, no_of_nodes));
            if (reset && no_of_nodes > 0) {
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

    options->n = (int) no_of_nodes;
    options->nev = 1;
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rnsolve */
    options->which[0] = 'L'; options->which[1] = 'M';
    options->start = 1;       /* no random start vector */

    directed = directed && igraph_is_directed(graph);

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Invalid length of weights vector when calculating PageRank scores.", IGRAPH_EINVAL);
        }

        /* Safe to call minmax, ecount == 0 case was caught earlier */
        IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));
        if (igraph_is_nan(min)) {
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
        double reset_sum, reset_min;
        reset_min = igraph_vector_min(reset);
        if (reset_min < 0) {
            IGRAPH_ERROR("The reset vector must not contain negative elements.", IGRAPH_EINVAL);
        }
        if (igraph_is_nan(reset_min)) {
            IGRAPH_ERROR("The reset vector must not contain NaN values.", IGRAPH_EINVAL);
        }
        reset_sum = igraph_vector_sum(reset);
        if (reset_sum == 0) {
            IGRAPH_ERROR("The sum of the elements in the reset vector must not be zero.", IGRAPH_EINVAL);
        }

        IGRAPH_CHECK(igraph_vector_copy(&normalized_reset, reset));
        IGRAPH_FINALLY(igraph_vector_destroy, &normalized_reset);

        igraph_vector_scale(&normalized_reset, 1.0 / reset_sum);
    }

    if (!weights) {

        igraph_adjlist_t adjlist;
        igraph_i_pagerank_data_t data;

        data.graph = graph;
        data.adjlist = &adjlist;
        data.damping = damping;
        data.outdegree = &outdegree;
        data.tmp = &tmp;
        data.reset = reset ? &normalized_reset : NULL;

        IGRAPH_CHECK(igraph_degree(graph, &outdegree, igraph_vss_all(),
                                   directed ? IGRAPH_OUT : IGRAPH_ALL, IGRAPH_LOOPS));
        IGRAPH_CHECK(igraph_degree(graph, &indegree, igraph_vss_all(),
                                   directed ? IGRAPH_IN : IGRAPH_ALL, IGRAPH_LOOPS));
        /* Set up an appropriate starting vector. We start from the in-degrees
         * plus some small random noise to avoid convergence problems */
        for (i = 0; i < options->n; i++) {
            if (VECTOR(indegree)[i]) {
                MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(-1e-4, 1e-4);
            } else {
                MATRIX(vectors, i, 0) = 1;
            }
        }

        IGRAPH_CHECK(igraph_adjlist_init(
            graph, &adjlist, dirmode, IGRAPH_LOOPS, IGRAPH_MULTIPLE
        ));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_pagerank,
                                           &data, options, 0, &values, &vectors));

        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);

    } else {

        igraph_inclist_t inclist;
        igraph_bool_t negative_weight_warned = 0;
        igraph_i_pagerank_data2_t data;

        data.graph = graph;
        data.inclist = &inclist;
        data.weights = weights;
        data.damping = damping;
        data.outdegree = &outdegree;
        data.tmp = &tmp;
        data.reset = reset ? &normalized_reset : NULL;

        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, dirmode, IGRAPH_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

        /* Weighted degree */
        for (i = 0; i < no_of_edges; i++) {
            long int from = IGRAPH_FROM(graph, i);
            long int to = IGRAPH_TO(graph, i);
            igraph_real_t weight = VECTOR(*weights)[i];
            if (weight < 0 && !negative_weight_warned) {
                IGRAPH_WARNING("Replacing negative weights with zeros during PageRank calculation.");
                weight = 0;
                negative_weight_warned = 1;
            }
            VECTOR(outdegree)[from] += weight;
            VECTOR(indegree) [to]   += weight;
            if (!directed) {
                VECTOR(outdegree)[to]   += weight;
                VECTOR(indegree) [from] += weight;
            }
        }
        /* Set up an appropriate starting vector. We start from the in-degrees
         * plus some small random noise to avoid convergence problems */
        for (i = 0; i < options->n; i++) {
            if (VECTOR(indegree)[i]) {
                MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(-1e-4, 1e-4);
            } else {
                MATRIX(vectors, i, 0) = 1;
            }
        }

        IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_pagerank2,
                                           &data, options, 0, &values, &vectors));

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
        long int i;
        igraph_vit_t vit;
        long int nodes_to_calc;
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
            VECTOR(*vector)[i] = MATRIX(vectors, (long int)IGRAPH_VIT_GET(vit), 0);
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
