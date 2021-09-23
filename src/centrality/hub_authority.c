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
#include "igraph_structural.h"
#include "igraph_blas.h"

#include "centrality/centrality_internal.h"

#include "config.h"

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

/* ARPACK auxiliary routine for the unweighted HITS algorithm */
static igraph_error_t igraph_i_kleinberg_unweighted(igraph_real_t *to,
                                         const igraph_real_t *from,
                                         int n, void *extra) {
    igraph_i_kleinberg_data_t *data = (igraph_i_kleinberg_data_t*)extra;
    igraph_adjlist_t *in = data->in;
    igraph_adjlist_t *out = data->out;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    igraph_integer_t i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(in, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_integer_t nei = VECTOR(*neis)[j];
            VECTOR(*tmp)[i] += from[nei];
        }
    }

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(out, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_integer_t nei = VECTOR(*neis)[j];
            to[i] += VECTOR(*tmp)[nei];
        }
    }

    return IGRAPH_SUCCESS;
}

/* ARPACK auxiliary routine for the weighted HITS algorithm */
static igraph_error_t igraph_i_kleinberg_weighted(igraph_real_t *to,
                                       const igraph_real_t *from,
                                       int n, void *extra) {

    igraph_i_kleinberg_data2_t *data = (igraph_i_kleinberg_data2_t*)extra;
    igraph_inclist_t *in = data->in;
    igraph_inclist_t *out = data->out;
    igraph_vector_t *tmp = data->tmp;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *g = data->graph;
    igraph_vector_int_t *neis;
    igraph_integer_t i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(in, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_integer_t nei_edge = VECTOR(*neis)[j];
            igraph_integer_t nei = IGRAPH_OTHER(g, nei_edge, i);
            VECTOR(*tmp)[i] += from[nei] * VECTOR(*weights)[nei_edge];
        }
    }

    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(out, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            igraph_integer_t nei_edge = VECTOR(*neis)[j];
            igraph_integer_t nei = IGRAPH_OTHER(g, nei_edge, i);
            to[i] += VECTOR(*tmp)[nei] * VECTOR(*weights)[nei_edge];
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hub_and_authority_scores
 * \brief Kleinerg's hub and authority scores.
 *
 * The hub and authority scores of the vertices are defined as the principal
 * eigenvectors of <code>A*A^T</code> and <code>A^T*A</code>, respectively, 
 * where <code>A</code> is the adjacency
 * matrix of the graph, <code>A^T</code> is its transposed.
 * </para><para>
 * See the following reference on the meaning of this score:
 * J. Kleinberg. Authoritative sources in a hyperlinked
 * environment. \emb Proc. 9th ACM-SIAM Symposium on Discrete
 * Algorithms, \eme 1998. Extended version in \emb Journal of the
 * ACM \eme 46(1999). Also appears as IBM Research Report RJ 10076, May
 * 1997.
 * \param graph The input graph. Can be directed and undirected.
 * \param hub_vector Pointer to an initialized vector, the hub scores are 
 *    stored here. If a null pointer then it is ignored.
 * \param authority_vector Pointer to an initialized vector, the authority scores are 
 *    stored here. If a null pointer then it is ignored.
 * \param value If not a null pointer then the eigenvalue
 *    corresponding to the calculated eigenvectors is stored here.
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
 * \sa \ref igraph_hub_score(), \ref igraph_authority_score()
 * for the separate calculations, 
 * \ref igraph_pagerank(), \ref igraph_personalized_pagerank(),
 * \ref igraph_eigenvector_centrality() for similar measures.
 */
igraph_error_t igraph_hub_and_authority_scores(const igraph_t *graph,
        igraph_vector_t *hub_vector, igraph_vector_t *authority_vector,
        igraph_real_t *value, igraph_bool_t scale,
        const igraph_vector_t *weights, igraph_arpack_options_t *options) {

    igraph_adjlist_t myinadjlist, myoutadjlist;
    igraph_inclist_t myininclist, myoutinclist;
    igraph_adjlist_t *inadjlist, *outadjlist;
    igraph_inclist_t *ininclist, *outinclist;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_t tmp;
    igraph_vector_t values;
    igraph_matrix_t vectors;
    igraph_i_kleinberg_data_t extra;
    igraph_i_kleinberg_data2_t extra2;
    igraph_integer_t i;
    igraph_vector_t *my_hub_vector_p;
    igraph_vector_t my_hub_vector;


    if (igraph_ecount(graph) == 0 || no_of_nodes == 1) {
        /* special case: empty graph or single vertex */
        if (value) {
            *value = igraph_ecount(graph) ? 1.0 : IGRAPH_NAN;
        }
        if (hub_vector) {
            igraph_vector_resize(hub_vector, no_of_nodes);
            igraph_vector_fill(hub_vector, 1.0);
        }
        if (authority_vector) {
            igraph_vector_resize(authority_vector, no_of_nodes);
            igraph_vector_fill(authority_vector, 1.0);
        }
        return IGRAPH_SUCCESS;
    }

    if (weights) {
        igraph_real_t min, max;

        if (igraph_vector_size(weights) != igraph_ecount(graph)) {
            IGRAPH_ERRORF(
                    "Weights vector length (%d) should match number of "
                    "edges (%d) when calculating "
                    "hub or authority scores.",
                    IGRAPH_EINVAL,
                    igraph_vector_size(weights),
                    igraph_ecount(graph));
        }
        /* Safe to call minmax, ecount == 0 case was caught earlier */
        IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));
        if (min == 0 && max == 0) {
            /* special case: all weights are zeros */
            if (value) {
                *value = IGRAPH_NAN;
            }
            if (hub_vector) {
                igraph_vector_resize(hub_vector, no_of_nodes);
                igraph_vector_fill(hub_vector, 1);
            }
            if (authority_vector) {
                igraph_vector_resize(authority_vector, no_of_nodes);
                igraph_vector_fill(authority_vector, 1);
            }
            return IGRAPH_SUCCESS;
        }
    }

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph has too many vertices for ARPACK", IGRAPH_EOVERFLOW);
    }

    options->n = (int) no_of_nodes;
    options->start = 1;   /* no random start vector */

    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);

    inadjlist = &myinadjlist;
    outadjlist = &myoutadjlist;
    ininclist = &myininclist;
    outinclist = &myoutinclist;

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

    IGRAPH_CHECK(igraph_strength(graph, &tmp, igraph_vss_all(), IGRAPH_ALL, 0, 0));
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

    if (hub_vector || authority_vector) {
        if (!hub_vector) {
            IGRAPH_VECTOR_INIT_FINALLY(&my_hub_vector, options->n);
            my_hub_vector_p = &my_hub_vector;
        } else {
            my_hub_vector_p = hub_vector;
        }
        igraph_real_t amax = 0;
        igraph_integer_t which = 0;
        igraph_integer_t i;
        IGRAPH_CHECK(igraph_vector_resize(my_hub_vector_p, options->n));
        for (i = 0; i < options->n; i++) {
            igraph_real_t tmp;
            VECTOR(*my_hub_vector_p)[i] = MATRIX(vectors, i, 0);
            tmp = fabs(VECTOR(*my_hub_vector_p)[i]);
            if (tmp > amax) {
                amax = tmp;
                which = i;
            }
        }
        if (scale && amax != 0) {
            igraph_vector_scale(my_hub_vector_p, 1 / VECTOR(*my_hub_vector_p)[which]);
        } else if (igraph_i_vector_mostly_negative(my_hub_vector_p)) {
            igraph_vector_scale(my_hub_vector_p, -1.0);
        }

        /* Correction for numeric inaccuracies (eliminating -0.0) */
        for (i = 0; i < options->n; i++) {
            if (VECTOR(*my_hub_vector_p)[i] < 0) {
                VECTOR(*my_hub_vector_p)[i] = 0;
            }
        }
    }

    if (options->info) {
        IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
    }
    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    if (authority_vector) {
        igraph_real_t norm;
        igraph_vector_resize(authority_vector, igraph_vcount(graph));
        igraph_vector_null(authority_vector);
        IGRAPH_CHECK(igraph_inclist_init(graph, outinclist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, outinclist);
        for (int i = 0; i < igraph_vcount(graph); i++) {
            igraph_vector_int_t *inc = igraph_inclist_get(outinclist, i);
            for (int j = 0; j < igraph_vector_int_size(inc); j++) {
                igraph_integer_t edge = VECTOR(*inc)[j];
                if (weights) {
                    VECTOR(*authority_vector)[IGRAPH_TO(graph, edge)] += VECTOR(*weights)[edge] * VECTOR(*my_hub_vector_p)[IGRAPH_FROM(graph, edge)];
                } else {
                    VECTOR(*authority_vector)[IGRAPH_TO(graph, edge)] += VECTOR(*my_hub_vector_p)[IGRAPH_FROM(graph, edge)];
                }
            }
        }
        if (!scale) {
            norm = 1.0 / igraph_blas_dnrm2(authority_vector);
        } else {
            norm = 1.0 / igraph_vector_max(authority_vector);
        }
        igraph_vector_scale(authority_vector, norm);
        igraph_inclist_destroy(outinclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (!hub_vector && authority_vector) {
        igraph_vector_destroy(&my_hub_vector);
        IGRAPH_FINALLY_CLEAN(1);
    }
    return IGRAPH_SUCCESS;
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

igraph_error_t igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
                     igraph_real_t *value, igraph_bool_t scale,
                     const igraph_vector_t *weights,
                     igraph_arpack_options_t *options) {
    return igraph_hub_and_authority_scores(graph, vector, NULL, value, scale, weights, options);
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

igraph_error_t igraph_authority_score(const igraph_t *graph, igraph_vector_t *vector,
                           igraph_real_t *value, igraph_bool_t scale,
                           const igraph_vector_t *weights,
                           igraph_arpack_options_t *options) {
    return igraph_hub_and_authority_scores(graph, NULL, vector, value, scale, weights, options);
}
