/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_structural.h"
#include "igraph_interface.h"

#include "math/safe_intop.h"

#include <math.h>

static igraph_error_t igraph_i_laplacian_validate_weights(
    const igraph_t* graph,  const igraph_vector_t* weights
) {
    igraph_integer_t no_of_edges;

    if (weights == NULL) {
        return IGRAPH_SUCCESS;
    }

    no_of_edges = igraph_ecount(graph);

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (no_of_edges > 0) {
        igraph_real_t minweight = igraph_vector_min(weights);
        if (minweight < 0) {
            IGRAPH_ERROR("Weight vector must be non-negative.", IGRAPH_EINVAL);
        } else if (isnan(minweight)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_laplacian
 * \brief Returns the Laplacian matrix of a graph.
 *
 * The Laplacian matrix \c L of a graph is defined as
 * <code>L_ij = - A_ij</code> when <code>i != j</code> and
 * <code>L_ii = d_i - A_ii</code>. Here \c A denotes the (possibly weighted)
 * adjacency matrix and <code>d_i</code> is the degree (or strength, if weighted)
 * of vertex \c i. In directed graphs, the \p mode parameter controls whether to use
 * out- or in-degrees. Correspondingly, the rows or columns will sum to zero.
 * In undirected graphs, <code>A_ii</code> is taken to be \em twice the number
 * (or total weight) of self-loops, ensuring that <code>d_i = \sum_j A_ij</code>.
 * Thus, the Laplacian of an undirected graph is the same as the Laplacian
 * of a directed one obtained by replacing each undirected edge with two reciprocal
 * directed ones.
 *
 * </para><para>
 * More compactly, <code>L = D - A</code> where the \c D is a diagonal matrix
 * containing the degrees. The Laplacian matrix can also be normalized, with several
 * conventional normalization methods. See \ref igraph_laplacian_normalization_t for
 * the methods available in igraph.
 *
 * </para><para>
 * The first version of this function was written by Vincent Matossian.
 *
 * \param graph Pointer to the graph to convert.
 * \param res Pointer to an initialized matrix object, the result is
 *        stored here. It will be resized if needed.
 * \param mode Controls whether to use out- or in-degrees in directed graphs.
 *        If set to \c IGRAPH_ALL, edge directions will be ignored.
 * \param normalization The normalization method to use when calculating the
 *        Laplacian matrix. See \ref igraph_laplacian_normalization_t for
 *        possible values.
 * \param weights An optional vector containing non-negative edge weights,
 *        to calculate the weighted Laplacian matrix. Set it to a null pointer to
 *        calculate the unweighted Laplacian.
 * \return Error code.
 *
 * Time complexity: O(|V|^2), |V| is the number of vertices in the graph.
 *
 * \example examples/simple/igraph_get_laplacian.c
 */

igraph_error_t igraph_get_laplacian(
    const igraph_t *graph, igraph_matrix_t *res, igraph_neimode_t mode,
    igraph_laplacian_normalization_t normalization,
    const igraph_vector_t *weights
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_t degree;
    igraph_integer_t i;

    IGRAPH_ASSERT(res != NULL);

    IGRAPH_CHECK(igraph_i_laplacian_validate_weights(graph, weights));

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_null(res);

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(), mode, IGRAPH_LOOPS, weights));

    /* Value of 'mode' is validated in igraph_strength() call above. */
    if (! directed) {
        mode = IGRAPH_ALL;
    } else if (mode == IGRAPH_ALL) {
        directed = 0;
    }

    for (i = 0; i < no_of_nodes; i++) {
        switch (normalization) {
        case IGRAPH_LAPLACIAN_UNNORMALIZED:
            MATRIX(*res, i, i) = VECTOR(degree)[i];
            break;

        case IGRAPH_LAPLACIAN_SYMMETRIC:
            if (VECTOR(degree)[i] > 0) {
                MATRIX(*res, i, i) = 1;
                VECTOR(degree)[i] = 1.0 / sqrt(VECTOR(degree)[i]);
            }
            break;

        case IGRAPH_LAPLACIAN_LEFT:
        case IGRAPH_LAPLACIAN_RIGHT:
            if (VECTOR(degree)[i] > 0) {
                MATRIX(*res, i, i) = 1;
                VECTOR(degree)[i] = 1.0 / VECTOR(degree)[i];
            }
            break;

        default:
            IGRAPH_ERROR("Invalid Laplacian normalization method.", IGRAPH_EINVAL);
        }
    }

    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = IGRAPH_FROM(graph, i);
        igraph_integer_t to   = IGRAPH_TO(graph, i);
        igraph_real_t weight  = weights ? VECTOR(*weights)[i] : 1.0;
        igraph_real_t norm;

        switch (normalization) {
        case IGRAPH_LAPLACIAN_UNNORMALIZED:
            MATRIX(*res, from, to) -= weight;
            if (!directed) {
                MATRIX(*res, to, from) -= weight;
            }
            break;

        case IGRAPH_LAPLACIAN_SYMMETRIC:
            norm = VECTOR(degree)[from] * VECTOR(degree)[to];
            if (norm == 0 && weight != 0) {
                IGRAPH_ERRORF(
                    "Found non-isolated vertex with zero %s-%s, "
                    "cannot perform symmetric normalization of Laplacian with '%s' mode.",
                    IGRAPH_EINVAL,
                    mode == IGRAPH_OUT ? "out" : "in", weights ? "strength" : "degree", mode == IGRAPH_OUT ? "out" : "in");
            }
            weight *= norm;
            MATRIX(*res, from, to) -= weight;
            if (!directed) {
                MATRIX(*res, to, from) -= weight;
            }
            break;

        case IGRAPH_LAPLACIAN_LEFT:
            norm = VECTOR(degree)[from];
            if (norm == 0 && weight != 0) {
                IGRAPH_ERRORF(
                    "Found non-isolated vertex with zero in-%s, "
                    "cannot perform left stochastic normalization of Laplacian with 'in' mode.",
                    IGRAPH_EINVAL,
                    weights ? "strength" : "degree");
            }
            MATRIX(*res, from, to) -= weight * norm;
            if (!directed) {
                /* no failure possible in undirected case, as zero degrees occur only for isolated vertices */
                MATRIX(*res, to, from) -= weight * VECTOR(degree)[to];
            }
            break;

        case IGRAPH_LAPLACIAN_RIGHT:
            norm = VECTOR(degree)[to];
            if (norm == 0 && weight != 0) {
                IGRAPH_ERRORF(
                    "Found non-isolated vertex with zero out-%s, "
                    "cannot perform right stochastic normalization of Laplacian with 'out' mode.",
                    IGRAPH_EINVAL,
                    weights ? "strength" : "degree");
            }
            MATRIX(*res, from, to) -= weight * norm;
            if (!directed) {
                /* no failure possible in undirected case, as zero degrees occur only for isolated vertices */
                MATRIX(*res, to, from) -= weight * VECTOR(degree)[from];
            }
            break;
        }
    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_get_laplacian_sparse
 * \brief Returns the Laplacian of a graph in a sparse matrix format.
 *
 * See \ref igraph_get_laplacian() for the definition of the Laplacian matrix.
 *
 * </para><para>
 * The first version of this function was written by Vincent Matossian.
 *
 * \param graph Pointer to the graph to convert.
 * \param sparseres Pointer to an initialized sparse matrix object, the
 *        result is stored here.
 * \param mode Controls whether to use out- or in-degrees in directed graphs.
 *        If set to \c IGRAPH_ALL, edge directions will be ignored.
 * \param normalization The normalization method to use when calculating the
 *        Laplacian matrix. See \ref igraph_laplacian_normalization_t for
 *        possible values.
 * \param weights An optional vector containing non-negative edge weights,
 *        to calculate the weighted Laplacian matrix. Set it to a null pointer to
 *        calculate the unweighted Laplacian.
 * \return Error code.
 *
 * Time complexity: O(|E|), |E| is the number of edges in the graph.
 *
 * \example examples/simple/igraph_get_laplacian_sparse.c
 */

igraph_error_t igraph_get_laplacian_sparse(
    const igraph_t *graph, igraph_sparsemat_t *sparseres, igraph_neimode_t mode,
    igraph_laplacian_normalization_t normalization,
    const igraph_vector_t *weights
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_t degree;
    igraph_integer_t i;
    igraph_integer_t nz;

    if (directed) {
        IGRAPH_SAFE_ADD(no_of_edges, no_of_nodes, &nz);
    } else {
        IGRAPH_SAFE_ADD(no_of_edges * 2, no_of_nodes, &nz);
    }

    IGRAPH_ASSERT(sparseres != NULL);

    IGRAPH_CHECK(igraph_i_laplacian_validate_weights(graph, weights));

    IGRAPH_CHECK(igraph_sparsemat_resize(sparseres, no_of_nodes, no_of_nodes, nz));

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(), mode, IGRAPH_LOOPS, weights));

    for (i = 0; i < no_of_nodes; i++) {
        switch (normalization) {
        case IGRAPH_LAPLACIAN_UNNORMALIZED:
            IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, VECTOR(degree)[i]));
            break;

        case IGRAPH_LAPLACIAN_SYMMETRIC:
            if (VECTOR(degree)[i] > 0) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, 1));
                VECTOR(degree)[i] = 1.0 / sqrt(VECTOR(degree)[i]);
            }
            break;

        case IGRAPH_LAPLACIAN_LEFT:
        case IGRAPH_LAPLACIAN_RIGHT:
            if (VECTOR(degree)[i] > 0) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, 1));
                VECTOR(degree)[i] = 1.0 / VECTOR(degree)[i];
            }
            break;

        default:
            IGRAPH_ERROR("Invalid Laplacian normalization method.", IGRAPH_EINVAL);
        }
    }

    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = IGRAPH_FROM(graph, i);
        igraph_integer_t to   = IGRAPH_TO(graph, i);
        igraph_real_t weight  = weights ? VECTOR(*weights)[i] : 1.0;
        igraph_real_t norm;

        switch (normalization) {
        case IGRAPH_LAPLACIAN_UNNORMALIZED:
            IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -weight));
            if (!directed) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -weight));
            }
            break;

        case IGRAPH_LAPLACIAN_SYMMETRIC:
            norm = VECTOR(degree)[from] * VECTOR(degree)[to];
            if (norm == 0 && weight != 0) {
                IGRAPH_ERRORF(
                    "Found non-isolated vertex with zero %s-%s, "
                    "cannot perform symmetric normalization of Laplacian with '%s' mode.",
                    IGRAPH_EINVAL,
                    mode == IGRAPH_OUT ? "out" : "in", weights ? "strength" : "degree", mode == IGRAPH_OUT ? "out" : "in");
            }
            weight *= norm;
            IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -weight));
            if (!directed) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -weight));
            }
            break;

        case IGRAPH_LAPLACIAN_LEFT:
            norm = VECTOR(degree)[from];
            if (norm == 0 && weight != 0) {
                IGRAPH_ERRORF(
                    "Found non-isolated vertex with zero in-%s, "
                    "cannot perform left stochastic normalization of Laplacian with 'in' mode.",
                    IGRAPH_EINVAL,
                    weights ? "strength" : "degree");
            }
            IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -weight * norm));
            if (!directed) {
                /* no failure possible in undirected case, as zero degrees occur only for isolated vertices */
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -weight * VECTOR(degree)[to]));
            }
            break;

        case IGRAPH_LAPLACIAN_RIGHT:
            norm = VECTOR(degree)[to];
            if (norm == 0 && weight != 0) {
                IGRAPH_ERRORF(
                    "Found non-isolated vertex with zero out-%s, "
                    "cannot perform right stochastic normalization of Laplacian with 'out' mode.",
                    IGRAPH_EINVAL,
                    weights ? "strength" : "degree");
            }
            IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -weight * norm));
            if (!directed) {
                /* no failure possible in undirected case, as zero degrees occur only for isolated vertices */
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -weight * VECTOR(degree)[from]));
            }
            break;
        }
    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_laplacian
 * \brief Returns the Laplacian matrix of a graph (deprecated).
 *
 * This function produces the Laplacian matrix of a graph in either dense or
 * sparse format. When \p normalized is set to true, the type of normalization
 * used depends on the directnedness of the graph: symmetric normalization
 * is used for undirected graphs and left stochastic normalization for
 * directed graphs.
 *
 * \param graph Pointer to the graph to convert.
 * \param res Pointer to an initialized matrix object or \c NULL. The dense matrix
 *        result will be stored here.
 * \param sparseres Pointer to an initialized sparse matrix object or \c NULL.
 *        The sparse matrix result will be stored here.
 * \param mode Controls whether to use out- or in-degrees in directed graphs.
 *        If set to \c IGRAPH_ALL, edge directions will be ignored.
 * \param normalized Boolean, whether to normalize the result.
 * \param weights An optional vector containing non-negative edge weights,
 *        to calculate the weighted Laplacian matrix. Set it to a null pointer to
 *        calculate the unweighted Laplacian.
 * \return Error code.
 *
 * \deprecated-by igraph_get_laplacian 0.10.0
 */

igraph_error_t igraph_laplacian(
    const igraph_t *graph, igraph_matrix_t *res, igraph_sparsemat_t *sparseres,
    igraph_bool_t normalized, const igraph_vector_t *weights
) {
    igraph_laplacian_normalization_t norm_method = IGRAPH_LAPLACIAN_UNNORMALIZED;

    if (!res && !sparseres) {
        IGRAPH_ERROR("Laplacian: specify at least one of 'res' or 'sparseres'",
                     IGRAPH_EINVAL);
    }

    if (normalized) {
        if (igraph_is_directed(graph)) {
            norm_method = IGRAPH_LAPLACIAN_LEFT;
        } else {
            norm_method = IGRAPH_LAPLACIAN_SYMMETRIC;
        }
    }

    if (res) {
        IGRAPH_CHECK(igraph_get_laplacian(graph, res, IGRAPH_OUT, norm_method, weights));
    }

    if (sparseres) {
        IGRAPH_CHECK(igraph_get_laplacian_sparse(graph, sparseres, IGRAPH_OUT, norm_method, weights));
    }

    return IGRAPH_SUCCESS;
}
