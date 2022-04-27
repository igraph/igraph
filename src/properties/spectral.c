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

#include <math.h>

static igraph_error_t igraph_i_get_laplacian_unweighted(
    const igraph_t *graph, igraph_matrix_t *res, igraph_bool_t normalized
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t from, to;
    igraph_vector_t degree;
    igraph_integer_t i;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_null(res);

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS, 0));

    if (directed) {

        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                MATRIX(*res, i, i) = VECTOR(degree)[i];
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    MATRIX(*res, from, to) -= 1;
                }
            }
        } else {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(degree)[i] > 0) {
                    MATRIX(*res, i, i) = 1;
                    VECTOR(degree)[i] = 1.0 / VECTOR(degree)[i];
                }
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    MATRIX(*res, from, to) -= VECTOR(degree)[from];
                }
            }
        }

    } else {

        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                MATRIX(*res, i, i) = VECTOR(degree)[i];
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    MATRIX(*res, to, from) -= 1;
                    MATRIX(*res, from, to) -= 1;
                }
            }
        } else {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(degree)[i] > 0) {
                    MATRIX(*res, i, i) = 1;
                    VECTOR(degree)[i] = sqrt(VECTOR(degree)[i]);
                }
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    double diff = 1.0 / (VECTOR(degree)[from] * VECTOR(degree)[to]);
                    MATRIX(*res, from, to) -= diff;
                    MATRIX(*res, to, from) -= diff;
                }
            }
        }

    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_get_laplacian_unweighted_sparse(
    const igraph_t *graph, igraph_sparsemat_t *sparseres, igraph_bool_t normalized
) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t from, to;
    igraph_vector_t degree;
    igraph_integer_t i;

    igraph_integer_t nz = directed ? no_of_edges + no_of_nodes :
                no_of_edges * 2 + no_of_nodes;
    IGRAPH_CHECK(igraph_sparsemat_resize(sparseres, no_of_nodes, no_of_nodes, nz));

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS, 0));

    if (directed) {
        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, VECTOR(degree)[i]));
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -1.0));
                }
            }
        } else {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(degree)[i] > 0) {
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, 1));
                    VECTOR(degree)[i] = 1.0 / VECTOR(degree)[i];
                }
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    IGRAPH_CHECK(
                        igraph_sparsemat_entry(sparseres, from, to, -VECTOR(degree)[from])
                    );
                }
            }
        }

    } else {

        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i,
                                                    VECTOR(degree)[i]));
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -1.0));
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -1.0));
                }
            }
        } else {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(degree)[i] > 0) {
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, 1));
                    VECTOR(degree)[i] = sqrt(VECTOR(degree)[i]);
                }
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, i, &from, &to);
                if (from != to) {
                    double diff = 1.0 / (VECTOR(degree)[from] * VECTOR(degree)[to]);
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -diff));
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -diff));
                }
            }
        }

    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_get_laplacian_weighted(
    const igraph_t *graph, igraph_matrix_t *res,
    igraph_bool_t normalized, const igraph_vector_t *weights
) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_t degree;
    igraph_integer_t i;

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid edge weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_null(res);

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS, weights));

    if (directed) {

        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                MATRIX(*res, i, i) = VECTOR(degree)[i];
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    MATRIX(*res, from, to) -= weight;
                }
            }
        } else { /* normalized */
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(degree)[i] > 0) {
                    MATRIX(*res, i, i) = 1;
                }
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    MATRIX(*res, from, to) -= weight / VECTOR(degree)[from];
                }
            }
        }

    } else { /* undirected */

        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                MATRIX(*res, i, i) = VECTOR(degree)[i];
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    MATRIX(*res, from, to) -= weight;
                    MATRIX(*res, to, from) -= weight;
                }
            }
        } else { /* normalized */
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(degree)[i] > 0) {
                    MATRIX(*res, i, i) = 1;
                    VECTOR(degree)[i] = sqrt(VECTOR(degree)[i]);
                }
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    igraph_real_t diff = weight / (VECTOR(degree)[from] * VECTOR(degree)[to]);
                    MATRIX(*res, from, to) -= diff;
                    MATRIX(*res, to, from) -= diff;
                }
            }
        }

    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_get_laplacian_weighted_sparse(
    const igraph_t *graph, igraph_sparsemat_t *sparseres,
    igraph_bool_t normalized, const igraph_vector_t *weights
) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_t degree;
    igraph_integer_t i;
    igraph_integer_t nz = directed ? no_of_edges + no_of_nodes : no_of_edges * 2 + no_of_nodes;

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid edge weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_sparsemat_resize(sparseres, no_of_nodes, no_of_nodes, nz));

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS, weights));

    if (directed) {

        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, VECTOR(degree)[i]));
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -weight));
                }
            }
        } else { /* normalized */
            for (i = 0; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, VECTOR(degree)[i] > 0 ? 1 : 0));
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    IGRAPH_CHECK(
                        igraph_sparsemat_entry(sparseres, from, to, -weight / VECTOR(degree)[from])
                    );
                }
            }
        }

    } else { /* undirected */

        if (!normalized) {
            for (i = 0; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, VECTOR(degree)[i]));
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -weight));
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -weight));
                }
            }
        } else { /* normalized */
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(degree)[i] > 0) {
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, i, i, 1));
                    VECTOR(degree)[i] = sqrt(VECTOR(degree)[i]);
                }
            }
            for (i = 0; i < no_of_edges; i++) {
                igraph_integer_t from = IGRAPH_FROM(graph, i);
                igraph_integer_t to   = IGRAPH_TO(graph, i);
                igraph_real_t weight  = VECTOR(*weights)[i];
                if (from != to) {
                    igraph_real_t diff = weight / (VECTOR(degree)[from] * VECTOR(degree)[to]);
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, from, to, -diff));
                    IGRAPH_CHECK(igraph_sparsemat_entry(sparseres, to, from, -diff));
                }
            }
        }

    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_laplacian
 * \brief Returns the Laplacian matrix of a graph.
 *
 * </para><para>
 * The graph Laplacian matrix is similar to an adjacency matrix but
 * contains -1's instead of 1's and the vertex degrees are included in
 * the diagonal. So the result for edge i--j is -1 if i!=j and is equal
 * to the degree of vertex i if i==j. This function will work on a
 * directed graph; in this case, the diagonal will contain the out-degrees.
 * Loop edges will be ignored.
 *
 * </para><para>
 * The normalized version of the Laplacian matrix has 1 in the diagonal and
 * -1/sqrt(d[i]d[j]) if there is an edge from i to j.
 *
 * </para><para>
 * The first version of this function was written by Vincent Matossian.
 *
 * \param graph Pointer to the graph to convert.
 * \param res Pointer to an initialized matrix object, the result is
 *        stored here. It will be resized if needed.
 * \param normalized Whether to create a normalized Laplacian matrix.
 * \param weights An optional vector containing edge weights, to calculate
 *        the weighted Laplacian matrix. Set it to a null pointer to
 *        calculate the unweighted Laplacian.
 * \return Error code.
 *
 * Time complexity: O(|V||V|), |V| is the number of vertices in the graph.
 *
 * \example examples/simple/igraph_get_laplacian.c
 */

igraph_error_t igraph_get_laplacian(
    const igraph_t *graph, igraph_matrix_t *res, igraph_bool_t normalized,
    const igraph_vector_t *weights
) {
    IGRAPH_ASSERT(res != NULL);
    if (weights) {
        return igraph_i_get_laplacian_weighted(graph, res, normalized, weights);
    } else {
        return igraph_i_get_laplacian_unweighted(graph, res, normalized);
    }
}


/**
 * \function igraph_get_laplacian_sparse
 * \brief Returns the Laplacian matrix of a graph in a sparse matrix format.
 *
 * </para><para>
 * See \ref igraph_get_laplacian() for the definition of the Laplacian matrix.
 *
 * </para><para>
 * The first version of this function was written by Vincent Matossian.
 * \param graph Pointer to the graph to convert.
 * \param sparseres Pointer to an initialized sparse matrix object, the
 *        result is stored here.
 * \param normalized Whether to create a normalized Laplacian matrix.
 * \param weights An optional vector containing edge weights, to calculate
 *        the weighted Laplacian matrix. Set it to a null pointer to
 *        calculate the unweighted Laplacian.
 * \return Error code.
 *
 * Time complexity: O(|E|), |E| is the number of edges in the graph.
 *
 * \example examples/simple/igraph_get_laplacian_sparse.c
 */

igraph_error_t igraph_get_laplacian_sparse(
    const igraph_t *graph, igraph_sparsemat_t *sparseres, igraph_bool_t normalized,
    const igraph_vector_t *weights
) {
    IGRAPH_ASSERT(sparseres != NULL);
    if (weights) {
        return igraph_i_get_laplacian_weighted_sparse(graph, sparseres, normalized, weights);
    } else {
        return igraph_i_get_laplacian_unweighted_sparse(graph, sparseres, normalized);
    }
}

/**
 * \function igraph_laplacian
 * \brief Returns the Laplacian matrix of a graph
 *
 * </para><para>
 * The graph Laplacian matrix is similar to an adjacency matrix but
 * contains -1's instead of 1's and the vertex degrees are included in
 * the diagonal. So the result for edge i--j is -1 if i!=j and is equal
 * to the degree of vertex i if i==j. igraph_laplacian will work on a
 * directed graph; in this case, the diagonal will contain the out-degrees.
 * Loop edges will be ignored.
 *
 * </para><para>
 * The normalized version of the Laplacian matrix has 1 in the diagonal and
 * -1/sqrt(d[i]d[j]) if there is an edge from i to j.
 *
 * </para><para>
 * The first version of this function was written by Vincent Matossian.
 * \param graph Pointer to the graph to convert.
 * \param res Pointer to an initialized matrix object, the result is
 *        stored here. It will be resized if needed.
 *        If it is a null pointer, then it is ignored.
 *        At least one of \p res and \p sparseres must be a non-null pointer.
 * \param sparseres Pointer to an initialized sparse matrix object, the
 *        result is stored here, if it is not a null pointer.
 *        At least one of \p res and \p sparseres must be a non-null pointer.
 * \param normalized Whether to create a normalized Laplacian matrix.
 * \param weights An optional vector containing edge weights, to calculate
 *        the weighted Laplacian matrix. Set it to a null pointer to
 *        calculate the unweighted Laplacian.
 * \return Error code.
 *
 * Time complexity: O(|V||V|),
 * |V| is the
 * number of vertices in the graph.
 *
 * \example examples/simple/igraph_laplacian.c
 */

igraph_error_t igraph_laplacian(const igraph_t *graph, igraph_matrix_t *res,
                     igraph_sparsemat_t *sparseres,
                     igraph_bool_t normalized,
                     const igraph_vector_t *weights) {
    if (!res && !sparseres) {
        IGRAPH_ERROR("Laplacian: specify at least one of `res' or `sparseres'",
                     IGRAPH_EINVAL);
    }

    if (res) {
        IGRAPH_CHECK(igraph_get_laplacian(graph, res, normalized, weights));
    }

    if (sparseres) {
        IGRAPH_CHECK(igraph_get_laplacian_sparse(graph, sparseres, normalized, weights));
    }

    return IGRAPH_SUCCESS;
}
