/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

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

#include "igraph_layout.h"

#include "igraph_blas.h"
#include "igraph_components.h"
#include "igraph_eigen.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_paths.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include <limits.h>

static igraph_error_t igraph_i_layout_mds_step(igraph_real_t *to, const igraph_real_t *from,
                                    int n, void *extra);

static igraph_error_t igraph_i_layout_mds_single(const igraph_t* graph, igraph_matrix_t *res,
                                      igraph_matrix_t *dist, igraph_integer_t dim);

static igraph_error_t igraph_i_layout_mds_step(igraph_real_t *to, const igraph_real_t *from,
                                    int n, void *extra) {
    igraph_matrix_t* matrix = (igraph_matrix_t*)extra;
    IGRAPH_UNUSED(n);
    IGRAPH_CHECK(igraph_blas_dgemv_array(0, 1, matrix, from, 0, to));
    return IGRAPH_SUCCESS;
}

/* MDS layout for a connected graph, with no error checking on the
 * input parameters. The distance matrix will be modified in-place. */
igraph_error_t igraph_i_layout_mds_single(const igraph_t* graph, igraph_matrix_t *res,
                               igraph_matrix_t *dist, igraph_integer_t dim) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t nev = dim;
    igraph_matrix_t vectors;
    igraph_vector_t values, row_means;
    igraph_real_t grand_mean;
    igraph_integer_t i, j, k;
    igraph_eigen_which_t which;

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph too large for eigenvector calculations", IGRAPH_EOVERFLOW);
    }

    if (nev > INT_MAX) {
        IGRAPH_ERROR("Dimensionality too large for eigenvector calculations", IGRAPH_EOVERFLOW);
    }

    /* Handle the trivial cases */
    if (no_of_nodes == 1) {
        IGRAPH_CHECK(igraph_matrix_resize(res, 1, dim));
        igraph_matrix_fill(res, 0);
        return IGRAPH_SUCCESS;
    }
    if (no_of_nodes == 2) {
        IGRAPH_CHECK(igraph_matrix_resize(res, 2, dim));
        igraph_matrix_fill(res, 0);
        for (j = 0; j < dim; j++) {
            MATRIX(*res, 1, j) = 1;
        }
        return IGRAPH_SUCCESS;
    }

    /* Initialize some stuff */
    IGRAPH_VECTOR_INIT_FINALLY(&values, no_of_nodes);
    IGRAPH_CHECK(igraph_matrix_init(&vectors, no_of_nodes, dim));
    IGRAPH_FINALLY(igraph_matrix_destroy, &vectors);

    /* Take the square of the distance matrix */
    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            MATRIX(*dist, i, j) *= MATRIX(*dist, i, j);
        }
    }

    /* Double centering of the distance matrix */
    IGRAPH_VECTOR_INIT_FINALLY(&row_means, no_of_nodes);
    igraph_vector_fill(&values, 1.0 / no_of_nodes);
    IGRAPH_CHECK(igraph_blas_dgemv(0, 1, dist, &values, 0, &row_means));
    grand_mean = igraph_vector_sum(&row_means) / no_of_nodes;
    igraph_matrix_add_constant(dist, grand_mean);
    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            MATRIX(*dist, i, j) -= VECTOR(row_means)[i] + VECTOR(row_means)[j];
            MATRIX(*dist, i, j) *= -0.5;
        }
    }
    igraph_vector_destroy(&row_means);
    IGRAPH_FINALLY_CLEAN(1);

    /* Calculate the top `dim` eigenvectors. */
    which.pos = IGRAPH_EIGEN_LA;
    which.howmany = (int) nev;
    IGRAPH_CHECK(igraph_eigen_matrix_symmetric(/*A=*/ 0, /*sA=*/ 0,
                 /*fun=*/ igraph_i_layout_mds_step,
                 /*n=*/ (int) no_of_nodes, /*extra=*/ dist,
                 /*algorithm=*/ IGRAPH_EIGEN_LAPACK,
                 &which, /*options=*/ 0, /*storage=*/ 0,
                 &values, &vectors));

    /* Calculate and normalize the final coordinates */
    for (j = 0; j < nev; j++) {
        VECTOR(values)[j] = sqrt(fabs(VECTOR(values)[j]));
    }
    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, dim));
    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0, k = nev - 1; j < nev; j++, k--) {
            MATRIX(*res, i, k) = VECTOR(values)[j] * MATRIX(vectors, i, j);
        }
    }

    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_mds
 * \brief Place the vertices on a plane using multidimensional scaling.
 *
 * </para><para>
 * This layout requires a distance matrix, where the intersection of
 * row i and column j specifies the desired distance between vertex i
 * and vertex j. The algorithm will try to place the vertices in a
 * space having a given number of dimensions in a way that approximates
 * the distance relations prescribed in the distance matrix. igraph
 * uses the classical multidimensional scaling by Torgerson; for more
 * details, see Cox &amp; Cox: Multidimensional Scaling (1994), Chapman
 * and Hall, London.
 *
 * </para><para>
 * If the input graph is disconnected, igraph will decompose it
 * first into its subgraphs, lay out the subgraphs one by one
 * using the appropriate submatrices of the distance matrix, and
 * then merge the layouts using \ref igraph_layout_merge_dla.
 * Since \ref igraph_layout_merge_dla works for 2D layouts only,
 * you cannot run the MDS layout on disconnected graphs for
 * more than two dimensions.
 *
 * </para><para>
 * Warning: if the graph is symmetric to the exchange of two vertices
 * (as is the case with leaves of a tree connecting to the same parent),
 * classical multidimensional scaling may assign the same coordinates to
 * these vertices.
 *
 * \param graph A graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized if needed.
 * \param dist The distance matrix. It must be symmetric and this
 *        function does not check whether the matrix is indeed
 *        symmetric. Results are unspecified if you pass a non-symmetric
 *        matrix here. You can set this parameter to null; in this
 *        case, the undirected shortest path lengths between vertices
 *        will be used as distances.
 * \param dim The number of dimensions in the embedding space. For
 *        2D layouts, supply 2 here.
 * \return Error code.
 *
 * Added in version 0.6.
 *
 * </para><para>
 * Time complexity: usually around O(|V|^2 dim).
 */

igraph_error_t igraph_layout_mds(const igraph_t* graph, igraph_matrix_t *res,
                      const igraph_matrix_t *dist, igraph_integer_t dim) {
    igraph_integer_t i, no_of_nodes = igraph_vcount(graph);
    igraph_matrix_t m;
    igraph_bool_t conn;

    RNG_BEGIN();

    /* Check the distance matrix */
    if (dist && (igraph_matrix_nrow(dist) != no_of_nodes ||
                 igraph_matrix_ncol(dist) != no_of_nodes)) {
        IGRAPH_ERROR("invalid distance matrix size", IGRAPH_EINVAL);
    }

    /* Check the number of dimensions */
    if (dim <= 1) {
        IGRAPH_ERROR("dim must be positive", IGRAPH_EINVAL);
    }
    if (no_of_nodes > 0 && dim > no_of_nodes) {
        IGRAPH_ERROR("dim must be less than the number of nodes", IGRAPH_EINVAL);
    }

    /* Copy or obtain the distance matrix */
    if (dist == 0) {
        IGRAPH_MATRIX_INIT_FINALLY(&m, no_of_nodes, no_of_nodes);
        IGRAPH_CHECK(igraph_distances(graph, &m, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL));
    } else {
        IGRAPH_CHECK(igraph_matrix_init_copy(&m, dist));
        IGRAPH_FINALLY(igraph_matrix_destroy, &m);
        /* Make sure that the diagonal contains zeroes only */
        for (i = 0; i < no_of_nodes; i++) {
            MATRIX(m, i, i) = 0.0;
        }
    }

    /* Check whether the graph is connected */
    IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
    if (conn) {
        /* Yes, it is, just do the MDS */
        IGRAPH_CHECK(igraph_i_layout_mds_single(graph, res, &m, dim));
    } else {
        /* The graph is not connected, lay out the components one by one */
        igraph_matrix_list_t layouts;
        igraph_vector_int_t vertex_order;
        igraph_vector_int_t comp;
        igraph_t subgraph;
        igraph_matrix_t layout;
        igraph_matrix_t dist_submatrix;
        igraph_bool_t *seen_vertices;
        igraph_integer_t j, n, processed_vertex_count = 0;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&comp, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_order, no_of_nodes);

        IGRAPH_MATRIX_LIST_INIT_FINALLY(&layouts, 0);
        IGRAPH_MATRIX_INIT_FINALLY(&layout, 0, 0);

        IGRAPH_CHECK(igraph_matrix_init(&dist_submatrix, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &dist_submatrix);

        seen_vertices = IGRAPH_CALLOC(no_of_nodes, igraph_bool_t);
        if (seen_vertices == 0) {
            IGRAPH_ERROR("cannot calculate MDS layout", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, seen_vertices);

        for (i = 0; i < no_of_nodes; i++) {
            if (seen_vertices[i]) {
                continue;
            }

            /* This is a vertex whose component we did not lay out so far */
            IGRAPH_CHECK(igraph_subcomponent(graph, &comp, i, IGRAPH_ALL));
            /* Take the subgraph */
            IGRAPH_CHECK(igraph_induced_subgraph(graph, &subgraph, igraph_vss_vector(&comp),
                                                 IGRAPH_SUBGRAPH_AUTO));
            IGRAPH_FINALLY(igraph_destroy, &subgraph);
            /* Calculate the submatrix of the distances */
            IGRAPH_CHECK(igraph_matrix_select_rows_cols(&m, &dist_submatrix, &comp, &comp));
            /* Lay out the subgraph */
            IGRAPH_CHECK(igraph_i_layout_mds_single(&subgraph, &layout, &dist_submatrix, dim));
            /* Store the layout */
            IGRAPH_CHECK(igraph_matrix_list_push_back_copy(&layouts, &layout));
            /* Free the newly created subgraph */
            igraph_destroy(&subgraph);
            IGRAPH_FINALLY_CLEAN(1);
            /* Mark all the vertices in the component as visited */
            n = igraph_vector_int_size(&comp);
            for (j = 0; j < n; j++) {
                seen_vertices[VECTOR(comp)[j]] = 1;
                VECTOR(vertex_order)[VECTOR(comp)[j]] = processed_vertex_count++;
            }
        }
        /* Merge the layouts - reusing dist_submatrix here */
        IGRAPH_CHECK(igraph_layout_merge_dla(0, &layouts, &dist_submatrix));
        /* Reordering the rows of res to match the original graph */
        IGRAPH_CHECK(igraph_matrix_select_rows(&dist_submatrix, res, &vertex_order));

        igraph_free(seen_vertices);
        igraph_matrix_destroy(&dist_submatrix);
        igraph_matrix_destroy(&layout);
        igraph_matrix_list_destroy(&layouts);
        igraph_vector_int_destroy(&vertex_order);
        igraph_vector_int_destroy(&comp);
        IGRAPH_FINALLY_CLEAN(6);
    }

    RNG_END();

    igraph_matrix_destroy(&m);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
