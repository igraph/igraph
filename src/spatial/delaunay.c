/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_spatial.h"

#include "igraph_bitset.h"
#include "igraph_constructors.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_memory.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "internal/utils.h"
#include "spatial/spatial_internal.h"

#include "qhull/libqhull_r/libqhull_r.h"

/**
 * Raises an error if a spatial point set is invalid.
 * The coordinate matrix must have at least one column and must not
 * contain NaN or infinities.
 *
 * \param points Matrix, each row is a spatial point.
 * \return Error code.
 */
igraph_error_t igraph_i_check_spatial_points(const igraph_matrix_t *points) {
    const igraph_int_t dim = igraph_matrix_ncol(points);
    const igraph_int_t n = igraph_matrix_nrow(points);

    /* Special case: we allow zero columns when there are zero rows, i.e. no points.
     * Some languages cannot represent size-zero matrices and length-zero point
     * lists may translate as a 0-by-0 matrix to igraph. */
    if (dim == 0 && n > 0) {
        IGRAPH_ERROR("Point sets must not be zero-dimensional.", IGRAPH_EINVAL);
    }

    if (!igraph_vector_is_all_finite(&points->data)) {
        IGRAPH_ERROR("Coordinates must not be NaN or infinite.", IGRAPH_EINVAL);
    }

    return IGRAPH_SUCCESS;
}

// Append an undirected clique of the indices in destination to source.
// Assumes that source is an initialized vector.
static igraph_error_t add_clique(igraph_vector_int_t *destination, const igraph_vector_int_t *source) {
    igraph_int_t num_points = igraph_vector_int_size(source);
    for (igraph_int_t a = 0; a < num_points - 1; a++) {
        for (igraph_int_t b = a + 1; b < num_points; b++) {
            IGRAPH_CHECK(igraph_vector_int_push_back(destination, VECTOR(*source)[a]));
            IGRAPH_CHECK(igraph_vector_int_push_back(destination, VECTOR(*source)[b]));
        }
    }
    return IGRAPH_SUCCESS;
}


// Helper that can go on the finally stack to free qhT Qhull data structure in case of an error.
static void destroy_qhull(qhT *qh) {
    int curlong, totlong;
    qh->NOerrexit = True; /* no more setjmp */
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);
}


/* In the 1D case we simply connect points in sorted order.
 * This function assumes that there is at least one point. */
static igraph_error_t delaunay_edges_1d(igraph_vector_int_t *edges, const igraph_matrix_t *points) {
    const igraph_int_t numpoints = igraph_matrix_nrow(points);
    const igraph_vector_t coords = igraph_vector_view(&MATRIX(*points, 0, 0), numpoints);
    igraph_vector_int_t order;

    IGRAPH_ASSERT(igraph_matrix_ncol(points) == 1);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, numpoints);

    IGRAPH_CHECK(igraph_vector_sort_ind(&coords, &order, IGRAPH_ASCENDING));

    IGRAPH_CHECK(igraph_vector_int_resize(edges, 2*(numpoints-1)));

    for (igraph_int_t i=0; i < numpoints-1; i++) {
        igraph_int_t from = VECTOR(order)[i];
        igraph_int_t to   = VECTOR(order)[i+1];
        VECTOR(*edges)[2*i] = from;
        VECTOR(*edges)[2*i + 1] = to;
        if (VECTOR(coords)[from] == VECTOR(coords)[to]) {
            IGRAPH_ERROR("Duplicate points for Delaunay triangulation.", IGRAPH_EINVAL);
        }
    }

    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_delaunay_edges(igraph_vector_int_t *edges, const igraph_matrix_t *points) {
    const igraph_int_t numpoints = igraph_matrix_nrow(points) ;
    const igraph_int_t dim = igraph_matrix_ncol(points);
    int exitcode;
    qhT qh_qh; /* Qhull's data structure. First argument of most Qhull calls. */
    qhT *qh = &qh_qh; /* Convenience pointer. */

    /* Error checks */

    /* Validate point set. */
    IGRAPH_CHECK(igraph_i_check_spatial_points(points));

    /* No edges for one or zero points. */
    if (numpoints <= 1) {
        igraph_vector_int_clear(edges);
        return IGRAPH_SUCCESS;
    }

    /* Qhull does not support the 1D case. */
    if (dim == 1) {
        return delaunay_edges_1d(edges, points);
    }

    if (dim >= numpoints) {
        IGRAPH_ERRORF("Not enough points to create simplex, need at least %" IGRAPH_PRId ".", IGRAPH_EINVAL, dim);
    }

    /* Prevent overflow in igraph_int_t -> int conversions below.
     * Note that Qhull will likely already fail for a much smaller point count. */
    if (numpoints > INT_MAX) {
        IGRAPH_ERROR("Too many points for Qhull.", IGRAPH_EOVERFLOW);
    }

    /* Prepare point set in row-major format for Qhull */

    coordT *qhull_points = IGRAPH_CALLOC(dim*numpoints, igraph_real_t);
    IGRAPH_CHECK_OOM(qhull_points, "Insufficient memory for constructing Delaunay graph.");
    IGRAPH_FINALLY(igraph_free, qhull_points);
    igraph_matrix_copy_to(points, qhull_points, IGRAPH_ROW_MAJOR);

    /* Call Qhull.
     *
     * This is mainly based on qdelaunay/qdelaun_r.c and qh_new_qhull() in user_r.c.
     *
     * Note that output routines in userprintf_r.c are patched to not print
     * when passing NULL as a file pointer, which is what we do here. */

    /* Check for compatible library. Not technically necessary, as igraph
     * vendors Qhull due to the need to override output routines. Would
     * become necessary if linking to an external Qhull. */
    QHULL_LIB_CHECK

    /* Initializes qh, sets qh->qhull_command. */
    qh_init_A(qh, NULL, NULL, NULL, 0, NULL);
    IGRAPH_FINALLY(destroy_qhull, qh);

    exitcode = setjmp(qh->errexit);
    if (!exitcode) {
        qh->NOerrexit = False;

        /* qh_option() does not change the operation of Qhull. It simply records options to be
         * output with error messages. Here we manually keep it in sync with the settings below. */
        qh_option(qh, "delaunay  Qz-infinity-point Q3-no-merge-vertices", NULL, NULL);

        qh->PROJECTdelaunay = True; // project points to parabola to calculate delaunay triangulation
        qh->DELAUNAY = True; // 'd'
        qh->ATinfinity = True; // 'Qz', required for cocircular points
        qh->MERGEvertices = False; // 'Q3', do not merge identical vertices

        qh_initflags(qh, qh->qhull_command);
        qh_init_B(qh, qhull_points, numpoints, dim, /*ismalloc=*/ False); // read points and project them to parabola.
        qh_qhull(qh); // do the triangulation
        qh_triangulate(qh); // this guarantees that everything is simplicial

        igraph_vector_int_t simplex;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&simplex, dim + 1); // a simplex in n dimensions has n+1 incident vertices.

        facetT *facet; // required for FORALLfacets
        vertexT *vertex, **vertexp; // required for FOREACHvertex_

        FORALLfacets {
            if (!facet->upperdelaunay) {
                igraph_int_t curr_vert = 0;
                FOREACHvertex_(facet->vertices) {
                    VECTOR(simplex)[curr_vert++] = qh_pointid(qh, vertex->point);
                }
                IGRAPH_CHECK(add_clique(edges, &simplex));
            }
        }
        igraph_i_simplify_edge_list(edges, true, true, false);

        /* Check if there are any points/vertices that do not appear in the edge list.
         * This happens when there are duplicate points, as Qhull ignores one of them.
         * We raise an error when there are duplicates. */
        igraph_int_t edges_size = igraph_vector_int_size(edges);
        igraph_bitset_t in_edge_list;

        IGRAPH_BITSET_INIT_FINALLY(&in_edge_list, numpoints);
        for (igraph_int_t i = 0; i < edges_size; i++) {
            IGRAPH_BIT_SET(in_edge_list, VECTOR(*edges)[i]);
        }
        if (igraph_bitset_is_any_zero(&in_edge_list)) {
            IGRAPH_ERROR("Duplicate points for Delaunay triangulation.", IGRAPH_EINVAL);
        }

        igraph_bitset_destroy(&in_edge_list);
        igraph_vector_int_destroy(&simplex);
        destroy_qhull(qh);
        igraph_free(qhull_points);
        IGRAPH_FINALLY_CLEAN(4);
    } else {
        switch (qh->last_errcode) {
        /* TODO: More specific error descriptions for common Qhull errors. */
        default:
            /* TODO: Report Qhull error text? */
            IGRAPH_ERRORF("Error while computing Delaunay triangulation, Qhull error code %d.", IGRAPH_EINVAL, qh->last_errcode);
        }
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_delaunay_graph
 * \brief Computes the Delaunay graph of a spatial point set.
 *
 * \experimental
 *
 * This function constructs the graph corresponding to the Delaunay triangulation
 * of an n-dimensional spatial point set.
 *
 * </para><para>
 * The current implementation uses Qhull.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Barber, C. Bradford, David P. Dobkin, and Hannu Huhdanpaa.
 * The Quickhull Algorithm for Convex Hulls.
 * ACM Transactions on Mathematical Software 22, no. 4 (1996): 469–83.
 * https://doi.org/10.1145/235815.235821.
 *
 * \param graph A pointer to the graph that will be created.
 * \param points A matrix containing the points that will be used to create the graph.
 *     Each row is a point, dimensionality is inferred from the column count.
 *     There must not be duplicate points.
 * \return Error code.
 *
 * Time complexity: According to Theorem 3.2 in the Qhull paper,
 * O(n log n) for d &lt;= 3 and O(n^floor(d/2) / floor(d/2)!) where
 * n is the number of points and d is the dimensionality of the point set.
 */
igraph_error_t igraph_delaunay_graph(igraph_t *graph, const igraph_matrix_t *points) {
    igraph_vector_int_t edges;
    const igraph_int_t numpoints = igraph_matrix_nrow(points);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_i_delaunay_edges(&edges, points));
    IGRAPH_CHECK(igraph_create(graph, &edges, numpoints, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
