/*
   IGraph library.
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

#include "igraph_constructors.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_qsort.h"
#include "igraph_vector.h"

#include "spatial/spatial_internal.h"

#include "qhull/libqhull_r/libqhull_r.h"
#include "qhull/libqhull_r/poly_r.h"

static void add_clique(igraph_vector_int_t *destination, const igraph_vector_int_t *source) {
    igraph_integer_t num_points = igraph_vector_int_size(source);
    for (igraph_integer_t a = 0; a < num_points - 1; a++) {
        for (igraph_integer_t b = a + 1; b < num_points; b++) {
            igraph_vector_int_push_back(destination, VECTOR(*source)[a]);
            igraph_vector_int_push_back(destination, VECTOR(*source)[b]);
        }
    }
}


static int edge_comparator(const void *a, const void *b) {
    igraph_integer_t * A = (igraph_integer_t *)a;
    igraph_integer_t * B = (igraph_integer_t *)b;
    if (A[0] < B[0]) {
        return -1;
    }
    if (A[0] > B[0]) {
        return  1;
    }
    // first are equal
    if (A[1] < B[1]) {
        return -1;
    }
    if (A[1] > B[1]) {
        return  1;
    }

    // second are equal
    return 0;
}

// Simplify an edge list in place
// except is the vertex "at infinity", and should not be included in the edge list.
static igraph_error_t simplify_edge_list(igraph_vector_int_t *in, igraph_integer_t except) {
    igraph_integer_t size = igraph_vector_int_size(in);
    if (size == 0) {
        return IGRAPH_SUCCESS;
    }
    // reorder
    igraph_integer_t temp;
    for (igraph_integer_t i = 0; i < size; i += 2) {
        if (VECTOR(*in)[i] > VECTOR(*in)[i + 1]) {
            temp = VECTOR(*in)[i];
            VECTOR(*in)[i] = VECTOR(*in)[i+1];
            VECTOR(*in)[i+1] = temp;
        }
    }
    // sort
    igraph_qsort(VECTOR(*in), size / 2, 2 * sizeof(igraph_integer_t), &edge_comparator);
    // nub and filter out point at infinity
    igraph_integer_t last_added = 0;
    for (igraph_integer_t i = 2; i < size; i += 2) {
        if ((VECTOR(*in)[i] != VECTOR(*in)[2 * last_added] || VECTOR(*in)[i + 1] != VECTOR(*in)[2*last_added + 1]) && VECTOR(*in)[i+1] != except) {
            last_added += 1;
            VECTOR(*in)[2*last_added]      = VECTOR(*in)[i];
            VECTOR(*in)[2*last_added + 1 ] = VECTOR(*in)[i+1];
        }
    }
    IGRAPH_CHECK(igraph_vector_int_resize(in, 2*last_added+2));
    return IGRAPH_SUCCESS;
}


// helper that can go on the finally stack to free qhull
static void destroy_qhull(qhT *qh) {
    int curlong, totlong;
    qh->NOerrexit = True; /* no more setjmp */
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);
}

// Make a transposed copy with space for added point at infinity.
static igraph_error_t copy_transpose(const igraph_matrix_t *in, igraph_matrix_t *out) {
    igraph_integer_t rows = igraph_matrix_nrow(in), columns = igraph_matrix_ncol(in);

    // add extra space for point at infinity
    IGRAPH_CHECK(igraph_matrix_init(out, columns, rows+1));

    for (igraph_integer_t row = 0; row < rows; row++) {
        for (igraph_integer_t col = 0; col < columns; col++) {
            MATRIX(*out, col, row) = MATRIX(*in, row, col);
        }
    }

    return IGRAPH_SUCCESS;
}


igraph_error_t igraph_i_delaunay_edges(igraph_vector_int_t *edges, const igraph_matrix_t *points) {
    int curlong, totlong; /* used !qh_NOmem */
    int exitcode;
    igraph_integer_t numpoints = igraph_matrix_nrow(points) ;
    igraph_integer_t dim = igraph_matrix_ncol(points);
    boolT ismalloc = False; // handle memory allocation of points explicitly
    qhT qh_qh;
    qhT *qh = &qh_qh;

    igraph_matrix_t int_points;

    if (numpoints > INT_MAX) {
        IGRAPH_ERROR("Too many points for Qhull.", IGRAPH_EOVERFLOW);
    }

    // For internal use, it returns a complete graph.
    // This is not nessecarily a delaunay graph, but for nondegenerate point sets it is.
    if (numpoints <= dim) {

    }

    if (dim == 0) {
        IGRAPH_ERROR("0 dimensional point sets are not supported.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(copy_transpose(points, &int_points));
    IGRAPH_FINALLY(igraph_matrix_destroy, &int_points);

    QHULL_LIB_CHECK; /* Check for compatible library */

    qh_init_A(qh, NULL, NULL, NULL, 0, NULL);  /* sets qh->qhull_command */
    IGRAPH_FINALLY(destroy_qhull, qh);
    exitcode = setjmp(qh->errexit); /* simple statement for CRAY J916 */
    if (!exitcode) {
        qh->NOerrexit = False;
        qh_option(qh, "delaunay  Qbbound-last QJ", NULL, NULL);
        qh->PROJECTdelaunay = True; // project points to parabola to calculate delaunay
        qh->USEstdout = False;
        qh->DELAUNAY = True;    /* 'd'   */
        qh->SCALElast = True;   /* 'Qbb' */
        qh->KEEPcoplanar = True; /* 'Qc', to keep coplanars in 'p' */
        qh->TRIangulate = True;  // Set flag to split non-simplicial points into
        qh_checkflags(qh, qh->qhull_command, "  ");
        qh_initflags(qh, qh->qhull_command);
        qh->ATinfinity = True; // required for cocircular points, should not mess anything else up.
        qh->MERGEvertices = False; // Do not merge identical vertices
        qh_init_B(qh, &MATRIX(int_points,0,0), numpoints, dim, ismalloc); // read points and project them to parabola.
        qh_qhull(qh); // Do the triangulation
        qh_triangulate(qh); // In some cases, the triangulate flag is not enough, this guarantees that everyhing is simplicial

        facetT *facet;
        vertexT *vertex, **vertexp;

        igraph_vector_int_t simplex;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&simplex, dim + 1); // a simplex in n dimensions has n+1 incident vertices.

        igraph_integer_t curr_vert;

        FORALLfacets {
            if (!facet->upperdelaunay) {
                curr_vert = 0;
                FOREACHvertex_(facet->vertices) {
                    VECTOR(simplex)[curr_vert++] = qh_pointid(qh, vertex->point);
                }
                add_clique(edges, &simplex);
            }
        }
        IGRAPH_CHECK(simplify_edge_list(edges, numpoints + 1));

        igraph_matrix_destroy(&int_points);
        igraph_vector_int_destroy(&simplex);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        switch (qh->last_errcode) {
            default:
                IGRAPH_ERRORF("Error while computing delaunay triangulation, Qhull error code %d.", IGRAPH_EINVAL, qh->last_errcode);
        }
    }

    destroy_qhull(qh);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_delaunay_triangulation
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
 * \param graph A pointer to the graph that will be created.
 * \param points A matrix containing the points that will be used to create the graph.
 *     Each row is a point, dimensionality is inferred from the column count.
 *
 * \return Error code.
 *
 */
igraph_error_t igraph_delaunay_graph(igraph_t *graph, const igraph_matrix_t *points) {
    igraph_vector_int_t edges;

    igraph_integer_t num_points = igraph_matrix_nrow(points);

    if (num_points < 2) {
        igraph_full(graph, num_points, false, false);
        return IGRAPH_SUCCESS;
    }

    if (num_points <= igraph_matrix_ncol(points)) {
        IGRAPH_ERROR("General point sets without enough points to create a simplex (dimension + 1) are not supported.", IGRAPH_UNIMPLEMENTED);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_i_delaunay_edges(&edges, points));
    IGRAPH_CHECK(igraph_create(graph, &edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
