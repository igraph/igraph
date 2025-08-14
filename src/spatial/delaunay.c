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

#include "igraph_bitset.h"
#include "igraph_spatial.h"

#include "igraph_constructors.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_memory.h"
#include "igraph_qsort.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "spatial/spatial_internal.h"

#include "qhull/libqhull_r/libqhull_r.h"
#include "qhull/libqhull_r/poly_r.h"

// Append an undirected clique of the indices in destination to source.
// Assumes that source is an initialized vector
static igraph_error_t add_clique(igraph_vector_int_t *destination, const igraph_vector_int_t *source) {
    igraph_integer_t num_points = igraph_vector_int_size(source);
    for (igraph_integer_t a = 0; a < num_points - 1; a++) {
        for (igraph_integer_t b = a + 1; b < num_points; b++) {
            IGRAPH_CHECK(igraph_vector_int_push_back(destination, VECTOR(*source)[a]));
            IGRAPH_CHECK(igraph_vector_int_push_back(destination, VECTOR(*source)[b]));
        }
    }
    return IGRAPH_SUCCESS;
}

// Lexicographic comparator, used for qsort in simplify_edge_list
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

    // second are equal, so the edges must be equal.
    return 0;
}

// Simplify an edge list in place
// except is the vertex "at infinity", and should not be included in the edge list.
static igraph_error_t simplify_edge_list(igraph_vector_int_t *in, igraph_bool_t self_loops, igraph_bool_t multi_edges, igraph_bool_t directed) {
    igraph_integer_t size = igraph_vector_int_size(in);
    if (size == 0) {
        return IGRAPH_SUCCESS;
    }
    // reorder
    if (!directed) {
        igraph_integer_t temp;
        for (igraph_integer_t i = 0; i < size; i += 2) {
            if (VECTOR(*in)[i] > VECTOR(*in)[i + 1]) {
                temp = VECTOR(*in)[i];
                VECTOR(*in)[i] = VECTOR(*in)[i+1];
                VECTOR(*in)[i+1] = temp;
            }
        }
    }

    // sort, not needed if multi edges are allowed
    if (!multi_edges) {
        igraph_qsort(VECTOR(*in), size / 2, 2 * sizeof(igraph_integer_t), &edge_comparator);
    }

    // remove duplicates
    igraph_integer_t last_added = 0;

    igraph_integer_t skipped = 0;

    if (!self_loops) {
        for (igraph_integer_t i = 0; i < size; i += 2) {
            if (VECTOR(*in)[i] == VECTOR(*in)[i+1]) {
                skipped += 1;
            } else {
                break;
            }
        }
    }
    for (igraph_integer_t i = skipped * 2 + 2 ; i < size; i += 2) {
        if ( !(
            (!multi_edges && VECTOR(*in)[i] == VECTOR(*in)[2 * last_added] && VECTOR(*in)[i + 1] == VECTOR(*in)[2*last_added + 1])
            || (!self_loops && VECTOR(*in)[i] == VECTOR(*in)[i+1]))) {
            last_added += 1;
            VECTOR(*in)[2*last_added]      = VECTOR(*in)[i];
            VECTOR(*in)[2*last_added + 1 ] = VECTOR(*in)[i+1];
        }
    }

    IGRAPH_CHECK(igraph_vector_int_resize(in, 2*last_added+2));
    return IGRAPH_SUCCESS;
}


// helper that can go on the finally stack to free qhull in case of an error
static void destroy_qhull(qhT *qh) {
    int curlong, totlong;
    qh->NOerrexit = True; /* no more setjmp */
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);
}


/* In the 1D case we simply connect points in sorted order.
 * This function assumes that there is at least one point. */
static igraph_error_t delaunay_edges_1d(igraph_vector_int_t *edges, const igraph_matrix_t *points) {
    const igraph_integer_t numpoints = igraph_matrix_nrow(points);
    igraph_vector_t coords;
    igraph_vector_int_t order;

    IGRAPH_ASSERT(igraph_matrix_ncol(points) == 1);

    igraph_vector_view(&coords, &MATRIX(*points, 0, 0), numpoints);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, numpoints);

    IGRAPH_CHECK(igraph_vector_sort_ind(&coords, &order, IGRAPH_ASCENDING));

    IGRAPH_CHECK(igraph_vector_int_resize(edges, 2*(numpoints-1)));

    for (igraph_integer_t i=0; i < numpoints-1; i++) {
        igraph_integer_t from = VECTOR(order)[i];
        igraph_integer_t to   = VECTOR(order)[i+1];
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
    const igraph_integer_t numpoints = igraph_matrix_nrow(points) ;
    const igraph_integer_t dim = igraph_matrix_ncol(points);
    int exitcode;
    boolT ismalloc = False; // handle memory allocation of points explicitly
    qhT qh_qh;
    qhT *qh = &qh_qh;
    coordT *qhull_points;

    /* The point matrix must have at least one column, unless it has zero rows. */
    if (dim == 0 && numpoints > 0) {
        IGRAPH_ERROR("0 dimensional point sets are not supported.", IGRAPH_EINVAL);
    }

    /* No edges for one or zero points. */
    if (numpoints <= 1) {
        igraph_vector_int_clear(edges);
        return IGRAPH_SUCCESS;
    }

    /* Qhull does not support the 1D case. */
    if (dim == 1) {
        return delaunay_edges_1d(edges, points);
    }

    /* Prevent overflow in igraph_integer_t -> int conversions below.
     * Note that Qhull will likely already fail for a much smaller point count. */
    if (numpoints > INT_MAX) {
        IGRAPH_ERROR("Too many points for Qhull.", IGRAPH_EOVERFLOW);
    }

    qhull_points = IGRAPH_CALLOC(dim*numpoints, igraph_real_t);
    IGRAPH_CHECK_OOM(qhull_points, "Insufficient memory for constructing Delaunay graph.");
    IGRAPH_FINALLY(igraph_free, qhull_points);
    igraph_matrix_copy_to(points, qhull_points, IGRAPH_ROW_MAJOR);

    QHULL_LIB_CHECK; /* Check for compatible library */

    qh_init_A(qh, NULL, NULL, NULL, 0, NULL);  /* sets qh->qhull_command */
    IGRAPH_FINALLY(destroy_qhull, qh);

    exitcode = setjmp(qh->errexit); /* simple statement for CRAY J916 */
    if (!exitcode) {
        // flag setting
        qh->NOerrexit = False;
        qh_option(qh, "delaunay  Qbbound-last", NULL, NULL);
        qh->PROJECTdelaunay = True; // project points to parabola to calculate delaunay triangulation
        qh->DELAUNAY = True;    /* 'd'   */
        qh->ATinfinity = True; // required for cocircular points, should not mess anything else up.
        qh->MERGEvertices = True; // Do not merge identical vertices
        qh_initflags(qh, qh->qhull_command);
        qh_init_B(qh, qhull_points, numpoints, dim, ismalloc); // read points and project them to parabola.
        qh_qhull(qh); // Do the triangulation
        qh_triangulate(qh); // this guarantees that everything is simplicial


        igraph_vector_int_t simplex;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&simplex, dim + 1); // a simplex in n dimensions has n+1 incident vertices.

        igraph_integer_t curr_vert;

        facetT *facet; // required for FORALLfacets
        vertexT *vertex, **vertexp; // required for FOREACHvertex_

        FORALLfacets {
            if (!facet->upperdelaunay) {
                curr_vert = 0;
                FOREACHvertex_(facet->vertices) {
                    VECTOR(simplex)[curr_vert++] = qh_pointid(qh, vertex->point);
                }
                IGRAPH_CHECK(add_clique(edges, &simplex));
            }
        }
        IGRAPH_CHECK(simplify_edge_list(edges, false, false, false));


        igraph_integer_t edge_size = igraph_vector_int_size(edges);
        igraph_bitset_t connected_verts;

        IGRAPH_BITSET_INIT_FINALLY(&connected_verts, numpoints);
        for (igraph_integer_t i = 0; i < edge_size; i++) {
            IGRAPH_BIT_SET(connected_verts, VECTOR(*edges)[i]);
        }
        if (igraph_bitset_is_any_zero(&connected_verts)) {
            IGRAPH_ERROR("Duplicate points for Delaunay triangulation.", IGRAPH_EINVAL);
        }
        igraph_bitset_destroy(&connected_verts);
        igraph_vector_int_destroy(&simplex);
        destroy_qhull(qh);
        igraph_free(qhull_points);
        IGRAPH_FINALLY_CLEAN(4);
    } else {
        switch (qh->last_errcode) {
            default:
                IGRAPH_ERRORF("Error while computing delaunay triangulation, Qhull error code %d.", IGRAPH_EINVAL, qh->last_errcode);
        }
    }

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
    const igraph_integer_t numpoints = igraph_matrix_nrow(points);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_i_delaunay_edges(&edges, points));
    IGRAPH_CHECK(igraph_create(graph, &edges, numpoints, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
