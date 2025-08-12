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

#include "qhull/libqhull_r/libqhull_r.h"
#include "qhull/libqhull_r/poly_r.h"

void add_clique(igraph_vector_int_t *destination, igraph_vector_int_t *source) {
    igraph_integer_t num_points = igraph_vector_int_size(source);
    for (igraph_integer_t a = 0; a < num_points - 1; a++) {
        for (igraph_integer_t b = a + 1; b < num_points; b++) {
            igraph_vector_int_push_back(destination, VECTOR(*source)[a]);
            igraph_vector_int_push_back(destination, VECTOR(*source)[b]);
        }
    }
}

int edge_comparator(const void *a, const void *b) {
    igraph_integer_t * A = (igraph_integer_t *)a;
    igraph_integer_t * B = (igraph_integer_t *)b;
    if (A[0] < B[0]) {
        return -1;
    }
    if (A[0] > B[0]) {
        return  1;
    }
    //first are equal
    if (A[1] < B[1]) {
        return -1;
    }
    if (A[1] > B[1]) {
        return  1;
    }

    // second are equal
    return 0;
}

// Simplify an edge list
void simplify_edge_list(igraph_vector_int_t *in, igraph_vector_int_t *out, igraph_integer_t except) {
    igraph_integer_t size = igraph_vector_int_size(in);
    if (size == 0) {
        return;
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

    igraph_vector_int_push_back(out, VECTOR(*in)[0]);
    igraph_vector_int_push_back(out, VECTOR(*in)[1]);
    // nub and filter out point at infinity
    for (igraph_integer_t i = 2; i < size; i += 2) {
        if ((VECTOR(*in)[i] != VECTOR(*in)[i - 2] || VECTOR(*in)[i + 1] != VECTOR(*in)[i - 1]) && VECTOR(*in)[i+1] != except) {
            igraph_vector_int_push_back(out, VECTOR(*in)[i]);
            igraph_vector_int_push_back(out, VECTOR(*in)[i + 1]);
        }
    }
}

igraph_error_t igraph_i_delaunay_edges(igraph_vector_int_t *edges, igraph_matrix_t *points) {
    int curlong, totlong; /* used !qh_NOmem */
    int exitcode;
    igraph_integer_t numpoints = igraph_matrix_nrow(points) ;
    igraph_integer_t dim = igraph_matrix_ncol(points);
    boolT ismalloc = False; // handle memory allocation of points explicitly
    qhT qh_qh;
    qhT *qh = &qh_qh;

    igraph_vector_int_t int_edges;

    igraph_vector_t int_points;

    if (numpoints > INT_MAX) {
        IGRAPH_ERROR("Too many points", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_transpose(points)); // get in row-major order

    IGRAPH_VECTOR_INIT_FINALLY(&int_points, numpoints * dim); //allocate space for extra ponint added due to

    igraph_matrix_copy_to(points, VECTOR(int_points));

    IGRAPH_CHECK(igraph_matrix_transpose(points)); // Transpose = Transpose ^ -1, so this returns the point set to the initial state.

    QHULL_LIB_CHECK; /* Check for compatible library */

    qh_init_A(qh, NULL, NULL, NULL, 0, NULL);  /* sets qh->qhull_command */
    exitcode = setjmp(qh->errexit); /* simple statement for CRAY J916 */
    if (!exitcode) {
        qh->NOerrexit = False;
        qh_option(qh, "delaunay  Qbbound-last QJ", NULL, NULL);
        qh->PROJECTdelaunay = True; // project points to parabola to calculate delaunay
        qh->USEstdout = False;
        qh->DELAUNAY = True;    /* 'd'   */
        qh->SCALElast = True;   /* 'Qbb' */
        qh->KEEPcoplanar = True; /* 'Qc', to keep coplanars in 'p' */
        qh->TRIangulate = True;
        qh_checkflags(qh, qh->qhull_command, "  ");
        qh_initflags(qh, qh->qhull_command);
        qh->ATinfinity = True;
        qh->MERGEvertices = False;
        qh_init_B(qh, VECTOR(int_points), numpoints, dim, ismalloc); // read points and fiddle with them a little.
        qh_qhull(qh);
        qh_triangulate(qh);

        facetT *facet;
        vertexT *vertex, **vertexp;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&int_edges, 0);

        igraph_vector_int_t simplex;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&simplex, dim + 1); // a simplex in n dimensions has n+1 incident vertices.

        igraph_integer_t curr_vert;

        FORALLfacets {
            if (!facet->upperdelaunay) {
                curr_vert = 0;
                FOREACHvertex_(facet->vertices) {
                    VECTOR(simplex)[curr_vert++] = qh_pointid(qh, vertex->point);
                }
                add_clique(&int_edges, &simplex);
            }
        }
        IGRAPH_CHECK(igraph_matrix_transpose(points));
        simplify_edge_list(&int_edges, edges, numpoints + 1);

        igraph_vector_destroy(&int_points);
        igraph_vector_int_destroy(&int_edges);
        igraph_vector_int_destroy(&simplex);
        IGRAPH_FINALLY_CLEAN(3);
    } else {
        qh->NOerrexit = True; /* no more setjmp */
        qh_freeqhull(qh, !qh_ALL);
        qh_memfreeshort(qh, &curlong, &totlong);
        switch (qh->last_errcode) {
            default:
                IGRAPH_ERRORF("Error while computing delaunay triangulation, qhull error code %" IGRAPH_PRId "", IGRAPH_EINVAL, (igraph_integer_t)qh->last_errcode);
        }
    }

    qh->NOerrexit = True; /* no more setjmp */
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_delaunay_triangulation
 * \brief Computes the delaunay triangulation for a spatial point set.
 *
 * This function constructs the delaunay triangulation of a given spatial point set.
 *
 * \param graph A pointer to the graph that will be created.
 * \param points a matrix containing the points that will be used to create the graph.
 *     Each row is a point, dimensionality is inferred from the column count.
 *
 * \return Error code, often a qhull error code, refer to qhull for more detail.
 *
 * */
igraph_error_t igraph_delaunay_triangulation(igraph_t *graph, igraph_matrix_t *points) {
    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    IGRAPH_CHECK(igraph_i_delaunay_edges(&edges, points));

    IGRAPH_CHECK(igraph_create(graph, &edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
