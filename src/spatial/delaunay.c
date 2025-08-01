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



#include "igraph_constructors.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_matrix.h"

#include "igraph_operators.h"
#include "igraph_vector.h"
#include "libqhull_r/io_r.h"
#include "libqhull_r/merge_r.h"
#include "qhull/libqhull_r/libqhull_r.h"
#include "qhull/libqhull_r/poly_r.h"

void add_clique(igraph_vector_int_t *destination, igraph_vector_int_t *source) {
    igraph_integer_t num_points = igraph_vector_int_size(source);
    for (igraph_integer_t a = 0; a < num_points - 1; a++) {
        for (igraph_integer_t b = a + 1; b < num_points; b++) {
            igraph_vector_int_push_back(destination, VECTOR(*source)[a]);
            igraph_vector_int_push_back(destination, VECTOR(*source)[b]);
            printf("%li->%li, ", VECTOR(*source)[a], VECTOR(*source)[b]);
        }
    }
    printf("\n");
}

igraph_error_t igraph_delaunay_triangulation(igraph_t *graph, igraph_matrix_t *points_) {
    int curlong, totlong; /* used !qh_NOmem */
    int exitcode;
    int numpoints = 5;
    int dim = 2;
    //coordT *points;
    boolT ismalloc = False; // handle memory allocation of points explicitly
    qhT qh_qh;
    qhT *qh = &qh_qh;
    igraph_real_t *points;

    IGRAPH_CHECK(igraph_matrix_transpose(points_));

    points = &MATRIX(*points_, 0, 0);

    QHULL_LIB_CHECK; /* Check for compatible library */

    qh_init_A(qh, stdin, stdout, stderr, 0, NULL);  /* sets qh->qhull_command */
    exitcode = setjmp(qh->errexit); /* simple statement for CRAY J916 */
    if (!exitcode) {
        qh->NOerrexit = False;
        qh_option(qh, "delaunay  Qbbound-last", NULL, NULL);
        qh->PROJECTdelaunay = True; // project points to parabola to calculate delaunay
        qh->PRINToptions1st = True;
        qh->DELAUNAY = True;    /* 'd'   */
        qh->SCALElast = True;   /* 'Qbb' */
        qh->KEEPcoplanar = True; /* 'Qc', to keep coplanars in 'p' */
        qh->TRIangulate = True;
        qh_checkflags(qh, qh->qhull_command, "  ");
        qh_initflags(qh, qh->qhull_command);
        //points = qh_readpoints(qh, &numpoints, &dim, &ismalloc); // read points from file
        qh_init_B(qh, points, numpoints, dim, ismalloc);
        qh_qhull(qh);
        for (int i = 0; i < 10; i++) {
            printf("%f, ", points[i]);
        }
        qh_check_output(qh);
        qh_produce_output(qh);
        if (qh->VERIFYoutput && !qh->FORCEoutput && !qh->STOPpoint && !qh->STOPcone) {
            qh_check_points(qh);
        }

        facetT *facet;
        vertexT *vertex, **vertexp;

        igraph_vector_int_t edges;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

        igraph_vector_int_t simplex;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&simplex, 3);

        igraph_integer_t curr_vert;



        FORALLfacets {
            if (!facet->upperdelaunay) {
                printf ("%d", qh_setsize (qh, facet->vertices));
                curr_vert = 0;
                FOREACHvertex_(facet->vertices) VECTOR(simplex)[curr_vert++] = qh_pointid(qh, vertex->point);
                add_clique(&edges, &simplex);
            }
        }
        igraph_create(graph, &edges, numpoints, false);
        igraph_simplify(graph, true, true, NULL);
        igraph_vector_int_destroy(&simplex);
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(2);
    }

    qh->NOerrexit = True; /* no more setjmp */
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);

    if (curlong || totlong)
        qh_fprintf_stderr(7079, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
                          totlong, curlong);

    return IGRAPH_SUCCESS;
}
