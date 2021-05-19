/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

#include <igraph.h>

#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_vector_t ev;
    igraph_t scg_graph;
    igraph_matrix_t scg_matrix;
    igraph_sparsemat_t scg_sparsemat;
    igraph_matrix_t L, R;
    igraph_sparsemat_t Lsparse, Rsparse;
    igraph_vector_t groups;
    igraph_vector_complex_t eval;
    igraph_matrix_complex_t evec;

    igraph_tree(&g, 10, /* children= */ 3, IGRAPH_TREE_UNDIRECTED);

    igraph_vector_init(&ev, 1);
    igraph_matrix_init(&L, 0, 0);
    igraph_matrix_init(&R, 0, 0);
    igraph_matrix_init(&scg_matrix, 0, 0);
    igraph_vector_init(&groups, 0);
    igraph_vector_complex_init(&eval, 0);
    igraph_matrix_complex_init(&evec, 0, 0);

#define CALLLAP() do {                           \
        igraph_vector_resize(&groups, 0);                    \
        igraph_vector_complex_resize(&eval, 0);              \
        igraph_matrix_complex_resize(&evec, 0, 0);               \
        igraph_scg_laplacian(&g, /*matrix=*/ 0, /*sparsemat=*/ 0, &ev,   \
                             /* intervals= */ 2, /* intervals_vector= */ 0, \
                             /* algorithm= */ IGRAPH_SCG_EXACT,     \
                             IGRAPH_SCG_NORM_ROW,               \
                             IGRAPH_SCG_DIRECTION_DEFAULT, &eval, &evec,    \
                             &groups, /* use_arpack= */ 0,          \
                             /* maxiter= */ 0, &scg_graph, &scg_matrix, \
                             &scg_sparsemat, &L, &R,            \
                             &Lsparse, &Rsparse);               \
    } while (0)

#define PRINTRES()                      \
    do {                              \
        printf("--------------------------------\n");       \
        igraph_vector_print(&groups);               \
        printf("---\n");                        \
        igraph_vector_complex_print(&eval);             \
        print_matrix_complex_first_row_positive(&evec);             \
        printf("---\n");                        \
        igraph_write_graph_edgelist(&scg_graph, stdout);        \
        printf("---\n");                        \
        igraph_sparsemat_print(&scg_sparsemat, stdout);     \
        printf("---\n");                        \
        igraph_sparsemat_print(&Lsparse, stdout);           \
        printf("---\n");                        \
        igraph_sparsemat_print(&Rsparse, stdout);           \
        printf("---\n");                        \
    } while (0)

    VECTOR(ev)[0] = 1;
    CALLLAP();
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    VECTOR(ev)[0] = 3;
    CALLLAP();
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    igraph_vector_resize(&ev, 2);
    VECTOR(ev)[0] = 1;
    VECTOR(ev)[1] = 3;
    CALLLAP();
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    igraph_matrix_complex_destroy(&evec);
    igraph_vector_complex_destroy(&eval);
    igraph_vector_destroy(&groups);
    igraph_matrix_destroy(&scg_matrix);
    igraph_matrix_destroy(&L);
    igraph_matrix_destroy(&R);
    igraph_vector_destroy(&ev);
    igraph_destroy(&g);

    /* -------------------------------------------------------------------- */

    VERIFY_FINALLY_STACK();

    return 0;
}
