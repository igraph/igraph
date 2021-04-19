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

int main() {

    igraph_t g;
    igraph_matrix_t L, R;
    igraph_sparsemat_t Lsparse, Rsparse;
    igraph_matrix_t adj, V;
    igraph_vector_t groups;
    igraph_eigen_which_t which;

    igraph_matrix_init(&L, 0, 0);
    igraph_matrix_init(&R, 0, 0);
    igraph_matrix_init(&adj, 0, 0);
    igraph_matrix_init(&V, 0, 0);
    igraph_vector_init(&groups, 0);

    igraph_tree(&g, 10, /* children= */ 3, IGRAPH_TREE_UNDIRECTED);

    igraph_get_adjacency(&g, &adj, IGRAPH_GET_ADJACENCY_BOTH, /*eids=*/ 0);

    which.pos = IGRAPH_EIGEN_LM;
    which.howmany = 1;
    igraph_eigen_matrix_symmetric(&adj, /*sparsemat=*/ 0, /*fun=*/ 0,
                                  igraph_vcount(&g), /*extra=*/ 0,
                                  /*algorithm=*/ IGRAPH_EIGEN_LAPACK,
                                  &which, /*options=*/ 0, /*storage=*/ 0,
                                  /*values=*/ 0, &V);

#define SEMI()                              \
    do {                                  \
        igraph_scg_semiprojectors(&groups, IGRAPH_SCG_SYMMETRIC, &L, &R,    \
                                  &Lsparse, &Rsparse, /*p=*/ 0,     \
                                  IGRAPH_SCG_NORM_ROW);         \
    } while(0)

#define PRINTRES()              \
    do {                      \
        printf("----------------------\n");     \
        igraph_matrix_print(&L);            \
        printf("---\n");                \
        igraph_matrix_print(&R);            \
        printf("---\n");                \
        igraph_sparsemat_destroy(&Lsparse);         \
        igraph_sparsemat_destroy(&Rsparse);         \
    } while (0)

    /* -------------- */

    igraph_scg_grouping(&V, &groups, /*intervals=*/ 3,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_SYMMETRIC,
                        IGRAPH_SCG_OPTIMUM, /*p=*/ 0, /*maxiter=*/ 10000);
    SEMI();
    PRINTRES();

    /* -------------- */

    igraph_scg_grouping(&V, &groups, /*intervals=*/ 2,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_SYMMETRIC,
                        IGRAPH_SCG_INTERV_KM, /*p=*/ 0, /*maxiter=*/ 10000);
    SEMI();
    PRINTRES();

    /* -------------- */

    igraph_scg_grouping(&V, &groups, /*intervals=*/ 2,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_SYMMETRIC,
                        IGRAPH_SCG_INTERV, /*p=*/ 0, /*maxiter=*/ 10000);
    SEMI();
    PRINTRES();

    /* -------------- */

    igraph_scg_grouping(&V, &groups, /*(ignored) intervals=*/ 0,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_SYMMETRIC,
                        IGRAPH_SCG_EXACT, /*p=*/ 0, /*maxiter=*/ 10000);
    SEMI();
    PRINTRES();

    /* -------------- */

    igraph_vector_destroy(&groups);
    igraph_matrix_destroy(&L);
    igraph_matrix_destroy(&R);
    igraph_matrix_destroy(&V);
    igraph_matrix_destroy(&adj);
    igraph_destroy(&g);

    return 0;
}

