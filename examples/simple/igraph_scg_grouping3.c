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

    const int nodes = 67;
    igraph_t g;
    igraph_matrix_t V, V3;
    igraph_matrix_complex_t V2;
    igraph_sparsemat_t stochastic, stochasticT;
    igraph_vector_t groups;
    igraph_eigen_which_t which;
    igraph_vector_t p, selcol;

    /* This is a tree with no non-trivial automorphisms */
    igraph_small(&g, nodes, IGRAPH_UNDIRECTED,
                     16, 59, 11, 16, 11, 22, 22, 60, 31, 60, 12, 31, 12, 20, 20, 47, 18,
                     47, 18, 23, 23, 24, 24, 38, 18, 21, 21, 56, 56, 65, 34, 65, 34, 64,
                     54, 64, 22, 25, 56, 63, 28, 63, 37, 63, 1, 47, 1, 5, 50, 65, 2, 50,
                     23, 53, 11, 32, 0, 32, 0, 48, 8, 48, 0, 27, 40, 50, 40, 41, 3, 5, 3,
                     51, 52, 60, 18, 55, 9, 40, 16, 42, 6, 21, 6, 14, 14, 43, 50, 66, 30,
                     38, 30, 62, 13, 14, 7, 37, 42, 44, 41, 49, 29, 60, 4, 29, 61, 65, 43,
                     46, 35, 65, 35, 45, 36, 43, 19, 36, 10, 19, 1, 58, 12, 39, 6, 17, 26,
                     44, 25, 57, 15, 66, 33, 36,
                     -1);

    igraph_matrix_complex_init(&V2, 0, 0);
    igraph_matrix_init(&V, 0, 0);
    igraph_matrix_init(&V3, 0, 0);
    igraph_vector_init(&groups, 0);
    igraph_vector_init(&p, 0);
    igraph_vector_init(&selcol, 1);

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_get_stochastic_sparsemat(&g, &stochastic, /*column-wise=*/ 0);
    igraph_sparsemat_transpose(&stochastic, &stochasticT, /*values=*/ 1);

    which.pos = IGRAPH_EIGEN_LR;
    which.howmany = 1;

    igraph_eigen_matrix(/*matrix=*/ 0, &stochasticT, /*fun=*/ 0, nodes,
                                    /*extra=*/ 0, /*1algorithm=*/ IGRAPH_EIGEN_LAPACK,
                                    &which, /*options=*/ 0, /*storage=*/ 0,
                                    /*values=*/ 0, &V2);
    igraph_matrix_complex_real(&V2, &V);

    /* `p' is always the eigenvector corresponding to the 1-eigenvalue */
    igraph_matrix_get_col(&V, &p, 0);
    igraph_vector_print(&p);

    which.howmany = 3;
    igraph_eigen_matrix(/*matrix=*/ 0, &stochastic, /*fun=*/ 0, nodes,
                                    /*extra=*/ 0, /*algorithm=*/ IGRAPH_EIGEN_LAPACK,
                                    &which, /*options=*/ 0, /*storage=*/ 0,
                                    /*values=*/ 0, &V2);
    igraph_matrix_complex_real(&V2, &V3);
    VECTOR(selcol)[0] = 2;
    igraph_matrix_select_cols(&V3, &V, &selcol);

    /* ------------ */

    igraph_scg_grouping(&V, &groups, /*intervals=*/ 3,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_STOCHASTIC,
                        IGRAPH_SCG_OPTIMUM, &p, /*maxiter=*/ 10000);
    igraph_vector_print(&groups);

    /* ------------ */

    igraph_scg_grouping(&V, &groups, /*intervals=*/ 3,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_STOCHASTIC,
                        IGRAPH_SCG_INTERV_KM, &p, /*maxiter=*/ 10000);
    igraph_vector_print(&groups);

    /* ------------ */

    igraph_scg_grouping(&V, &groups, /*intervals=*/ 3,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_STOCHASTIC,
                        IGRAPH_SCG_INTERV, &p, /*maxiter=*/ 10000);
    igraph_vector_print(&groups);

    /* ------------ */

    igraph_scg_grouping(&V, &groups, /*(ignored) intervals=*/ 0,
                        /*intervals_vector=*/ 0, IGRAPH_SCG_STOCHASTIC,
                        IGRAPH_SCG_EXACT, &p, /*maxiter=*/ 10000);
    igraph_vector_print(&groups);

    /* ------------ */

    igraph_vector_destroy(&p);
    igraph_vector_destroy(&selcol);
    igraph_vector_destroy(&groups);
    igraph_matrix_destroy(&V);
    igraph_matrix_destroy(&V3);
    igraph_matrix_complex_destroy(&V2);
    igraph_sparsemat_destroy(&stochasticT);
    igraph_sparsemat_destroy(&stochastic);
    igraph_destroy(&g);

    return 0;
}

