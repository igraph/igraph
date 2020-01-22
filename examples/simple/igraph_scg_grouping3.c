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

    const int nodes = 10;
    igraph_t g;
    igraph_matrix_t V, V3;
    igraph_matrix_complex_t V2;
    igraph_sparsemat_t stochastic, stochasticT;
    igraph_vector_t groups;
    igraph_eigen_which_t which;
    igraph_vector_t p, selcol;

    igraph_tree(&g, nodes, /* children= */ 3, IGRAPH_TREE_UNDIRECTED);

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

