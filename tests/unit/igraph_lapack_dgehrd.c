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

    int nodes = 10;
    igraph_t tree;
    igraph_matrix_t sto;
    igraph_matrix_t hess;
    igraph_matrix_complex_t evec1, evec2;
    igraph_vector_complex_t eval1, eval2;
    igraph_eigen_which_t which;
    int i;

    igraph_tree(&tree, nodes, /* children= */ 3, IGRAPH_TREE_UNDIRECTED);

    igraph_matrix_init(&sto, nodes, nodes);
    igraph_get_stochastic(&tree, &sto, /*column_wise=*/ 0);
    igraph_matrix_transpose(&sto);

    igraph_matrix_init(&hess, nodes, nodes);
    igraph_lapack_dgehrd(&sto, 1, nodes, &hess);

    igraph_matrix_complex_init(&evec1, 0, 0);
    igraph_vector_complex_init(&eval1, 0);
    which.pos = IGRAPH_EIGEN_ALL;
    igraph_eigen_matrix(&sto, 0, 0, nodes, 0, IGRAPH_EIGEN_LAPACK, &which, 0, 0,
                        &eval1, &evec1);

    igraph_matrix_complex_init(&evec2, 0, 0);
    igraph_vector_complex_init(&eval2, 0);
    igraph_eigen_matrix(&hess, 0, 0, nodes, 0, IGRAPH_EIGEN_LAPACK, &which, 0,
                        0, &eval2, &evec2);

    for (i = 0; i < nodes; i++) {
        igraph_real_t d = igraph_complex_abs(igraph_complex_sub(VECTOR(eval1)[i],
                                             VECTOR(eval2)[i]));
        if (d > 1e-14) {
            printf("Difference: %g\n", d);
            return 1;
        }
    }

    igraph_matrix_complex_destroy(&evec2);
    igraph_vector_complex_destroy(&eval2);

    igraph_matrix_complex_destroy(&evec1);
    igraph_vector_complex_destroy(&eval1);

    igraph_matrix_destroy(&hess);
    igraph_matrix_destroy(&sto);
    igraph_destroy(&tree);

    VERIFY_FINALLY_STACK();

    return 0;
}
