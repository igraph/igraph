/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#define DIM 10

int check_ev(const igraph_matrix_t *A, const igraph_vector_t *values,
             const igraph_matrix_t *vectors) {

    int i, n = igraph_matrix_nrow(A);
    int ne = igraph_matrix_ncol(vectors);
    igraph_vector_t v, lhs, rhs;

    if (ne != igraph_vector_size(values)) {
        printf("'values' and 'vectors' sizes do not match\n");
        exit(1);
    }

    igraph_vector_init(&lhs, n);
    igraph_vector_init(&rhs, n);

    for (i = 0; i < ne; i++) {
        igraph_vector_view(&v, &MATRIX(*vectors, 0, i), n);
        igraph_blas_dgemv(/*transpose=*/ 0, /*alpha=*/ 1, A, &v,
                                         /*beta=*/ 0, &lhs);
        igraph_vector_update(&rhs, &v);
        igraph_vector_scale(&rhs, VECTOR(*values)[i]);
        if (igraph_vector_maxdifference(&lhs, &rhs) > 1e-10) {
            printf("LHS: ");
            igraph_vector_print(&lhs);
            printf("RHS: ");
            igraph_vector_print(&rhs);
            exit(2);
        }
    }

    igraph_vector_destroy(&rhs);
    igraph_vector_destroy(&lhs);

    return 0;
}

int main() {

    igraph_matrix_t A;
    igraph_vector_t values;
    igraph_matrix_t vectors;
    int i, j;
    igraph_eigen_which_t which;

    igraph_rng_seed(igraph_rng_default(), 42 * 42);

    igraph_matrix_init(&A, DIM, DIM);
    igraph_matrix_init(&vectors, 0, 0);
    igraph_vector_init(&values, 0);

    /* All eigenvalues and eigenvectors */

    for (i = 0; i < DIM; i++) {
        for (j = i; j < DIM; j++) {
            MATRIX(A, i, j) = MATRIX(A, j, i) =
                                  igraph_rng_get_integer(igraph_rng_default(), 1, 10);
        }
    }

    which.pos = IGRAPH_EIGEN_LM;
    which.howmany = 5;
    igraph_eigen_matrix_symmetric(&A, /*sA=*/ 0, /*fun=*/ 0, DIM, /*extra=*/ 0,
                                  IGRAPH_EIGEN_LAPACK, &which, /*options=*/ 0,
                                  /*storage=*/ 0, &values, &vectors);
    igraph_vector_print(&values);
    check_ev(&A, &values, &vectors);

    which.howmany = 8;
    igraph_eigen_matrix_symmetric(&A, /*sA=*/ 0, /*fun=*/ 0, DIM, /*extra=*/ 0,
                                  IGRAPH_EIGEN_LAPACK, &which, /*options=*/ 0,
                                  /*storage=*/ 0, &values, &vectors);
    igraph_vector_print(&values);
    check_ev(&A, &values, &vectors);

    which.pos = IGRAPH_EIGEN_BE;
    which.howmany = 5;
    igraph_eigen_matrix_symmetric(&A, /*sA=*/ 0, /*fun=*/ 0, DIM, /*extra=*/ 0,
                                  IGRAPH_EIGEN_LAPACK, &which, /*options=*/ 0,
                                  /*storage=*/ 0, &values, &vectors);
    igraph_vector_print(&values);
    check_ev(&A, &values, &vectors);

    which.pos = IGRAPH_EIGEN_SM;
    which.howmany = 5;
    igraph_eigen_matrix_symmetric(&A, /*sA=*/ 0, /*fun=*/ 0, DIM, /*extra=*/ 0,
                                  IGRAPH_EIGEN_LAPACK, &which, /*options=*/ 0,
                                  /*storage=*/ 0, &values, &vectors);
    igraph_vector_print(&values);
    check_ev(&A, &values, &vectors);

    igraph_vector_destroy(&values);
    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&A);

    VERIFY_FINALLY_STACK();

    return 0;
}
