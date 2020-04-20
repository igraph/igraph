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
#include <stdio.h>

#define DIM 10

void igraph_print_warning(const char *reason, const char *file,
                          int line, int igraph_errno) {
    printf("Warning: %s\n", reason);
}

int main() {

    igraph_matrix_t A, B, RHS;
    int info;
    int i, j;

    /* Identity matrix, you have to start somewhere */

    igraph_matrix_init(&A, DIM, DIM);
    igraph_matrix_init(&B, DIM, 1);
    for (i = 0; i < DIM; i++) {
        MATRIX(A, i, i) = 1.0;
        MATRIX(B, i, 0) = i + 1;
    }

    igraph_matrix_copy(&RHS, &B);
    igraph_lapack_dgesv(&A, /*ipiv=*/ 0, &RHS, &info);

    if (info != 0) {
        return 1;
    }
    if (!igraph_matrix_all_e(&B, &RHS)) {
        return 2;
    }

    igraph_matrix_destroy(&A);
    igraph_matrix_destroy(&B);
    igraph_matrix_destroy(&RHS);

    /* Diagonal matrix */

    igraph_matrix_init(&A, DIM, DIM);
    igraph_matrix_init(&RHS, DIM, 1);
    for (i = 0; i < DIM; i++) {
        MATRIX(A, i, i) = i + 1;
        MATRIX(RHS, i, 0) = i + 1;
    }

    igraph_lapack_dgesv(&A, /*ipiv=*/ 0, &RHS, &info);

    if (info != 0) {
        return 3;
    }
    for (i = 0; i < DIM; i++) {
        if (MATRIX(RHS, i, 0) != 1.0) {
            return 4;
        }
    }

    igraph_matrix_destroy(&A);
    igraph_matrix_destroy(&RHS);

    /* A general matrix */

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_matrix_init(&A, DIM, DIM);
    igraph_matrix_init(&B, DIM, 1);
    igraph_matrix_init(&RHS, DIM, 1);
    for (i = 0; i < DIM; i++) {
        int j;
        MATRIX(B, i, 0) = igraph_rng_get_integer(igraph_rng_default(), 1, 10);
        for (j = 0; j < DIM; j++) {
            MATRIX(A, i, j) = igraph_rng_get_integer(igraph_rng_default(), 1, 10);
        }
    }
    igraph_blas_dgemv_array(/*transpose=*/ 0, /*alpha=*/ 1.0, /*a=*/ &A,
                                           /*x-*/ &MATRIX(B, 0, 0), /*beta=*/ 0,
                                           /*y=*/ &MATRIX(RHS, 0, 0));

    igraph_lapack_dgesv(&A, /*ipiv=*/ 0, &RHS, &info);
    if (info != 0) {
        return 5;
    }
    for (i = 0; i < DIM; i++) {
        if (fabs(MATRIX(B, i, 0) - MATRIX(RHS, i, 0)) > 1e-13) {
            return 6;
        }
    }

    igraph_matrix_destroy(&A);
    igraph_matrix_destroy(&B);
    igraph_matrix_destroy(&RHS);

    /* A singular matrix */

    igraph_matrix_init(&A, DIM, DIM);
    igraph_matrix_init(&B, DIM, 1);
    igraph_matrix_init(&RHS, DIM, 1);
    for (i = 0; i < DIM; i++) {
        MATRIX(B, i, 0) = igraph_rng_get_integer(igraph_rng_default(), 1, 10);
        for (j = 0; j < DIM; j++) {
            MATRIX(A, i, j) = i == j ? 1 : 0;
        }
    }
    for (i = 0; i < DIM; i++) {
        MATRIX(A, DIM - 1, i) = MATRIX(A, 0, i);
    }

    igraph_blas_dgemv_array(/*transpose=*/ 0, /*alpha=*/ 1.0, /*a=*/ &A,
                                           /*x-*/ &MATRIX(B, 0, 0), /*beta=*/ 0,
                                           /*y=*/ &MATRIX(RHS, 0, 0));

    igraph_set_warning_handler(igraph_print_warning);
    igraph_lapack_dgesv(&A, /*ipiv=*/ 0, &RHS, &info);
    if (info != 10) {
        printf("LAPACK returned info = %d, should have been 10", info);
        return 7;
    }

    igraph_matrix_destroy(&A);
    igraph_matrix_destroy(&B);
    igraph_matrix_destroy(&RHS);

    return 0;
}
