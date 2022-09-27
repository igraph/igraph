/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#define NCOMPLEX  /* to make it compile with MSVC on Windows */

#include <igraph.h>

#include "test_utilities.h"

igraph_error_t my_gaxpy(const igraph_matrix_t *m,
             const igraph_vector_t *v,
             igraph_vector_t *res) {
    return igraph_blas_dgemv(false, 1.0, m, v, 0.0, res);
}

igraph_bool_t check_same(const igraph_sparsemat_t *A,
                         const igraph_matrix_t *M) {
    igraph_matrix_t A_dense;
    igraph_bool_t result;

    igraph_matrix_init(&A_dense, 1, 1);
    igraph_sparsemat_as_matrix(&A_dense, A);
    result = igraph_matrix_all_e(&A_dense, M);
    igraph_matrix_destroy(&A_dense);

    return result;
}

int main(void) {

    igraph_sparsemat_t A, B, C, D;
    igraph_vector_t v, w, x, y;
    igraph_matrix_t M, N, O;
    igraph_integer_t i;

    RNG_BEGIN();

    /* Matrix-vector product */
#define NROW 10
#define NCOL 5
#define EDGES NROW*NCOL/3
    igraph_matrix_init(&M, NROW, NCOL);
    igraph_sparsemat_init(&A, NROW, NCOL, EDGES);
    for (i = 0; i < EDGES; i++) {
        igraph_integer_t r = RNG_INTEGER(0, NROW - 1);
        igraph_integer_t c = RNG_INTEGER(0, NCOL - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        MATRIX(M, r, c) = MATRIX(M, r, c) + value;
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_vector_init(&v, NCOL);
    igraph_vector_init(&w, NCOL);
    for (i = 0; i < NCOL; i++) {
        VECTOR(v)[i] = VECTOR(w)[i] = RNG_INTEGER(1, 5);
    }

    igraph_vector_init(&x, NROW);
    igraph_vector_init(&y, NROW);
    my_gaxpy(&M, &v, &x);
    igraph_vector_null(&y);
    igraph_sparsemat_gaxpy(&B, &w, &y);

    if (!igraph_vector_all_e(&x, &y)) {
        return 1;
    }

    igraph_vector_destroy(&x);
    igraph_vector_destroy(&y);
    igraph_vector_destroy(&v);
    igraph_vector_destroy(&w);
    igraph_sparsemat_destroy(&B);
    igraph_matrix_destroy(&M);

#undef NROW
#undef NCOL
#undef EDGES

    /* Matrix-matrix product */
#define NROW_A 10
#define NCOL_A 7
#define EDGES_A NROW_A*NCOL_A/3
#define NROW_B 7
#define NCOL_B 9
#define EDGES_B NROW_B*NCOL_B/3
    igraph_matrix_init(&M, NROW_A, NCOL_A);
    igraph_sparsemat_init(&A, NROW_A, NCOL_A, EDGES_A);
    for (i = 0; i < EDGES_A; i++) {
        igraph_integer_t r = RNG_INTEGER(0, NROW_A - 1);
        igraph_integer_t c = RNG_INTEGER(0, NCOL_A - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        MATRIX(M, r, c) = MATRIX(M, r, c) + value;
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &C);
    igraph_sparsemat_destroy(&A);

    igraph_matrix_init(&N, NROW_B, NCOL_B);
    igraph_sparsemat_init(&B, NROW_B, NCOL_B, EDGES_B);
    for (i = 0; i < EDGES_B; i++) {
        igraph_integer_t r = RNG_INTEGER(0, NROW_B - 1);
        igraph_integer_t c = RNG_INTEGER(0, NCOL_B - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        MATRIX(N, r, c) = MATRIX(N, r, c) + value;
        igraph_sparsemat_entry(&B, r, c, value);
    }
    igraph_sparsemat_compress(&B, &D);
    igraph_sparsemat_destroy(&B);

    igraph_matrix_init(&O, 0, 0);
    igraph_blas_dgemm(false, false, 1.0, &M, &N, 0.0, &O);
    igraph_sparsemat_multiply(&C, &D, &A);

    if (! check_same(&A, &O)) {
        return 2;
    }

    igraph_sparsemat_destroy(&C);
    igraph_sparsemat_destroy(&D);
    igraph_sparsemat_destroy(&A);
    igraph_matrix_destroy(&M);
    igraph_matrix_destroy(&N);
    igraph_matrix_destroy(&O);

    VERIFY_FINALLY_STACK();

    RNG_END();

    return 0;
}
