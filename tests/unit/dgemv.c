/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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

#include <igraph.h>

#include "test_utilities.h"

/* Matrix-vector multiplication: y = A.x */
void matmul(const igraph_matrix_t *A, const igraph_vector_t *x, igraph_vector_t *y, igraph_real_t beta) {
    igraph_integer_t i, j, nr = igraph_matrix_nrow(A), nc = igraph_matrix_ncol(A);

    IGRAPH_ASSERT(nc == igraph_vector_size(x));
    IGRAPH_ASSERT(nr == igraph_vector_size(y));

    igraph_vector_scale(y, beta);

    for (i=0; i < nr; ++i) {
        for (j=0; j < nc; ++j) {
            VECTOR(*y)[i] += MATRIX(*A, i, j) * VECTOR(*x)[j];
        }
    }
}

int main(void) {
    igraph_matrix_t A;
    igraph_vector_t x, y1, y2;
    igraph_integer_t i, j;
    const igraph_integer_t nr = 5, nc = 8;

    igraph_rng_seed(igraph_rng_default(), 54632);

    igraph_matrix_init(&A, nr, nc);
    igraph_vector_init(&x, nc);

    /* Fill with arbitrary values. Should be zeroes by beta. */
    igraph_vector_init_range(&y1, 1, nr + 1);
    igraph_vector_init_copy(&y2, &y1);

    for (i=0; i < nr; ++i) {
        for (j=0; j < nc; ++j) {
            MATRIX(A, i, j) = (igraph_real_t) RNG_INTEGER(-10, 10);
        }
    }

    for (j=0; j < nc; ++j) {
        VECTOR(x)[j] = (igraph_real_t) RNG_INTEGER(-10, 10);
    }

    printf("Input matrix A:\n");
    print_matrix(&A);

    printf("\nInput vector x:\n");
    print_vector(&x);

    igraph_blas_dgemv(0, 1, &A, &x, 0, &y1);
    matmul(&A, &x, &y2, 0);

    printf("\nResult vector DGEMV:\n");
    print_vector(&y1);

    printf("\nResult vector naive:\n");
    print_vector(&y2);

    /* Results should be exact since all values are integers */
    IGRAPH_ASSERT(igraph_vector_all_e(&y1, &y2));

    printf("\nAdding to previous result with beta=2:\n");

    igraph_blas_dgemv(0, 1, &A, &x, 2, &y1);
    matmul(&A, &x, &y2, 2);

    printf("\nResult vector DGEMV:\n");
    print_vector(&y1);

    printf("\nResult vector naive:\n");
    print_vector(&y2);

    IGRAPH_ASSERT(igraph_vector_all_e(&y1, &y2));

    igraph_vector_destroy(&y2);
    igraph_vector_destroy(&y1);
    igraph_vector_destroy(&x);
    igraph_matrix_destroy(&A);

    VERIFY_FINALLY_STACK();

    return 0;
}
