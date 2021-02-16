/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

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

#define DIM1 10
#define DIM2 5
#define DIM3 6

#define INT(a) (igraph_rng_get_integer(igraph_rng_default(), 0, (a)))
#define REAL() (igraph_rng_get_normal(igraph_rng_default(), 0, 1))

int main() {
    igraph_sparsemat_t sA, sB, sC;
    igraph_matrix_t A1, A2, A3, B, C;
    int i;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_sparsemat_init(&sA, DIM1, DIM2, 20);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&sA, INT(DIM1 - 1), INT(DIM2 - 1), REAL());
    }
    igraph_sparsemat_compress(&sA, &sB);
    igraph_sparsemat_destroy(&sA);

    igraph_sparsemat_init(&sA, DIM2, DIM3, 20);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&sA, INT(DIM2 - 1), INT(DIM3 - 1), REAL());
    }
    igraph_sparsemat_compress(&sA, &sC);
    igraph_sparsemat_destroy(&sA);

    igraph_matrix_init(&B, 0, 0);
    igraph_sparsemat_as_matrix(&B, &sB);
    igraph_matrix_init(&C, 0, 0);
    igraph_sparsemat_as_matrix(&C, &sC);

    /* All possible products */
    igraph_sparsemat_multiply(&sB, &sC, &sA);
    igraph_matrix_init(&A1, 0, 0);
    igraph_sparsemat_as_matrix(&A1, &sA);
    igraph_matrix_init(&A2, 0, 0);
    igraph_sparsemat_dense_multiply(&B, &sC, &A2);
    igraph_matrix_init(&A3, 0, 0);
    igraph_sparsemat_multiply_by_dense(&sB, &C, &A3);

    if (igraph_matrix_maxdifference(&A1, &A2) > 1e-10 ||
        igraph_matrix_maxdifference(&A2, &A3) > 1e-10) {
        return 1;
    }

    igraph_sparsemat_destroy(&sA);
    igraph_sparsemat_destroy(&sB);
    igraph_sparsemat_destroy(&sC);

    igraph_matrix_destroy(&A1);
    igraph_matrix_destroy(&A2);
    igraph_matrix_destroy(&A3);
    igraph_matrix_destroy(&B);
    igraph_matrix_destroy(&C);

    VERIFY_FINALLY_STACK();

    return 0;
}
