/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "test_utilities.h"

#define N  10
#define M  20
#define NZ 50

#define MIN 0
#define MAX 10

typedef igraph_error_t fun(igraph_sparsemat_t *A, igraph_vector_t *res);

int doit(int which) {

    igraph_integer_t  i;
    igraph_sparsemat_t A, A2;
    igraph_vector_t vec;
    fun *colfun, *rowfun;

    if (which == MIN) {
        colfun = igraph_sparsemat_colmins;
        rowfun = igraph_sparsemat_rowmins;
    } else {
        colfun = igraph_sparsemat_colmaxs;
        rowfun = igraph_sparsemat_rowmaxs;
    }

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Triplet diagonal matrix */

    printf("Triplet diagonal matrix\n");
    igraph_vector_init(&vec, N);
    for (i = 0; i < N; i++) {
        VECTOR(vec)[i] = i;
    }
    igraph_sparsemat_init_diag(&A, /*nzmax=*/ N, /*values=*/ &vec, /*compress=*/ 0);

    igraph_vector_null(&vec);
    rowfun(&A, &vec);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 1;
        }
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 2;
        }
    }

    igraph_vector_destroy(&vec);
    igraph_sparsemat_destroy(&A);

    /* Compressed diagonal matrix */

    printf("Compressed diagonal matrix\n");
    igraph_vector_init(&vec, N);
    for (i = 0; i < N; i++) {
        VECTOR(vec)[i] = i;
    }
    igraph_sparsemat_init_diag(&A, /*nzmax=*/ N, /*values=*/ &vec, /*compress=*/ 1);

    igraph_vector_null(&vec);
    rowfun(&A, &vec);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 3;
        }
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 4;
        }
    }

    igraph_vector_destroy(&vec);
    igraph_sparsemat_destroy(&A);


    /* Random triplet matrix */

    printf("Random triplet matrix\n");
    igraph_sparsemat_init(&A, /*rows=*/ N, /*cols=*/ M, /*nzmax=*/ NZ + 5);
    for (i = 0; i < NZ; i++) {
        int r = igraph_rng_get_integer(igraph_rng_default(), 0, N - 1);
        int c = igraph_rng_get_integer(igraph_rng_default(), 0, M - 1);
        igraph_real_t x = igraph_rng_get_integer(igraph_rng_default(),
                          -10, 10);
        IGRAPH_ASSERT(x >= -10 && x <= 10);
        igraph_sparsemat_entry(&A, r, c, x);
    }
    if (which == MAX) {
        igraph_sparsemat_scale(&A, -1.0);
    }

    igraph_vector_init(&vec, 0);
    colfun(&A, &vec);
    igraph_vector_print(&vec);

    igraph_vector_null(&vec);
    rowfun(&A, &vec);
    igraph_vector_print(&vec);

    /* Random compresssed matrix */

    printf("Random compressed matrix\n");
    igraph_sparsemat_compress(&A, &A2);

    igraph_vector_null(&vec);
    colfun(&A2, &vec);
    igraph_vector_print(&vec);

    igraph_vector_null(&vec);
    rowfun(&A2, &vec);
    igraph_vector_print(&vec);

    igraph_vector_destroy(&vec);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&A2);

    /* Matrix with zero rows, triplet */

    printf("Matrix with zero rows, triplet\n");
    igraph_sparsemat_init(&A, /*rows=*/ 0, /*cols=*/ M, /*nzmax=*/ NZ);
    if (which == MAX) {
        igraph_sparsemat_scale(&A, -1.0);
    }

    igraph_vector_init(&vec, 5);
    rowfun(&A, &vec);
    if (igraph_vector_size(&vec) != 0) {
        return which + 5;
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec);
    igraph_vector_print(&vec);

    /* Matrix with zero rows, compressed */

    printf("Matrix with zero rows, compressed\n");
    igraph_sparsemat_compress(&A, &A2);

    igraph_vector_null(&vec);
    rowfun(&A, &vec);
    if (igraph_vector_size(&vec) != 0) {
        return which + 6;
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec);
    igraph_vector_print(&vec);

    igraph_vector_destroy(&vec);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&A2);

    /* Matrix with zero columns, triplet */

    printf("Matrix with zero columns, triplet\n");
    igraph_sparsemat_init(&A, /*rows=*/ N, /*cols=*/ 0, /*nzmax=*/ NZ);
    if (which == MAX) {
        igraph_sparsemat_scale(&A, -1.0);
    }

    igraph_vector_init(&vec, 5);
    colfun(&A, &vec);
    if (igraph_vector_size(&vec) != 0) {
        return which + 7;
    }

    igraph_vector_null(&vec);
    rowfun(&A, &vec);
    igraph_vector_print(&vec);

    /* Matrix with zero columns, compressed */

    printf("Matrix with zero columns, compressed\n");
    igraph_sparsemat_compress(&A, &A2);

    igraph_vector_null(&vec);
    colfun(&A, &vec);
    if (igraph_vector_size(&vec) != 0) {
        return which + 8;
    }

    igraph_vector_null(&vec);
    rowfun(&A, &vec);
    igraph_vector_print(&vec);

    igraph_vector_destroy(&vec);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&A2);

    return 0;
}

int main(void) {
    int res;

    res = doit(/*which=*/ MIN);
    if (res) {
        return res;
    }

    VERIFY_FINALLY_STACK();

    res = doit(/*which=*/ MAX);
    if (res) {
        return res;
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
