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

#include "test_utilities.inc"

#define N  10
#define M  20
#define NZ 50

#define MIN 0
#define MAX 10

typedef int fun(igraph_sparsemat_t *A, igraph_vector_t *res,
                igraph_vector_int_t *pos);

int doit(int which) {

    int i;
    igraph_sparsemat_t A, A2;
    igraph_vector_t vec;
    igraph_vector_int_t pos;
    fun *colfun, *rowfun;

    if (which == MIN) {
        colfun = igraph_sparsemat_which_min_cols;
        rowfun = igraph_sparsemat_which_min_rows;
    } else {
        /* colfun = */ /* TODO */
        /* rowfun = */ /* TODO */
    }

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Triplet diagonal matrix */

    igraph_vector_init(&vec, N);
    igraph_vector_int_init(&pos, N);
    for (i = 0; i < N; i++) {
        VECTOR(vec)[i] = i;
    }
    igraph_sparsemat_diag(&A, /*nzmax=*/ N, /*values=*/ &vec,
                          /*compress=*/ 0);

    igraph_vector_null(&vec);
    igraph_vector_int_null(&pos);
    rowfun(&A, &vec, &pos);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 1;
        }
    }
    for (i = 0; i < N; i++) {
        if (VECTOR(pos)[i] != i) {
            return which + 2;
        }
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec, &pos);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 3;
        }
    }
    for (i = 0; i < N; i++) {
        if (VECTOR(pos)[i] != i) {
            return which + 4;
        }
    }

    igraph_vector_destroy(&vec);
    igraph_vector_int_destroy(&pos);
    igraph_sparsemat_destroy(&A);

    /* Compressed diagonal matrix */

    igraph_vector_init(&vec, N);
    igraph_vector_int_init(&pos, N);
    for (i = 0; i < N; i++) {
        VECTOR(vec)[i] = i;
    }
    igraph_sparsemat_diag(&A, /*nzmax=*/ N, /*values=*/ &vec,
                          /*compress=*/ 1);

    igraph_vector_null(&vec);
    rowfun(&A, &vec, &pos);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 5;
        }
    }
    for (i = 0; i < N; i++) {
        if (VECTOR(pos)[i] != i) {
            return which + 6;
        }
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec, &pos);
    for (i = 0; i < N; i++) {
        if (VECTOR(vec)[i] != i) {
            return which + 7;
        }
    }
    for (i = 0; i < N; i++) {
        if (VECTOR(pos)[i] != i) {
            return which + 8;
        }
    }

    igraph_vector_destroy(&vec);
    igraph_vector_int_destroy(&pos);
    igraph_sparsemat_destroy(&A);


    /* Random triplet matrix */

    igraph_sparsemat_init(&A, /*rows=*/ N, /*cols=*/ M, /*nzmax=*/ NZ + 5);
    for (i = 0; i < NZ; i++) {
        int r = igraph_rng_get_integer(igraph_rng_default(), 0, N - 1);
        int c = igraph_rng_get_integer(igraph_rng_default(), 0, M - 1);
        igraph_real_t x = igraph_rng_get_integer(igraph_rng_default(),
                          -10, 10);
        igraph_sparsemat_entry(&A, r, c, x);
    }
    if (which == MAX) {
        igraph_sparsemat_scale(&A, -1.0);
    }

    igraph_vector_init(&vec, 0);
    igraph_vector_int_init(&pos, 0);
    colfun(&A, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    igraph_vector_null(&vec);
    rowfun(&A, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    /* Random compresssed matrix */

    igraph_sparsemat_compress(&A, &A2);

    igraph_vector_null(&vec);
    colfun(&A2, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    igraph_vector_null(&vec);
    rowfun(&A2, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    igraph_vector_destroy(&vec);
    igraph_vector_int_destroy(&pos);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&A2);

    /* Matrix with zero rows, triplet */

    igraph_sparsemat_init(&A, /*rows=*/ 0, /*cols=*/ M, /*nzmax=*/ NZ);
    if (which == MAX) {
        igraph_sparsemat_scale(&A, -1.0);
    }

    igraph_vector_init(&vec, 5);
    igraph_vector_int_init(&pos, 5);
    rowfun(&A, &vec, &pos);
    if (igraph_vector_size(&vec) != 0) {
        return which + 5;
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    /* Matrix with zero rows, compressed */

    igraph_sparsemat_compress(&A, &A2);

    igraph_vector_null(&vec);
    rowfun(&A, &vec, &pos);
    if (igraph_vector_size(&vec) != 0) {
        return which + 6;
    }

    igraph_vector_null(&vec);
    colfun(&A, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    igraph_vector_destroy(&vec);
    igraph_vector_int_destroy(&pos);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&A2);

    /* Matrix with zero columns, triplet */

    igraph_sparsemat_init(&A, /*rows=*/ N, /*cols=*/ 0, /*nzmax=*/ NZ);
    if (which == MAX) {
        igraph_sparsemat_scale(&A, -1.0);
    }

    igraph_vector_init(&vec, 5);
    igraph_vector_int_init(&pos, 5);
    colfun(&A, &vec, &pos);
    if (igraph_vector_size(&vec) != 0) {
        return which + 7;
    }

    igraph_vector_null(&vec);
    rowfun(&A, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    /* Matrix with zero columns, compressed */

    igraph_sparsemat_compress(&A, &A2);

    igraph_vector_null(&vec);
    colfun(&A, &vec, &pos);
    if (igraph_vector_size(&vec) != 0) {
        return which + 8;
    }

    igraph_vector_null(&vec);
    rowfun(&A, &vec, &pos);
    igraph_vector_print(&vec);
    igraph_vector_int_print(&pos);

    igraph_vector_destroy(&vec);
    igraph_vector_int_destroy(&pos);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&A2);

    return 0;
}

int main() {
    int res;

    res = doit(/*which=*/ MIN);
    if (res) {
        return res;
    }

    /* res = doit(/\*which=*\/ MAX); */
    /* if (res) { return res; } */

    VERIFY_FINALLY_STACK();

    return 0;
}
