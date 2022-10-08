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

#include <igraph.h>

void permute(const igraph_matrix_t *M,
            const igraph_vector_int_t *p,
            const igraph_vector_int_t *q,
            igraph_matrix_t *res) {

    igraph_integer_t nrow = igraph_vector_int_size(p);
    igraph_integer_t ncol = igraph_vector_int_size(q);
    igraph_integer_t i, j;

    igraph_matrix_resize(res, nrow, ncol);

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            igraph_integer_t ii = VECTOR(*p)[i];
            igraph_integer_t jj = VECTOR(*q)[j];
            MATRIX(*res, i, j) = MATRIX(*M, ii, jj);
        }
    }
}

void permute_rows(const igraph_matrix_t *M,
                 const igraph_vector_int_t *p,
                 igraph_matrix_t *res) {

    igraph_integer_t nrow = igraph_vector_int_size(p);
    igraph_integer_t ncol = igraph_matrix_ncol(M);
    igraph_integer_t i, j;

    igraph_matrix_resize(res, nrow, ncol);

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            igraph_integer_t ii = VECTOR(*p)[i];
            MATRIX(*res, i, j) = MATRIX(*M, ii, j);
        }
    }
}

void permute_cols(const igraph_matrix_t *M,
                 const igraph_vector_int_t *q,
                 igraph_matrix_t *res) {

    igraph_integer_t nrow = igraph_matrix_nrow(M);
    igraph_integer_t ncol = igraph_vector_int_size(q);
    igraph_integer_t i, j;

    igraph_matrix_resize(res, nrow, ncol);

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            igraph_integer_t jj = VECTOR(*q)[j];
            MATRIX(*res, i, j) = MATRIX(*M, i, jj);
        }
    }
}

void random_permutation(igraph_vector_int_t *vec) {
    /* We just do size(vec) * 2 swaps */
    igraph_integer_t one, two, i, n = igraph_vector_int_size(vec);
    igraph_integer_t tmp;
    for (i = 0; i < 2 * n; i++) {
        one = RNG_INTEGER(0, n - 1);
        two = RNG_INTEGER(0, n - 1);
        tmp = VECTOR(*vec)[one];
        VECTOR(*vec)[one] = VECTOR(*vec)[two];
        VECTOR(*vec)[two] = tmp;
    }
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

    igraph_sparsemat_t A, B;
    igraph_matrix_t M, N;
    igraph_vector_int_t p, q;
    igraph_integer_t i;

    RNG_BEGIN();

    /* Permutation of a matrix */

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

    igraph_vector_int_init_range(&p, 0, NROW);
    igraph_vector_int_init_range(&q, 0, NCOL);

    /* Identity */

    igraph_matrix_init(&N, 0, 0);
    permute(&M, &p, &q, &N);

    igraph_sparsemat_permute(&B, &p, &q, &A);
    igraph_sparsemat_dupl(&A);

    if (! check_same(&A, &N)) {
        return 1;
    }

    /* Random permutation */
    random_permutation(&p);
    random_permutation(&q);

    permute(&M, &p, &q, &N);

    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_permute(&B, &p, &q, &A);
    igraph_sparsemat_dupl(&A);

    if (! check_same(&A, &N)) {
        return 2;
    }

    igraph_vector_int_destroy(&p);
    igraph_vector_int_destroy(&q);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);
    igraph_matrix_destroy(&M);
    igraph_matrix_destroy(&N);

#undef NROW
#undef NCOL
#undef EDGES

    /* Indexing */

#define NROW 10
#define NCOL 5
#define EDGES NROW*NCOL/3
#define I_NROW 6
#define I_NCOL 3
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

    igraph_vector_int_init(&p, I_NROW);
    igraph_vector_int_init(&q, I_NCOL);

    for (i = 0; i < I_NROW; i++) {
        VECTOR(p)[i] = RNG_INTEGER(0, I_NROW - 1);
    }
    for (i = 0; i < I_NCOL; i++) {
        VECTOR(p)[i] = RNG_INTEGER(0, I_NCOL - 1);
    }

    igraph_matrix_init(&N, 0, 0);
    permute(&M, &p, &q, &N);

    igraph_sparsemat_index(&B, &p, &q, &A, 0);

    if (! check_same(&A, &N)) {
        return 3;
    }

    igraph_sparsemat_destroy(&A);

    /* Getting single elements with index() */

    igraph_vector_int_resize(&p, 1);
    igraph_vector_int_resize(&q, 1);

    for (i = 0; i < 100; i++) {
        igraph_real_t value;
        VECTOR(p)[0] = RNG_INTEGER(0, NROW - 1);
        VECTOR(q)[0] = RNG_INTEGER(0, NCOL - 1);
        igraph_sparsemat_index(&B, &p, &q, /*res=*/ 0, &value);
        if (value != MATRIX(M, VECTOR(p)[0], VECTOR(q)[0])) {
            return 4;
        }
    }

    /* Getting single elements with get() */

    igraph_vector_int_resize(&p, 1);
    igraph_vector_int_resize(&q, 1);

    for (i = 0; i < 100; i++) {
        igraph_integer_t row = RNG_INTEGER(0, NROW - 1);
        igraph_integer_t col = RNG_INTEGER(0, NCOL - 1);
        if (igraph_sparsemat_get(&B, row, col) != MATRIX(M, row, col)) {
            return 4;
        }
    }

    /* Getting submatrices with index() */

    for (i = 0; i < 100; i++) {
        igraph_real_t value;
        VECTOR(p)[0] = RNG_INTEGER(0, NROW - 1);
        VECTOR(q)[0] = RNG_INTEGER(0, NCOL - 1);
        igraph_sparsemat_index(&B, &p, &q, /*res=*/ &A, &value);
        igraph_sparsemat_destroy(&A);
        if (value != MATRIX(M, VECTOR(p)[0], VECTOR(q)[0])) {
            return 4;
        }
    }

    igraph_vector_int_destroy(&p);
    igraph_vector_int_destroy(&q);
    igraph_sparsemat_destroy(&B);
    igraph_matrix_destroy(&M);
    igraph_matrix_destroy(&N);

    /* Indexing only the rows or the columns */

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

    igraph_vector_int_init(&p, I_NROW);
    igraph_vector_int_init(&q, I_NCOL);

    for (i = 0; i < I_NROW; i++) {
        VECTOR(p)[i] = RNG_INTEGER(0, I_NROW - 1);
    }
    for (i = 0; i < I_NCOL; i++) {
        VECTOR(p)[i] = RNG_INTEGER(0, I_NCOL - 1);
    }

    igraph_matrix_init(&N, 0, 0);
    permute_rows(&M, &p, &N);

    igraph_sparsemat_index(&B, &p, 0, &A, 0);

    if (! check_same(&A, &N)) {
        return 5;
    }

    permute_cols(&M, &q, &N);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_index(&B, 0, &q, &A, 0);

    if (! check_same(&A, &N)) {
        return 6;
    }

    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);
    igraph_vector_int_destroy(&p);
    igraph_vector_int_destroy(&q);
    igraph_matrix_destroy(&M);
    igraph_matrix_destroy(&N);

    RNG_END();

    return 0;
}
