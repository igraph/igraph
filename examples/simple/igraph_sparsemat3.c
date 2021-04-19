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

#include <cs/cs.h>
#include <igraph.h>

int permute(const igraph_matrix_t *M,
            const igraph_vector_int_t *p,
            const igraph_vector_int_t *q,
            igraph_matrix_t *res) {

    long int nrow = igraph_vector_int_size(p);
    long int ncol = igraph_vector_int_size(q);
    long int i, j;

    igraph_matrix_resize(res, nrow, ncol);

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            int ii = VECTOR(*p)[i];
            int jj = VECTOR(*q)[j];
            MATRIX(*res, i, j) = MATRIX(*M, ii, jj);
        }
    }

    return 0;
}

int permute_rows(const igraph_matrix_t *M,
                 const igraph_vector_int_t *p,
                 igraph_matrix_t *res) {

    long int nrow = igraph_vector_int_size(p);
    long int ncol = igraph_matrix_ncol(M);
    long int i, j;

    igraph_matrix_resize(res, nrow, ncol);

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            int ii = VECTOR(*p)[i];
            MATRIX(*res, i, j) = MATRIX(*M, ii, j);
        }
    }

    return 0;
}

int permute_cols(const igraph_matrix_t *M,
                 const igraph_vector_int_t *q,
                 igraph_matrix_t *res) {

    long int nrow = igraph_matrix_nrow(M);
    long int ncol = igraph_vector_int_size(q);
    long int i, j;

    igraph_matrix_resize(res, nrow, ncol);

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            int jj = VECTOR(*q)[j];
            MATRIX(*res, i, j) = MATRIX(*M, i, jj);
        }
    }

    return 0;
}

int random_permutation(igraph_vector_int_t *vec) {
    /* We just do size(vec) * 2 swaps */
    long int one, two, i, n = igraph_vector_int_size(vec);
    int tmp;
    for (i = 0; i < 2 * n; i++) {
        one = RNG_INTEGER(0, n - 1);
        two = RNG_INTEGER(0, n - 1);
        tmp = VECTOR(*vec)[one];
        VECTOR(*vec)[one] = VECTOR(*vec)[two];
        VECTOR(*vec)[two] = tmp;
    }
    return 0;
}

igraph_bool_t check_same(const igraph_sparsemat_t *A,
                         const igraph_matrix_t *M) {

    long int nrow = igraph_sparsemat_nrow(A);
    long int ncol = igraph_sparsemat_ncol(A);
    long int j, p, nzero = 0;

    if (nrow != igraph_matrix_nrow(M) ||
        ncol != igraph_matrix_ncol(M)) {
        return 0;
    }

    for (j = 0; j < A->cs->n; j++) {
        for (p = A->cs->p[j]; p < A->cs->p[j + 1]; p++) {
            long int to = A->cs->i[p];
            igraph_real_t value = A->cs->x[p];
            if (value != MATRIX(*M, to, j)) {
                return 0;
            }
            nzero += 1;
        }
    }

    for (j = 0; j < nrow; j++) {
        for (p = 0; p < ncol; p++) {
            if (MATRIX(*M, j, p) != 0) {
                nzero -= 1;
            }
        }
    }

    return nzero == 0;
}

int main() {

    igraph_sparsemat_t A, B;
    igraph_matrix_t M, N;
    igraph_vector_int_t p, q;
    long int i;

    /* Permutation of a matrix */

#define NROW 10
#define NCOL 5
#define EDGES NROW*NCOL/3
    igraph_matrix_init(&M, NROW, NCOL);
    igraph_sparsemat_init(&A, NROW, NCOL, EDGES);
    for (i = 0; i < EDGES; i++) {
        long int r = RNG_INTEGER(0, NROW - 1);
        long int c = RNG_INTEGER(0, NCOL - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        MATRIX(M, r, c) = MATRIX(M, r, c) + value;
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_vector_int_init_seq(&p, 0, NROW - 1);
    igraph_vector_int_init_seq(&q, 0, NCOL - 1);

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
        long int r = RNG_INTEGER(0, NROW - 1);
        long int c = RNG_INTEGER(0, NCOL - 1);
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

    /* A single element */

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

    igraph_sparsemat_destroy(&A);

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
        long int r = RNG_INTEGER(0, NROW - 1);
        long int c = RNG_INTEGER(0, NCOL - 1);
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

    return 0;
}
