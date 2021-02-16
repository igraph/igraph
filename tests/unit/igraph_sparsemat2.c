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
#include "linalg/blas_internal.h"
#include "linalg/arpack_internal.h"

#include "test_utilities.inc"

int igraph_matrix_dgemv(const igraph_matrix_t *m,
                        const igraph_vector_t *v,
                        igraph_vector_t *res,
                        igraph_real_t alpha,
                        igraph_real_t beta,
                        igraph_bool_t transpose_m) {

    int nrow = igraph_matrix_nrow(m);
    int ncol = igraph_matrix_ncol(m);
    long int vlen = igraph_vector_size(v);
    int one = 1;
    char t = transpose_m ? 't' : 'n';
    long int input_len  = transpose_m ? nrow : ncol;
    long int output_len = transpose_m ? ncol : nrow;

    if (vlen != input_len) {
        IGRAPH_ERROR("Matrix and vector sizes are incompatible", IGRAPH_EINVAL);
    }

    if (beta != 0 && igraph_vector_size(res) != output_len) {
        IGRAPH_ERROR("Non-zero beta and bad `res' vector size, possible mistake",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_resize(res, output_len));

    igraphdgemv_(&t, &nrow, &ncol, &alpha, &MATRIX(*m, 0, 0),
                 &nrow, VECTOR(*v), &one, &beta, VECTOR(*res), &one);

    return 0;
}

int igraph_matrix_vector_prod(const igraph_matrix_t *m,
                              const igraph_vector_t *v,
                              igraph_vector_t *res) {
    return igraph_matrix_dgemv(m, v, res, 1.0, 0.0, /*transpose=*/ 0);
}

int my_dgemv(const igraph_matrix_t *m,
             const igraph_vector_t *v,
             igraph_vector_t *res,
             igraph_real_t alpha,
             igraph_real_t beta,
             igraph_bool_t transpose_m) {

    int nrow = igraph_matrix_nrow(m);
    int ncol = igraph_matrix_ncol(m);
    long int vlen = igraph_vector_size(v);
    int one = 1;
    char t = transpose_m ? 't' : 'n';
    long int input_len  = transpose_m ? nrow : ncol;
    long int output_len = transpose_m ? ncol : nrow;

    if (vlen != input_len) {
        IGRAPH_ERROR("Matrix and vector sizes are incompatible", IGRAPH_EINVAL);
    }

    if (beta != 0 && igraph_vector_size(res) != output_len) {
        IGRAPH_ERROR("Non-zero beta and bad `res' vector size, possible mistake",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_resize(res, output_len));

    igraphdgemv_(&t, &nrow, &ncol, &alpha, &MATRIX(*m, 0, 0),
                 &nrow, VECTOR(*v), &one, &beta, VECTOR(*res), &one);

    return 0;
}

int my_gaxpy(const igraph_matrix_t *m,
             const igraph_vector_t *v,
             igraph_vector_t *res) {
    return my_dgemv(m, v, res, 1.0, 0.0, /*transpose=*/ 0);
}

int my_dgemm(const igraph_matrix_t *m1,
             const igraph_matrix_t *m2,
             igraph_matrix_t *res) {

    long int m1_r = igraph_matrix_nrow(m1);
    long int m1_c = igraph_matrix_ncol(m1);
    long int m2_r = igraph_matrix_nrow(m2);
    long int m2_c = igraph_matrix_ncol(m2);
    long int i, j, k;

    if (m1_c != m2_r) {
        IGRAPH_ERROR("Cannot multiply matrices, invalid dimensions", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, m1_r, m2_c));
    igraph_matrix_null(res);

    for (i = 0; i < m1_r; i++) {
        for (j = 0; j < m2_c; j++) {
            for (k = 0; k < m1_c /* which is also m2_r*/; k++) {
                MATRIX(*res, i, j) += MATRIX(*m1, i, k) * MATRIX(*m2, k, j);
            }
        }
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

    igraph_sparsemat_t A, B, C, D;
    igraph_vector_t v, w, x, y;
    igraph_matrix_t M, N, O;
    long int i;

    srand(1);

    /* Matrix-vector product */
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
        long int r = RNG_INTEGER(0, NROW_A - 1);
        long int c = RNG_INTEGER(0, NCOL_A - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        MATRIX(M, r, c) = MATRIX(M, r, c) + value;
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &C);
    igraph_sparsemat_destroy(&A);

    igraph_matrix_init(&N, NROW_B, NCOL_B);
    igraph_sparsemat_init(&B, NROW_B, NCOL_B, EDGES_B);
    for (i = 0; i < EDGES_B; i++) {
        long int r = RNG_INTEGER(0, NROW_B - 1);
        long int c = RNG_INTEGER(0, NCOL_B - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        MATRIX(N, r, c) = MATRIX(N, r, c) + value;
        igraph_sparsemat_entry(&B, r, c, value);
    }
    igraph_sparsemat_compress(&B, &D);
    igraph_sparsemat_destroy(&B);

    igraph_matrix_init(&O, 0, 0);
    my_dgemm(&M, &N, &O);
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

    return 0;
}
