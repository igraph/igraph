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

int real_cplx_mult(const igraph_matrix_t *A,
                   const igraph_vector_t *v_real,
                   const igraph_vector_t *v_imag,
                   igraph_vector_t *res_real,
                   igraph_vector_t *res_imag) {

    int n = igraph_vector_size(v_real);
    int r, c;

    if (igraph_matrix_nrow(A) != n ||
        igraph_matrix_ncol(A) != n ||
        igraph_vector_size(v_imag) != n) {
        printf("Wrong matrix or vector size");
        return 1;
    }

    igraph_vector_resize(res_real, n);
    igraph_vector_resize(res_imag, n);

    for (r = 0; r < n; r++) {
        igraph_real_t s_real = 0.0;
        igraph_real_t s_imag = 0.0;
        for (c = 0; c < n; c++) {
            s_real += MATRIX(*A, r, c) * VECTOR(*v_real)[c];
            s_imag += MATRIX(*A, r, c) * VECTOR(*v_imag)[c];
        }
        VECTOR(*res_real)[r] = s_real;
        VECTOR(*res_imag)[r] = s_imag;
    }

    return 0;
}

int sc_cplx_cplx_mult(igraph_real_t lambda_real,
                      igraph_real_t lambda_imag,
                      const igraph_vector_t *v_real,
                      const igraph_vector_t *v_imag,
                      igraph_vector_t *res_real,
                      igraph_vector_t *res_imag) {

    int r;
    int n = igraph_vector_size(v_real);

    if (igraph_vector_size(v_imag) != n) {
        printf("Wrong vector sizes");
        return 1;
    }

    igraph_vector_resize(res_real, n);
    igraph_vector_resize(res_imag, n);

    for (r = 0; r < n; r++) {
        VECTOR(*res_real)[r] = (lambda_real * VECTOR(*v_real)[r] -
                                lambda_imag * VECTOR(*v_imag)[r]);
        VECTOR(*res_imag)[r] = (lambda_imag * VECTOR(*v_real)[r] +
                                lambda_real * VECTOR(*v_imag)[r]);
    }

    return 0;
}

igraph_bool_t check_ev(const igraph_matrix_t *A,
                       const igraph_vector_t *values_real,
                       const igraph_vector_t *values_imag,
                       const igraph_matrix_t *vectors_left,
                       const igraph_matrix_t *vectors_right,
                       igraph_real_t tol) {

    int i, n = igraph_matrix_nrow(A);
    igraph_vector_t v_real, v_imag;
    igraph_vector_t AV_real, AV_imag, lv_real, lv_imag;
    igraph_vector_t null;

    if (igraph_matrix_ncol(A)             != n) {
        return 1;
    }
    if (igraph_vector_size(values_real)   != n) {
        return 1;
    }
    if (igraph_vector_size(values_imag)   != n) {
        return 1;
    }
    if (igraph_matrix_nrow(vectors_left)  != n) {
        return 1;
    }
    if (igraph_matrix_ncol(vectors_left)  != n) {
        return 1;
    }
    if (igraph_matrix_nrow(vectors_right) != n) {
        return 1;
    }
    if (igraph_matrix_ncol(vectors_right) != n) {
        return 1;
    }

    igraph_vector_init(&AV_real, n);
    igraph_vector_init(&AV_imag, n);
    igraph_vector_init(&lv_real, n);
    igraph_vector_init(&lv_imag, n);
    igraph_vector_init(&null, n);
    igraph_vector_null(&null);

    for (i = 0; i < n; i++) {
        if (VECTOR(*values_imag)[i] == 0.0) {
            igraph_vector_view(&v_real, &MATRIX(*vectors_right, 0, i), n);
            igraph_vector_view(&v_imag, VECTOR(null), n);
        } else if (VECTOR(*values_imag)[i] > 0.0) {
            igraph_vector_view(&v_real, &MATRIX(*vectors_right, 0, i), n);
            igraph_vector_view(&v_imag, &MATRIX(*vectors_right, 0, i + 1), n);
        } else if (VECTOR(*values_imag)[i] < 0.0) {
            igraph_vector_view(&v_real, &MATRIX(*vectors_right, 0, i - 1), n);
            igraph_vector_view(&v_imag, &MATRIX(*vectors_right, 0, i), n);
            igraph_vector_scale(&v_imag, -1.0);
        }
        real_cplx_mult(A, &v_real, &v_imag, &AV_real, &AV_imag);
        sc_cplx_cplx_mult(VECTOR(*values_real)[i], VECTOR(*values_imag)[i],
                          &v_real, &v_imag, &lv_real, &lv_imag);

        if (igraph_vector_maxdifference(&AV_real, &lv_real) > tol ||
            igraph_vector_maxdifference(&AV_imag, &lv_imag) > tol) {
            printf("ERROR:\n");
            igraph_vector_print(&AV_real);
            igraph_vector_print(&AV_imag);
            igraph_vector_print(&lv_real);
            igraph_vector_print(&lv_imag);
            return 1;
        }
    }

    igraph_vector_destroy(&null);
    igraph_vector_destroy(&AV_imag);
    igraph_vector_destroy(&AV_real);
    igraph_vector_destroy(&lv_imag);
    igraph_vector_destroy(&lv_real);

    return 0;
}

int main() {

    igraph_matrix_t A;
    igraph_matrix_t vectors_left, vectors_right;
    igraph_vector_t values_real, values_imag;
    int i, j;
    int info = 1;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_matrix_init(&A, DIM, DIM);
    igraph_matrix_init(&vectors_left, 0, 0);
    igraph_matrix_init(&vectors_right, 0, 0);
    igraph_vector_init(&values_real, 0);
    igraph_vector_init(&values_imag, 0);

    for (i = 0; i < DIM; i++) {
        for (j = 0; j < DIM; j++) {
            MATRIX(A, i, j) = igraph_rng_get_integer(igraph_rng_default(), 1, 10);
        }
    }

    igraph_lapack_dgeev(&A, &values_real, &values_imag,
                        &vectors_left, &vectors_right, &info);

    if (check_ev(&A, &values_real, &values_imag,
                 &vectors_left, &vectors_right, /*tol=*/ 1e-8)) {
        return 1;
    }

    /* ------------------------------------------------------- */

    /* igraph_matrix_resize(&A, 10, 10); */
    /* igraph_matrix_null(&A); */
    /* for (i=0; i<10; i++) { MATRIX(A, i, i) = 1.0; }  */
    /* MATRIX(A,0,1) = 1.0; */

    /* igraph_lapack_dgeev(&A, &values_real, &values_imag,  */
    /*              &vectors_left, &vectors_right, &info);   */

    /* if (check_ev(&A, &values_real, &values_imag, */
    /*           &vectors_left, &vectors_right, /\*tol=*\/ 1e-8)) { */
    /*   return 2; */
    /* } */

    /* ------------------------------------------------------- */

    igraph_matrix_resize(&A, 10, 10);
    igraph_matrix_null(&A);
    MATRIX(A, 0, 1) = MATRIX(A, 0, 2) = MATRIX(A, 0, 3) = 1 / 3.0;
    MATRIX(A, 1, 0) = MATRIX(A, 1, 4) = MATRIX(A, 1, 5) = MATRIX(A, 1, 6) = 1 / 4.0;
    MATRIX(A, 2, 0) = MATRIX(A, 2, 7) = MATRIX(A, 2, 8) = MATRIX(A, 2, 9) = 1 / 4.0;
    MATRIX(A, 3, 0) = 1.0;
    MATRIX(A, 4, 1) = 1.0;
    MATRIX(A, 5, 1) = 1.0;
    MATRIX(A, 6, 1) = 1.0;
    MATRIX(A, 7, 2) = 1.0;
    MATRIX(A, 8, 2) = 1.0;
    MATRIX(A, 9, 2) = 1.0;

    info = 0;
    igraph_lapack_dgeev(&A, &values_real, &values_imag,
                        &vectors_left, &vectors_right, &info);

    /* igraph_matrix_print(&A); */
    /* printf("---\n"); */
    /* igraph_vector_print(&values_real); */
    /* igraph_vector_print(&values_imag); */
    /* igraph_matrix_print(&vectors_left); */

    if (check_ev(&A, &values_real, &values_imag,
                 &vectors_left, &vectors_right, /*tol=*/ 1e-8)) {
        return 3;
    }

    igraph_vector_destroy(&values_imag);
    igraph_vector_destroy(&values_real);
    igraph_matrix_destroy(&vectors_right);
    igraph_matrix_destroy(&vectors_left);
    igraph_matrix_destroy(&A);

    return 0;
}
