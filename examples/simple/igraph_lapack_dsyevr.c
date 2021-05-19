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

#define DIM 10

igraph_bool_t check_ev(const igraph_matrix_t *A,
                       const igraph_vector_t *values,
                       const igraph_matrix_t *vectors, igraph_real_t tol) {
    igraph_vector_t v, y;
    int i, j;
    int m = igraph_matrix_ncol(vectors);
    int n = igraph_matrix_nrow(A);

    if (igraph_matrix_ncol(A) != n)       {
        return 1;
    }
    if (igraph_vector_size(values) != m)  {
        return 1;
    }
    if (igraph_matrix_nrow(vectors) != n) {
        return 1;
    }

    igraph_vector_init(&y, n);

    for (i = 0; i < m; i++) {
        igraph_vector_view(&v, &MATRIX(*vectors, 0, i), n);
        igraph_vector_update(&y, &v);
        igraph_blas_dgemv(/*transpose=*/ 0, /*alpha=*/ 1.0, A, &v,
                                         /*beta=*/ -VECTOR(*values)[i], &y);
        for (j = 0; j < n; j++) {
            if (fabs(VECTOR(y)[i]) > tol) {
                printf("Matrix:\n");
                igraph_matrix_print(A);
                printf("lambda= %g\n", VECTOR(*values)[i]);
                printf("v= ");
                igraph_vector_print(&v);
                printf("residual: ");
                igraph_vector_print(&y);
                return 1;
            }
        }
    }

    igraph_vector_destroy(&y);
    return 0;
}

int main() {

    igraph_matrix_t A;
    igraph_matrix_t vectors, vectors2;
    igraph_vector_t values, values2;
    int i, j;
    int il, iu;
    igraph_real_t vl, vu;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_matrix_init(&A, DIM, DIM);
    igraph_matrix_init(&vectors, 0, 0);
    igraph_vector_init(&values, 0);

    /* All eigenvalues and eigenvectors */

    for (i = 0; i < DIM; i++) {
        for (j = i; j < DIM; j++) {
            MATRIX(A, i, j) = MATRIX(A, j, i) =
                                  igraph_rng_get_integer(igraph_rng_default(), 1, 10);
        }
    }

    igraph_lapack_dsyevr(&A, IGRAPH_LAPACK_DSYEV_ALL, /*vl=*/ 0, /*vu=*/ 0,
                         /*vestimate=*/ 0, /*il=*/ 0, /*iu=*/ 0,
                         /*abstol=*/ 1e-10, &values, &vectors, /*support=*/ 0);

    if (igraph_vector_size(&values) != DIM) {
        return 1;
    }
    if (igraph_matrix_nrow(&vectors) != DIM ||
        igraph_matrix_ncol(&vectors) != DIM) {
        return 2;
    }
    if (check_ev(&A, &values, &vectors, /*tol=*/ 1e-8)) {
        return 3;
    }

    /* Only a subset */

    igraph_matrix_init(&vectors2, 0, 0);
    igraph_vector_init(&values2, 0);

    il = 2;
    iu = 5;
    igraph_lapack_dsyevr(&A, IGRAPH_LAPACK_DSYEV_SELECT, /*vl=*/ 0, /*vu=*/ 0,
                         /*vestimate=*/ 0, /*il=*/ il, /*iu=*/ iu,
                         /*abstol=*/ 1e-10, &values2, &vectors2,
                         /*support=*/ 0);

    if (igraph_vector_size(&values2) != iu - il + 1) {
        return 4;
    }
    if (igraph_matrix_nrow(&vectors2) != DIM ||
        igraph_matrix_ncol(&vectors2) != iu - il + 1) {
        return 5;
    }
    for (i = 0; i < iu - il + 1; i++) {
        igraph_real_t m1 = 1.0;

        if (fabs(VECTOR(values)[il + i - 1] - VECTOR(values2)[i]) > 1e-8) {
            printf("Full:   ");
            igraph_vector_print(&values);
            printf("Subset: ");
            igraph_vector_print(&values2);
            return 6;
        }

        if (MATRIX(vectors, 0, il + i - 1) * MATRIX(vectors2, 0, i) < 0) {
            m1 = -1.0;
        } else {
            m1 = 1.0;
        }

        for (j = 0; j < DIM; j++) {
            if (fabs(MATRIX(vectors, j, il + i - 1) -
                     m1 * MATRIX(vectors2, j, i)) > 1e-8) {
                printf("Full:\n");
                igraph_matrix_print(&vectors);
                printf("Subset:\n");
                igraph_matrix_print(&vectors2);
                return 7;
            }
        }
    }

    igraph_vector_destroy(&values2);
    igraph_matrix_destroy(&vectors2);

    /* Subset based on an interval */

    igraph_matrix_init(&vectors2, 0, 0);
    igraph_vector_init(&values2, 0);

    il = 2;
    iu = 5;
    vl = (VECTOR(values)[il - 1] + VECTOR(values)[il - 2]) / 2.0;
    vu = (VECTOR(values)[iu] + VECTOR(values)[iu - 1]) / 2.0;

    igraph_lapack_dsyevr(&A, IGRAPH_LAPACK_DSYEV_INTERVAL, vl, vu,
                         /*vestimate=*/ iu - il + 1, /*il=*/ 0, /*iu=*/ 0,
                         /*abstol=*/ 1e-10, &values2, &vectors2,
                         /*support=*/ 0);

    if (igraph_vector_size(&values2) != iu - il + 1) {
        return 4;
    }
    if (igraph_matrix_nrow(&vectors2) != DIM ||
        igraph_matrix_ncol(&vectors2) != iu - il + 1) {
        return 5;
    }
    for (i = 0; i < iu - il + 1; i++) {
        igraph_real_t m1 = 1.0;

        if (fabs(VECTOR(values)[il + i - 1] - VECTOR(values2)[i]) > 1e-8) {
            printf("Full:   ");
            igraph_vector_print(&values);
            printf("Subset: ");
            igraph_vector_print(&values2);
            return 6;
        }

        if (MATRIX(vectors, 0, il + i - 1) * MATRIX(vectors2, 0, i) < 0) {
            m1 = -1.0;
        } else {
            m1 = 1.0;
        }

        for (j = 0; j < DIM; j++) {
            if (fabs(MATRIX(vectors, j, il + i - 1) -
                     m1 * MATRIX(vectors2, j, i)) > 1e-8) {
                printf("Full:\n");
                igraph_matrix_print(&vectors);
                printf("Subset:\n");
                igraph_matrix_print(&vectors2);
                return 7;
            }
        }
    }

    igraph_vector_destroy(&values2);
    igraph_matrix_destroy(&vectors2);

    igraph_vector_destroy(&values);
    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&A);

    return 0;
}
