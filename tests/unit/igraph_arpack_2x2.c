/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

/* Tests for the special-cased 2x2 problems in igraph_arpack_rssolve()
 * and igraph_arpack_rnsolve(). */

#include <igraph.h>

#include "test_utilities.h"

#include <float.h>

#define EPS 1e-12

/* Multiplies by a 2x2 matrix given in row-major order. */
static igraph_error_t matvec(igraph_real_t *to, const igraph_real_t *from, int n, void *extra) {
    const igraph_real_t *m = extra;
    IGRAPH_UNUSED(n);
    to[0] = m[0] * from[0] + m[1] * from[1];
    to[1] = m[2] * from[0] + m[3] * from[1];
    return IGRAPH_SUCCESS;
}

/* Checks that (re + i*im, (x0r + i*x0i, x1r + i*x1i)) is an eigenpair of m
 * with a unit-length eigenvector. The residual is tested relative to the
 * largest entry of m. */
static void check_eigenpair(const igraph_real_t *m, igraph_real_t re, igraph_real_t im,
                            igraph_real_t x0r, igraph_real_t x0i,
                            igraph_real_t x1r, igraph_real_t x1i) {
    const igraph_real_t tol = EPS * fmax(fmax(fabs(m[0]), fabs(m[1])), fmax(fabs(m[2]), fabs(m[3])));
    /* residual m x - lambda x, real and imaginary parts */
    const igraph_real_t r0r = m[0] * x0r + m[1] * x1r - (re * x0r - im * x0i);
    const igraph_real_t r0i = m[0] * x0i + m[1] * x1i - (re * x0i + im * x0r);
    const igraph_real_t r1r = m[2] * x0r + m[3] * x1r - (re * x1r - im * x1i);
    const igraph_real_t r1i = m[2] * x0i + m[3] * x1i - (re * x1i + im * x1r);

    IGRAPH_ASSERT(fabs(r0r) < tol && fabs(r0i) < tol && fabs(r1r) < tol && fabs(r1i) < tol);
    IGRAPH_ASSERT(fabs(x0r * x0r + x0i * x0i + x1r * x1r + x1i * x1i - 1) < EPS);
}

/* Runs igraph_arpack_rssolve() with nev = 2 and checks the eigenvalues
 * (in the expected order) and that the eigenvectors are orthonormal. */
static void check_rs(igraph_real_t *m, const char *which, igraph_real_t l1, igraph_real_t l2) {
    igraph_arpack_options_t options;
    igraph_vector_t values;
    igraph_matrix_t vectors;

    igraph_arpack_options_init(&options);
    options.n = 2;
    options.nev = 2;
    options.which[0] = which[0]; options.which[1] = which[1];

    igraph_vector_init(&values, 0);
    igraph_matrix_init(&vectors, 0, 0);
    igraph_arpack_rssolve(matvec, m, &options, NULL, &values, &vectors);

    IGRAPH_ASSERT(fabs(VECTOR(values)[0] - l1) < EPS);
    IGRAPH_ASSERT(fabs(VECTOR(values)[1] - l2) < EPS);
    for (igraph_int_t j = 0; j < 2; j++) {
        check_eigenpair(m, VECTOR(values)[j], 0,
                        MATRIX(vectors, 0, j), 0, MATRIX(vectors, 1, j), 0);
    }
    IGRAPH_ASSERT(fabs(MATRIX(vectors, 0, 0) * MATRIX(vectors, 0, 1) +
                       MATRIX(vectors, 1, 0) * MATRIX(vectors, 1, 1)) < EPS);

    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
}

/* Same for igraph_arpack_rnsolve(), with real eigenvalues. */
static void check_rn_real(igraph_real_t *m, const char *which, igraph_real_t l1, igraph_real_t l2) {
    igraph_arpack_options_t options;
    igraph_matrix_t values, vectors;

    igraph_arpack_options_init(&options);
    options.n = 2;
    options.nev = 2;
    options.which[0] = which[0]; options.which[1] = which[1];

    igraph_matrix_init(&values, 0, 0);
    igraph_matrix_init(&vectors, 0, 0);
    igraph_arpack_rnsolve(matvec, m, &options, NULL, &values, &vectors);

    IGRAPH_ASSERT(fabs(MATRIX(values, 0, 0) - l1) < EPS && MATRIX(values, 0, 1) == 0);
    IGRAPH_ASSERT(fabs(MATRIX(values, 1, 0) - l2) < EPS && MATRIX(values, 1, 1) == 0);
    for (igraph_int_t j = 0; j < 2; j++) {
        check_eigenpair(m, MATRIX(values, j, 0), 0,
                        MATRIX(vectors, 0, j), 0, MATRIX(vectors, 1, j), 0);
    }

    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&values);
}

int main(void) {
    /* Diagonal matrices: the eigenvector of the larger eigenvalue is (0, 1). */
    igraph_real_t diag[] = { 0, 0, 0, 1 };

    /* Nearly equal diagonal entries. Used to give the same eigenvector twice,
     * which the orthogonality check catches, or a spurious error. */
    igraph_real_t near_diag[] = { 1 + DBL_EPSILON, 0, 0, 1 };

    /* Eigenvalues -3 and -1: LM and SM differ from LA and SA. */
    igraph_real_t sym[] = { -2, 1, 1, -2 };

    /* Non-symmetric, real eigenvalues 3 and 1. */
    igraph_real_t nonsym[] = { 1, 2, 0, 3 };

    /* Eigenvalues 1 and 0, prone to cancellation in the discriminant. */
    igraph_real_t cancel[] = { 1e8, -1e8, 1e8 - 1, 1 - 1e8 };

    /* Eigenvalues +-i. */
    igraph_real_t rot[] = { 0, -1, 1, 0 };

    /* Symmetric */
    check_rs(diag, "LA", 1, 0);
    check_rs(diag, "SA", 0, 1);
    check_rs(near_diag, "LA", 1 + DBL_EPSILON, 1);
    check_rs(sym, "LA", -1, -3);
    check_rs(sym, "LM", -3, -1);
    check_rs(sym, "SM", -1, -3);

    /* Non-symmetric, real eigenvalues */
    check_rn_real(diag, "LR", 1, 0);
    check_rn_real(diag, "SR", 0, 1);
    check_rn_real(nonsym, "LM", 3, 1);
    check_rn_real(nonsym, "SM", 1, 3);
    check_rn_real(cancel, "LR", 1, 0);

    /* Non-symmetric, complex eigenvalues */
    {
        igraph_arpack_options_t options;
        igraph_matrix_t values, vectors;

        igraph_arpack_options_init(&options);
        options.n = 2;
        options.nev = 2;
        options.which[0] = 'S'; options.which[1] = 'I';

        igraph_matrix_init(&values, 0, 0);
        igraph_matrix_init(&vectors, 0, 0);
        igraph_arpack_rnsolve(matvec, rot, &options, NULL, &values, &vectors);

        /* A conjugate pair, positive imaginary part first, and the eigenvector
         * for the first one, with real and imaginary parts in two columns. */
        IGRAPH_ASSERT(MATRIX(values, 0, 0) == 0 && fabs(MATRIX(values, 0, 1) - 1) < EPS);
        IGRAPH_ASSERT(MATRIX(values, 1, 0) == 0 && fabs(MATRIX(values, 1, 1) + 1) < EPS);
        check_eigenpair(rot, MATRIX(values, 0, 0), MATRIX(values, 0, 1),
                        MATRIX(vectors, 0, 0), MATRIX(vectors, 0, 1),
                        MATRIX(vectors, 1, 0), MATRIX(vectors, 1, 1));

        igraph_matrix_destroy(&vectors);
        igraph_matrix_destroy(&values);
    }

    /* nev = 1 and invalid input */
    {
        igraph_arpack_options_t options;
        igraph_vector_t values;
        igraph_matrix_t vectors;
        igraph_real_t nan_mat[] = { 1, 0, 0, IGRAPH_NAN };

        igraph_arpack_options_init(&options);
        options.n = 2;
        options.nev = 1;
        options.which[0] = 'L'; options.which[1] = 'A';

        igraph_vector_init(&values, 0);
        igraph_matrix_init(&vectors, 0, 0);

        igraph_arpack_rssolve(matvec, diag, &options, NULL, &values, &vectors);
        IGRAPH_ASSERT(igraph_vector_size(&values) == 1 && VECTOR(values)[0] == 1);
        IGRAPH_ASSERT(igraph_matrix_ncol(&vectors) == 1 && fabs(MATRIX(vectors, 1, 0)) == 1);

        /* LR is valid only for non-symmetric problems */
        options.which[0] = 'L'; options.which[1] = 'R';
        CHECK_ERROR(igraph_arpack_rssolve(matvec, diag, &options, NULL, &values, &vectors), IGRAPH_EARPACK);

        options.which[0] = 'L'; options.which[1] = 'A';
        CHECK_ERROR(igraph_arpack_rssolve(matvec, nan_mat, &options, NULL, &values, &vectors), IGRAPH_EINVAL);

        igraph_matrix_destroy(&vectors);
        igraph_vector_destroy(&values);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
