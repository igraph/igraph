/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "test_utilities.inc"

typedef struct cb2_data_t {
    igraph_matrix_t *A;
} cb2_data_t;

int cb2(igraph_real_t *to, const igraph_real_t *from, int n, void *extra) {
    cb2_data_t *data = (cb2_data_t*) extra;
    igraph_blas_dgemv_array(/*transpose=*/ 0, /*alpha=*/ 1.0,
                                           data->A, from, /*beta=*/ 0.0, to);
    return 0;
}

int check_eigenvector(
    const char* test_name,
    igraph_matrix_t* A, igraph_matrix_t* values, igraph_matrix_t* vectors,
    int eval_idx, int evec_col_idx
) {
    igraph_complex_t eval, prod;
    igraph_complex_t *evec;
    int i, j, n = igraph_matrix_nrow(A);

    eval = igraph_complex(MATRIX(*values, eval_idx, 0), MATRIX(*values, eval_idx, 1));
    evec = (igraph_complex_t*) calloc(n, sizeof(igraph_complex_t));
    if (IGRAPH_IMAG(eval) == 0) {
        /* Real eigenvalue, so we have a real eigenvector */
        for (i = 0; i < n; i++) {
            evec[i] = igraph_complex(MATRIX(*vectors, i, evec_col_idx), 0);
        }
    } else {
        /* Complex eigenvalue pair, so we have a complex eigenvector pair */
        /* ARPACK always stores the eigenvector corresponding to the eigenvalue
         * with a positive imaginary part. If the imaginary part is negative, we
         * need to multiply the imaginary part of the eigenvector by -1 */
        for (i = 0; i < n; i++) {
            evec[i] = igraph_complex(
                          MATRIX(*vectors, i, evec_col_idx),
                          MATRIX(*vectors, i, evec_col_idx + 1) * (
                              IGRAPH_IMAG(eval) < 0 ? -1 : 1
                          )
                      );
        }
    }

    /* Multiply matrix with eigenvector */
    for (i = 0; i < n; i++) {
        prod = igraph_complex(0, 0);
        for (j = 0; j < n; j++) {
            prod = igraph_complex_add(
                       igraph_complex_mul_real(evec[j], MATRIX(*A, i, j)),
                       prod
                   );
        }
        prod = igraph_complex_div(prod, eval);
        if (!igraph_complex_eq_tol(prod, evec[i], 1e-6)) {
            prod = igraph_complex_sub(prod, evec[i]);
            printf("%s: vector corresponding to eigenvalue (%.4f + %.4f*i) is not an "
                   "eigenvector, coordinate %d differs by %.4f + %.4f*i\n",
                   test_name, IGRAPH_REAL(eval), IGRAPH_IMAG(eval),
                   i, IGRAPH_REAL(prod), IGRAPH_IMAG(prod));
            return 1;
        }
    }

    /* Free stuff */
    free(evec);

    return 0;
}

int check_eigenvectors(
    const char* test_name,
    igraph_matrix_t* A, igraph_matrix_t* values, igraph_matrix_t* vectors
) {
    int i, j;
    int nev = igraph_matrix_nrow(values);
    int errors = 0;
    igraph_bool_t conjugate_pair_will_come = 0;

    for (i = 0, j = 0; i < nev; i++) {
        errors += check_eigenvector(test_name, A, values, vectors, i, j);
        if (MATRIX(*values, i, 1) != 0) {
            /* Complex eigenvalue */
            if (conjugate_pair_will_come) {
                j += 2;
                conjugate_pair_will_come = 0;
            } else {
                conjugate_pair_will_come = 1;
            }
        } else {
            /* Real eigenvalue */
            j++;
        }
    }
    return (errors > 0) ? 1 : 0;
}

void print_debug_output(
    igraph_matrix_t* values, igraph_matrix_t* vectors
) {
    printf("---\n");
    igraph_matrix_print(values);
    printf("---\n");
    igraph_matrix_print(vectors);
    printf("---\n");
}

#define DIM 10

int main() {
    igraph_matrix_t A;
    igraph_matrix_t values, vectors;
    igraph_arpack_options_t options;
    cb2_data_t data = { &A };
    int i, j;

    igraph_rng_seed(igraph_rng_default(), 42 * 42);

    igraph_matrix_init(&A, DIM, DIM);

    for (i = 0; i < DIM; i++) {
        for (j = 0; j < DIM; j++) {
            MATRIX(A, i, j) = igraph_rng_get_integer(igraph_rng_default(), -10, 10);
        }
    }

    igraph_matrix_print(&A);
    printf("===\n");

    igraph_arpack_options_init(&options);
    options.n = DIM;
    options.start = 0;
    options.nev = 4;
    options.ncv = 9;
    options.which[0] = 'L' ;
    options.which[1] = 'M';

    igraph_matrix_init(&values, 0, 0);
    igraph_matrix_init(&vectors, options.n, 1);

    igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0,
                          &values, &vectors);
    if (check_eigenvectors("LM #1", &A, &values, &vectors)) {
        print_debug_output(&values, &vectors);
    }

    /* -------------- */

    options.nev = 3;
    options.which[0] = 'L' ;
    options.which[1] = 'M';

    igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0,
                          &values, &vectors);
    if (check_eigenvectors("LM #2", &A, &values, &vectors)) {
        print_debug_output(&values, &vectors);
    }

    /* -------------- */

    options.nev = 3;
    options.which[0] = 'S' ;
    options.which[1] = 'R';

    igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0,
                          &values, &vectors);
    if (check_eigenvectors("SR", &A, &values, &vectors)) {
        print_debug_output(&values, &vectors);
    }

    /* -------------- */

    options.nev = 3;
    options.which[0] = 'L' ;
    options.which[1] = 'I';

    igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0,
                          &values, &vectors);
    if (check_eigenvectors("LI", &A, &values, &vectors)) {
        print_debug_output(&values, &vectors);
    }

    /* -------------- */

    igraph_matrix_destroy(&values);
    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&A);

    VERIFY_FINALLY_STACK();

    return 0;
}
