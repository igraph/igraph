/*
   igraph library.
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

igraph_bool_t check_solution(const igraph_sparsemat_t *A,
                             const igraph_vector_t *x,
                             const igraph_vector_t *b) {

    igraph_vector_t res;
    igraph_real_t min, max;
    igraph_bool_t success;
    igraph_sparsemat_iterator_t it;

    igraph_vector_init_copy(&res, b);

    igraph_sparsemat_iterator_init(&it, (igraph_sparsemat_t*) A);
    while (!igraph_sparsemat_iterator_end(&it)) {
        igraph_int_t row = igraph_sparsemat_iterator_row(&it);
        igraph_int_t col = igraph_sparsemat_iterator_col(&it);
        igraph_real_t value = igraph_sparsemat_iterator_get(&it);
        VECTOR(res)[row] -= VECTOR(*x)[col] * value;
        igraph_sparsemat_iterator_next(&it);
    }

    igraph_vector_minmax(&res, &min, &max);

    success = fabs(min) < 1e-12 && fabs(max) < 1e-12;

    if (!success) {
        printf("Incorrect solution.\n\n");
        printf("A =\n"); igraph_sparsemat_print(A, stdout); printf("\n");
        printf("x =\n"); igraph_vector_print(x); printf("\n");
        printf("b =\n"); igraph_vector_print(b); printf("\n");
        printf("difference between A*x and b =\n");
        igraph_vector_print(&res); printf("\n\n");
    }

    igraph_vector_destroy(&res);

    return success;
}

int main(void) {

    igraph_sparsemat_t A, B, C;
    igraph_vector_t b, x;
    igraph_int_t i;

    /* Initialize the library. */
    igraph_setup();

    /* Seed the RNG for reproducible results */
    igraph_rng_seed(igraph_rng_default(), 123);

    /* lsolve */

#define DIM 10
#define EDGES (DIM*DIM/6)
    igraph_sparsemat_init(&A, DIM, DIM, EDGES + DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, RNG_INTEGER(1, 3));
    }
    for (i = 0; i < EDGES; i++) {
        igraph_int_t r = RNG_INTEGER(0, DIM - 1);
        igraph_int_t c = RNG_INTEGER(0, r);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);

    igraph_vector_init(&b, DIM);
    for (i = 0; i < DIM; i++) {
        VECTOR(b)[i] = RNG_INTEGER(1, 10);
    }

    igraph_vector_init(&x, DIM);
    igraph_sparsemat_lsolve(&B, &b, &x);

    if (! check_solution(&B, &x, &b)) {
        return 1;
    }

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&x);
    igraph_sparsemat_destroy(&B);

#undef DIM
#undef EDGES

    /* ltsolve */

#define DIM 10
#define EDGES (DIM*DIM/6)
    igraph_sparsemat_init(&A, DIM, DIM, EDGES + DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, RNG_INTEGER(1, 3));
    }
    for (i = 0; i < EDGES; i++) {
        igraph_int_t r = RNG_INTEGER(0, DIM - 1);
        igraph_int_t c = RNG_INTEGER(0, r);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);

    igraph_vector_init(&b, DIM);
    for (i = 0; i < DIM; i++) {
        VECTOR(b)[i] = RNG_INTEGER(1, 10);
    }

    igraph_vector_init(&x, DIM);
    igraph_sparsemat_ltsolve(&B, &b, &x);

    igraph_sparsemat_transpose(&B, &A);
    if (! check_solution(&A, &x, &b)) {
        return 2;
    }

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&x);
    igraph_sparsemat_destroy(&B);
    igraph_sparsemat_destroy(&A);

#undef DIM
#undef EDGES

    /* usolve */

#define DIM 10
#define EDGES (DIM*DIM/6)
    igraph_sparsemat_init(&A, DIM, DIM, EDGES + DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, RNG_INTEGER(1, 3));
    }
    for (i = 0; i < EDGES; i++) {
        igraph_int_t r = RNG_INTEGER(0, DIM - 1);
        igraph_int_t c = RNG_INTEGER(0, r);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);
    igraph_sparsemat_transpose(&B, &A);

    igraph_vector_init(&b, DIM);
    for (i = 0; i < DIM; i++) {
        VECTOR(b)[i] = RNG_INTEGER(1, 10);
    }

    igraph_vector_init(&x, DIM);
    igraph_sparsemat_usolve(&A, &b, &x);

    if (! check_solution(&A, &x, &b)) {
        return 3;
    }

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&x);
    igraph_sparsemat_destroy(&B);
    igraph_sparsemat_destroy(&A);

#undef DIM
#undef EDGES

    /* utsolve */

#define DIM 10
#define EDGES (DIM*DIM/6)
    igraph_sparsemat_init(&A, DIM, DIM, EDGES + DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, RNG_INTEGER(1, 3));
    }
    for (i = 0; i < EDGES; i++) {
        igraph_int_t r = RNG_INTEGER(0, DIM - 1);
        igraph_int_t c = RNG_INTEGER(0, r);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);
    igraph_sparsemat_transpose(&B, &A);
    igraph_sparsemat_destroy(&B);

    igraph_vector_init(&b, DIM);
    for (i = 0; i < DIM; i++) {
        VECTOR(b)[i] = RNG_INTEGER(1, 10);
    }

    igraph_vector_init(&x, DIM);
    igraph_sparsemat_utsolve(&A, &b, &x);

    igraph_sparsemat_transpose(&A, &B);
    if (! check_solution(&B, &x, &b)) {
        return 4;
    }

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&x);
    igraph_sparsemat_destroy(&B);
    igraph_sparsemat_destroy(&A);

#undef DIM
#undef EDGES

    /* cholsol */
    /* We need a positive definite matrix, so we create a full-rank
       matrix first and then calculate A'A, which will be positive
       definite. */

#define DIM 10
#define EDGES (DIM*DIM/6)
    igraph_sparsemat_init(&A, DIM, DIM, EDGES + DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, RNG_INTEGER(1, 3));
    }
    for (i = 0; i < EDGES; i++) {
        igraph_int_t from = RNG_INTEGER(0, DIM - 1);
        igraph_int_t to = RNG_INTEGER(0, DIM - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, from, to, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);
    igraph_sparsemat_transpose(&B, &A);
    igraph_sparsemat_multiply(&A, &B, &C);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);

    igraph_vector_init(&b, DIM);
    for (i = 0; i < DIM; i++) {
        VECTOR(b)[i] = RNG_INTEGER(1, 10);
    }

    igraph_vector_init(&x, DIM);
    igraph_sparsemat_cholsol(&C, &b, &x, /*order=*/ 0);

    if (! check_solution(&C, &x, &b)) {
        return 5;
    }

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&x);
    igraph_sparsemat_destroy(&C);

#undef DIM
#undef EDGES

    /* lusol */

#define DIM 10
#define EDGES (DIM*DIM/4)
    igraph_sparsemat_init(&A, DIM, DIM, EDGES + DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, RNG_INTEGER(1, 3));
    }
    for (i = 0; i < EDGES; i++) {
        igraph_int_t from = RNG_INTEGER(0, DIM - 1);
        igraph_int_t to = RNG_INTEGER(0, DIM - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, from, to, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);

    igraph_vector_init(&b, DIM);
    for (i = 0; i < DIM; i++) {
        VECTOR(b)[i] = RNG_INTEGER(1, 10);
    }

    igraph_vector_init(&x, DIM);
    igraph_sparsemat_lusol(&B, &b, &x, /*order=*/ 0, /*tol=*/ 1e-10);

    if (! check_solution(&B, &x, &b)) {
        return 6;
    }

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&x);
    igraph_sparsemat_destroy(&B);

#undef DIM
#undef EDGES

    return 0;
}
