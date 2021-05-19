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

igraph_bool_t check_solution(const igraph_sparsemat_t *A,
                             const igraph_vector_t *x,
                             const igraph_vector_t *b) {

    long int dim = igraph_vector_size(x);
    igraph_vector_t res;
    int j, p;
    igraph_real_t min, max;

    igraph_vector_copy(&res, b);

    for (j = 0; j < dim; j++) {
        for (p = A->cs->p[j]; p < A->cs->p[j + 1]; p++) {
            long int from = A->cs->i[p];
            igraph_real_t value = A->cs->x[p];
            VECTOR(res)[from] -= VECTOR(*x)[j] * value;
        }
    }

    igraph_vector_minmax(&res, &min, &max);
    igraph_vector_destroy(&res);

    return fabs(min) < 1e-15 && fabs(max) < 1e-15;
}

int main() {

    igraph_sparsemat_t A, B, C;
    igraph_vector_t b, x;
    long int i;

    /* lsolve */

#define DIM 10
#define EDGES (DIM*DIM/6)
    igraph_sparsemat_init(&A, DIM, DIM, EDGES + DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, RNG_INTEGER(1, 3));
    }
    for (i = 0; i < EDGES; i++) {
        long int r = RNG_INTEGER(0, DIM - 1);
        long int c = RNG_INTEGER(0, r);
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
        long int r = RNG_INTEGER(0, DIM - 1);
        long int c = RNG_INTEGER(0, r);
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

    igraph_sparsemat_transpose(&B, &A, /*values=*/ 1);
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
        long int r = RNG_INTEGER(0, DIM - 1);
        long int c = RNG_INTEGER(0, r);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);
    igraph_sparsemat_transpose(&B, &A, /*values=*/ 1);

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
        long int r = RNG_INTEGER(0, DIM - 1);
        long int c = RNG_INTEGER(0, r);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, r, c, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);
    igraph_sparsemat_transpose(&B, &A, /*values=*/ 1);
    igraph_sparsemat_destroy(&B);

    igraph_vector_init(&b, DIM);
    for (i = 0; i < DIM; i++) {
        VECTOR(b)[i] = RNG_INTEGER(1, 10);
    }

    igraph_vector_init(&x, DIM);
    igraph_sparsemat_utsolve(&A, &b, &x);

    igraph_sparsemat_transpose(&A, &B, /*values=*/ 1);
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
        long int from = RNG_INTEGER(0, DIM - 1);
        long int to = RNG_INTEGER(0, DIM - 1);
        igraph_real_t value = RNG_INTEGER(1, 5);
        igraph_sparsemat_entry(&A, from, to, value);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_dupl(&B);
    igraph_sparsemat_transpose(&B, &A, /*values=*/ 1);
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
        long int from = RNG_INTEGER(0, DIM - 1);
        long int to = RNG_INTEGER(0, DIM - 1);
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
