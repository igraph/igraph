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

#define EPS 1e-13


/* Generic test for 1x1 matrices */
void test_1x1(igraph_real_t value) {
    igraph_sparsemat_t A, B;
    igraph_matrix_t values, vectors;
    igraph_vector_t values2;
    igraph_arpack_options_t options;

    igraph_arpack_options_init(&options);

    igraph_sparsemat_init(&A, 1, 1, 1);
    igraph_sparsemat_entry(&A, 0, 0, value);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_matrix_init(&values, 0, 0);
    igraph_matrix_init(&vectors, 0, 0);
    options.mode = 1;
    igraph_sparsemat_arpack_rnsolve(&B, &options, /*storage=*/ 0,
                                    &values, &vectors);
    printf("rnsolve:\n  - eigenvalues:\n    ");
    igraph_matrix_print(&values);
    printf("  - eigenvectors:\n    ");
    igraph_matrix_print(&vectors);
    igraph_matrix_destroy(&values);
    igraph_matrix_destroy(&vectors);

    igraph_vector_init(&values2, 0);
    igraph_matrix_init(&vectors, 0, 0);
    options.mode = 1;
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values2, &vectors, IGRAPH_SPARSEMAT_SOLVE_LU);
    printf("rssolve:\n  - eigenvalues:\n    ");
    igraph_vector_print(&values2);
    printf("  - eigenvectors:\n    ");
    igraph_matrix_print(&vectors);
    igraph_vector_destroy(&values2);
    igraph_matrix_destroy(&vectors);

    igraph_sparsemat_destroy(&B);
}

/* Generic test for 2x2 matrices */
void test_2x2(igraph_real_t a, igraph_real_t b, igraph_real_t c, igraph_real_t d) {
    igraph_sparsemat_t A, B;
    igraph_matrix_t values, vectors;
    igraph_vector_t values2;
    igraph_arpack_options_t options;

    igraph_arpack_options_init(&options);
    options.mode = 1;
    options.nev = 2;

    igraph_sparsemat_init(&A, 2, 2, 4);
    igraph_sparsemat_entry(&A, 0, 0, a);
    igraph_sparsemat_entry(&A, 0, 1, b);
    igraph_sparsemat_entry(&A, 1, 0, c);
    igraph_sparsemat_entry(&A, 1, 1, d);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_matrix_init(&values, 0, 0);
    igraph_matrix_init(&vectors, 0, 0);
    igraph_sparsemat_arpack_rnsolve(&B, &options, /*storage=*/ 0,
                                    &values, &vectors);
    printf("rnsolve:\n  - eigenvalues:\n    ");
    igraph_matrix_print(&values);
    printf("  - eigenvectors:\n    ");
    igraph_matrix_print(&vectors);
    igraph_matrix_destroy(&values);
    igraph_matrix_destroy(&vectors);

    if (b == c) {
        igraph_vector_init(&values2, 0);
        igraph_matrix_init(&vectors, 0, 0);
        igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                        &values2, &vectors, IGRAPH_SPARSEMAT_SOLVE_QR);
        printf("rssolve:\n  - eigenvalues:\n    ");
        igraph_vector_print(&values2);
        printf("  - eigenvectors:\n    ");
        igraph_matrix_print(&vectors);
        igraph_vector_destroy(&values2);
        igraph_matrix_destroy(&vectors);
    }

    igraph_sparsemat_destroy(&B);
}

int main() {

    igraph_sparsemat_t A, B;
    igraph_matrix_t vectors, values2;
    igraph_vector_t values;
    long int i;
    igraph_arpack_options_t options;
    igraph_real_t min, max;
    igraph_t g1, g2, g3;

    /***********************************************************************/

    /* Identity matrix */
#define DIM 10
    igraph_sparsemat_init(&A, DIM, DIM, DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, 1.0);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_vector_init(&values, 0);
    igraph_arpack_options_init(&options);

    options.mode = 1;
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, /*vectors=*/ 0, /*solvemethod=*/0);
    if (VECTOR(values)[0] != 1.0) {
        return 1;
    }

    options.mode = 3;
    options.sigma = 2;
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, /*vectors=*/ 0,
                                    IGRAPH_SPARSEMAT_SOLVE_LU);
    if (VECTOR(values)[0] != 1.0) {
        return 21;
    }
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, /*vectors=*/ 0,
                                    IGRAPH_SPARSEMAT_SOLVE_QR);
    if (VECTOR(values)[0] != 1.0) {
        return 31;
    }

    igraph_vector_destroy(&values);
    igraph_sparsemat_destroy(&B);

#undef DIM

    /***********************************************************************/

    /* Diagonal matrix */
#define DIM 10
    igraph_sparsemat_init(&A, DIM, DIM, DIM);
    for (i = 0; i < DIM; i++) {
        igraph_sparsemat_entry(&A, i, i, i + 1.0);
    }
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_vector_init(&values, 0);
    igraph_matrix_init(&vectors, 0, 0);

    options.mode = 1;
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, /*vectors=*/ &vectors,
                                    /*solvemethod=*/ 0);
    if ( fabs(VECTOR(values)[0] - DIM) > EPS ) {
        printf("VECTOR(values)[0] numerical precision is only %g, should be %g",
               fabs((double)VECTOR(values)[0] - DIM), EPS);
        return 2;
    }

    if ( fabs(fabs(MATRIX(vectors, DIM - 1, 0)) - 1.0) > EPS) {
        return 3;
    }
    MATRIX(vectors, DIM - 1, 0) = 0.0;
    igraph_matrix_minmax(&vectors, &min, &max);
    if (fabs(min) > EPS) {
        return 3;
    }
    if (fabs(max) > EPS) {
        return 3;
    }

    options.mode = 3;
    options.sigma = 11;
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, /*vectors=*/ &vectors,
                                    IGRAPH_SPARSEMAT_SOLVE_LU);
    if ( fabs(VECTOR(values)[0] - DIM) > EPS ) {
        printf("VECTOR(values)[0] numerical precision is only %g, should be %g",
               fabs((double)VECTOR(values)[0] - DIM), EPS);
        return 22;
    }
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, /*vectors=*/ &vectors,
                                    IGRAPH_SPARSEMAT_SOLVE_QR);
    if ( fabs(VECTOR(values)[0] - DIM) > EPS ) {
        printf("VECTOR(values)[0] numerical precision is only %g, should be %g",
               fabs((double)VECTOR(values)[0] - DIM), EPS);
        return 32;
    }

    if ( fabs(fabs(MATRIX(vectors, DIM - 1, 0)) - 1.0) > EPS) {
        return 23;
    }
    MATRIX(vectors, DIM - 1, 0) = 0.0;
    igraph_matrix_minmax(&vectors, &min, &max);
    if (fabs(min) > EPS) {
        return 23;
    }
    if (fabs(max) > EPS) {
        return 23;
    }

    igraph_vector_destroy(&values);
    igraph_matrix_destroy(&vectors);
    igraph_sparsemat_destroy(&B);
#undef DIM

    /***********************************************************************/

    /* A tree, plus a ring */
#define DIM 10
    igraph_tree(&g1, DIM, /*children=*/ 2, IGRAPH_TREE_UNDIRECTED);
    igraph_ring(&g2, DIM, IGRAPH_UNDIRECTED, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_union(&g3, &g1, &g2, /*edge_map1=*/ 0, /*edge_map1=*/ 0);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    igraph_get_sparsemat(&g3, &A);
    igraph_destroy(&g3);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_vector_init(&values, 0);
    igraph_matrix_init(&vectors, 0, 0);

    options.mode = 1;
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, &vectors, /*solvemethod=*/ 0);

    if (MATRIX(vectors, 0, 0) < 0.0) {
        igraph_matrix_scale(&vectors, -1.0);
    }

    igraph_vector_print(&values);
    igraph_matrix_print(&vectors);

    options.mode = 3;
    options.sigma = VECTOR(values)[0] * 1.1;
    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, &vectors,
                                    IGRAPH_SPARSEMAT_SOLVE_LU);

    if (MATRIX(vectors, 0, 0) < 0.0) {
        igraph_matrix_scale(&vectors, -1.0);
    }
    igraph_vector_print(&values);
    igraph_matrix_print(&vectors);

    igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
                                    &values, &vectors,
                                    IGRAPH_SPARSEMAT_SOLVE_QR);
    if (MATRIX(vectors, 0, 0) < 0.0) {
        igraph_matrix_scale(&vectors, -1.0);
    }
    igraph_vector_print(&values);
    igraph_matrix_print(&vectors);

    igraph_vector_destroy(&values);
    igraph_matrix_destroy(&vectors);
    igraph_sparsemat_destroy(&B);
#undef DIM

    printf("--\n");

    /***********************************************************************/

    /* A directed tree and a directed, mutual ring */
#define DIM 10
    igraph_tree(&g1, DIM, /*children=*/ 2, IGRAPH_TREE_OUT);
    igraph_ring(&g2, DIM, IGRAPH_DIRECTED, /*mutual=*/ 1, /*circular=*/ 1);
    igraph_union(&g3, &g1, &g2, /*edge_map1=*/ 0, /*edge_map2=*/ 0);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    igraph_get_sparsemat(&g3, &A);
    igraph_destroy(&g3);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_matrix_init(&values2, 0, 0);
    igraph_matrix_init(&vectors, 0, 0);

    options.mode = 1;
    igraph_sparsemat_arpack_rnsolve(&B, &options, /*storage=*/ 0,
                                    &values2, &vectors);

    if (MATRIX(vectors, 0, 0) < 0.0) {
        igraph_matrix_scale(&vectors, -1.0);
    }

    igraph_matrix_print(&values2);
    igraph_matrix_print(&vectors);

    igraph_matrix_destroy(&values2);
    igraph_matrix_destroy(&vectors);
    igraph_sparsemat_destroy(&B);
#undef DIM

    /***********************************************************************/

    /* A small test graph */

    igraph_small(&g1, 11, IGRAPH_DIRECTED,
                 0, 1, 1, 3, 1, 8, 2, 10, 3, 6, 3, 10, 4, 2, 5, 4,
                 6, 1, 6, 4, 7, 9, 8, 5, 8, 7, 9, 8, 10, 0,
                 -1);

    igraph_get_sparsemat(&g1, &A);
    igraph_destroy(&g1);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&A);

    igraph_matrix_init(&values2, 0, 0);
    igraph_matrix_init(&vectors, 0, 0);

    options.mode = 1;
    igraph_sparsemat_arpack_rnsolve(&B, &options, /*storage=*/ 0,
                                    &values2, &vectors);

    if (MATRIX(vectors, 0, 0) < 0.0) {
        igraph_matrix_scale(&vectors, -1.0);
    }

    igraph_matrix_destroy(&values2);
    igraph_matrix_destroy(&vectors);
    igraph_sparsemat_destroy(&B);

    /***********************************************************************/

    /* Testing the special case solver for 1x1 matrices */
    printf("--\n");
    test_1x1(2);
    test_1x1(0);
    test_1x1(-3);

    /***********************************************************************/

    /* Testing the special case solver for 2x2 matrices */
    printf("--\n");
    test_2x2(1, 2, 2, 4);      /* symmetric */
    test_2x2(1, 2, 3, 4);      /* non-symmetric, real eigenvalues */
    test_2x2(1, -5, 10, 4);    /* non-symmetric, complex eigenvalues */
    test_2x2(0, 0, 0, 0);      /* symmetric, pathological */

    return 0;
}
