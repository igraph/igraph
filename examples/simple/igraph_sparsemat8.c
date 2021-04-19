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

#define DIM1 10
#define DIM2 5

#define INT(a) (igraph_rng_get_integer(igraph_rng_default(), 0, (a)))

int main() {
    igraph_matrix_t mat, mat2;
    igraph_sparsemat_t spmat, spmat2;
    int i, j, nz1, nz2;
    igraph_vector_t sums1, sums2;


    igraph_rng_seed(igraph_rng_default(), 42);

    /* COPY */

    igraph_sparsemat_init(&spmat, DIM1, DIM2, 20);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&spmat, INT(DIM1 - 1), INT(DIM2 - 1), 1.0);
    }
    igraph_sparsemat_copy(&spmat2, &spmat);

    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_as_matrix(&mat, &spmat);
    igraph_matrix_init(&mat2, 0, 0);
    igraph_sparsemat_as_matrix(&mat2, &spmat2);
    if (!igraph_matrix_all_e(&mat, &mat2)) {
        return 1;
    }

    igraph_matrix_destroy(&mat2);
    igraph_sparsemat_destroy(&spmat2);

    igraph_sparsemat_compress(&spmat, &spmat2);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_copy(&spmat, &spmat2);

    igraph_matrix_init(&mat2, 0, 0);
    igraph_sparsemat_as_matrix(&mat2, &spmat);
    if (!igraph_matrix_all_e(&mat, &mat2)) {
        return 2;
    }

    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat2);
    igraph_matrix_destroy(&mat);
    igraph_matrix_destroy(&mat2);

    /* COLSUMS, ROWSUMS */

    igraph_sparsemat_init(&spmat, DIM1, DIM2, 20);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&spmat, INT(DIM1 - 1), INT(DIM2 - 1), 1.0);
    }
    igraph_sparsemat_compress(&spmat, &spmat2);

    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_as_matrix(&mat, &spmat);
    igraph_vector_init(&sums1, 0);
    igraph_vector_init(&sums2, 0);
    igraph_sparsemat_colsums(&spmat, &sums1);
    igraph_matrix_colsum(&mat, &sums2);
    if (!igraph_vector_all_e(&sums1, &sums2)) {
        return 3;
    }
    igraph_sparsemat_colsums(&spmat2, &sums1);
    if (!igraph_vector_all_e(&sums1, &sums2)) {
        return 4;
    }

    igraph_sparsemat_rowsums(&spmat, &sums1);
    igraph_matrix_rowsum(&mat, &sums2);
    if (!igraph_vector_all_e(&sums1, &sums2)) {
        return 5;
    }
    igraph_sparsemat_rowsums(&spmat2, &sums1);
    if (!igraph_vector_all_e(&sums1, &sums2)) {
        return 6;
    }

    igraph_matrix_destroy(&mat);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat2);
    igraph_vector_destroy(&sums1);
    igraph_vector_destroy(&sums2);

    /* COUNT_NONZERO, COUNT_NONZEROTOL */

    igraph_sparsemat_init(&spmat, DIM1, DIM2, 20);
    igraph_sparsemat_entry(&spmat, 1, 2, 1.0);
    igraph_sparsemat_entry(&spmat, 1, 2, 1.0);
    igraph_sparsemat_entry(&spmat, 1, 3, 1e-12);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&spmat, INT(DIM1 - 1), INT(DIM2 - 1), 1.0);
    }
    igraph_sparsemat_compress(&spmat, &spmat2);

    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_as_matrix(&mat, &spmat2);

    nz1 = igraph_sparsemat_count_nonzero(&spmat2);
    for (nz2 = 0, i = 0; i < igraph_matrix_nrow(&mat); i++) {
        for (j = 0; j < igraph_matrix_ncol(&mat); j++) {
            if (MATRIX(mat, i, j) != 0) {
                nz2++;
            }
        }
    }
    if (nz1 != nz2) {
        printf("%i %i\n", nz1, nz2);
        return 7;
    }

    nz1 = igraph_sparsemat_count_nonzerotol(&spmat2, 1e-10);
    for (nz2 = 0, i = 0; i < igraph_matrix_nrow(&mat); i++) {
        for (j = 0; j < igraph_matrix_ncol(&mat); j++) {
            if (fabs(MATRIX(mat, i, j)) >= 1e-10) {
                nz2++;
            }
        }
    }
    if (nz1 != nz2) {
        printf("%i %i\n", nz1, nz2);
        return 8;
    }

    igraph_matrix_destroy(&mat);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat2);

    /* SCALE */

    igraph_sparsemat_init(&spmat, DIM1, DIM2, 20);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&spmat, INT(DIM1 - 1), INT(DIM2 - 1), 1.0);
    }
    igraph_sparsemat_compress(&spmat, &spmat2);

    igraph_sparsemat_scale(&spmat, 2.0);
    igraph_sparsemat_scale(&spmat2, 2.0);
    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_as_matrix(&mat, &spmat);
    igraph_matrix_init(&mat2, 0, 0);
    igraph_sparsemat_as_matrix(&mat2, &spmat2);
    igraph_matrix_scale(&mat, 1.0 / 2.0);
    igraph_matrix_scale(&mat2, 1.0 / 2.0);
    if (!igraph_matrix_all_e(&mat, &mat2)) {
        return 9;
    }

    igraph_matrix_destroy(&mat);
    igraph_matrix_destroy(&mat2);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat2);

    /* ADDROWS, ADDCOLS */

    igraph_sparsemat_init(&spmat, DIM1, DIM2, 20);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&spmat, INT(DIM1 - 1), INT(DIM2 - 1), 1.0);
    }
    igraph_sparsemat_compress(&spmat, &spmat2);

    igraph_sparsemat_add_rows(&spmat, 3);
    igraph_sparsemat_add_cols(&spmat, 2);

    igraph_sparsemat_add_rows(&spmat2, 3);
    igraph_sparsemat_add_cols(&spmat2, 2);

    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_as_matrix(&mat, &spmat);
    igraph_matrix_init(&mat2, 0, 0);
    igraph_sparsemat_as_matrix(&mat2, &spmat2);
    if (!igraph_matrix_all_e(&mat, &mat2)) {
        return 10;
    }

    igraph_matrix_destroy(&mat);
    igraph_matrix_destroy(&mat2);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat2);

    return 0;
}
