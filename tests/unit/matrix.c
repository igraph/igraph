/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_matrix_t m, m1;
    long int i, j, k;

    /* igraph_matrix_init, igraph_matrix_destroy */
    igraph_matrix_init(&m, 10, 10);
    igraph_matrix_destroy(&m);

    igraph_matrix_init(&m, 0, 0);
    igraph_matrix_destroy(&m);

    /* igraph_matrix_ncol, igraph_matrix_nrow */
    igraph_matrix_init(&m, 10, 5);
    if (igraph_matrix_nrow(&m) != 10) {
        return 1;
    }
    if (igraph_matrix_ncol(&m) != 5) {
        return 2;
    }

    /* igraph_matrix_size, igraph_matrix_resize */
    igraph_matrix_resize(&m, 6, 5);
    if (igraph_matrix_size(&m) != 30) {
        return 3;
    }
    if (igraph_matrix_nrow(&m) != 6) {
        return 4;
    }
    if (igraph_matrix_ncol(&m) != 5) {
        return 5;
    }
    igraph_matrix_resize(&m, 2, 4);
    if (igraph_matrix_nrow(&m) != 2) {
        return 6;
    }
    if (igraph_matrix_ncol(&m) != 4) {
        return 7;
    }
    igraph_matrix_destroy(&m);

    /* MATRIX, igraph_matrix_null */
    igraph_matrix_init(&m, 3, 4);
    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        for (j = 0; j < igraph_matrix_ncol(&m); j++) {
            MATRIX(m, i, j) = i + 1;
        }
    }
    print_matrix(&m);
    igraph_matrix_null(&m);
    print_matrix(&m);
    igraph_matrix_destroy(&m);

    /* igraph_matrix_add_cols, igraph_matrix_add_rows */
    igraph_matrix_init(&m, 4, 3);
    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        for (j = 0; j < igraph_matrix_ncol(&m); j++) {
            MATRIX(m, i, j) = (i + 1) * (j + 1);
        }
    }
    igraph_matrix_add_cols(&m, 2);
    igraph_matrix_add_rows(&m, 2);
    if (igraph_matrix_ncol(&m) != 5) {
        return 8;
    }
    if (igraph_matrix_nrow(&m) != 6) {
        return 9;
    }
    igraph_matrix_destroy(&m);

    /* igraph_matrix_remove_col */
    igraph_matrix_init(&m, 5, 3);
    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        for (j = 0; j < igraph_matrix_ncol(&m); j++) {
            MATRIX(m, i, j) = (i + 1) * (j + 1);
        }
    }
    igraph_matrix_remove_col(&m, 0);
    print_matrix(&m);
    igraph_matrix_remove_col(&m, 1);
    print_matrix(&m);
    igraph_matrix_destroy(&m);

    /* TODO: igraph_matrix_permdelete_rows */
    /* TODO: igraph_matrix_delete_rows_neg */

    /* igraph_matrix_copy */
    igraph_matrix_init(&m, 2, 3);
    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        for (j = 0; j < igraph_matrix_ncol(&m); j++) {
            MATRIX(m, i, j) = (i + 1) * (j + 1);
        }
    }
    igraph_matrix_copy(&m1, &m);
    print_matrix(&m1);
    igraph_matrix_destroy(&m);
    igraph_matrix_destroy(&m1);

    /* in-place transpose */
    igraph_matrix_init(&m, 5, 2);
    k = 0;
    for (i = 0; i < igraph_matrix_ncol(&m); i++) {
        for (j = 0; j < igraph_matrix_nrow(&m); j++) {
            MATRIX(m, j, i) = k++;
        }
    }
    print_matrix(&m);
    igraph_matrix_transpose(&m);
    print_matrix(&m);
    igraph_matrix_destroy(&m);

    igraph_matrix_init(&m, 5, 1);
    k = 0;
    for (i = 0; i < igraph_matrix_ncol(&m); i++) {
        for (j = 0; j < igraph_matrix_nrow(&m); j++) {
            MATRIX(m, j, i) = k++;
        }
    }
    print_matrix(&m);
    igraph_matrix_transpose(&m);
    print_matrix(&m);
    igraph_matrix_destroy(&m);

    igraph_matrix_init(&m, 1, 5);
    k = 0;
    for (i = 0; i < igraph_matrix_ncol(&m); i++) {
        for (j = 0; j < igraph_matrix_nrow(&m); j++) {
            MATRIX(m, j, i) = k++;
        }
    }
    print_matrix(&m);
    igraph_matrix_transpose(&m);
    print_matrix(&m);
    igraph_matrix_destroy(&m);

    VERIFY_FINALLY_STACK();

    return 0;
}
