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

void print_result(igraph_spmatrix_t *m, FILE *f) {
    print_spmatrix(m);
    fprintf(f, "==================================================\n");
}

int main() {
    igraph_spmatrix_t m, m1;
    igraph_spmatrix_iter_t mit;
    igraph_real_t arr[12];
    igraph_vector_t v;
    long int i, j;
    int order[] = { 1, 5, 8, 4, 0, 9, 6, 10, 11, 2, 3, 7 };

    /* igraph_spmatrix_init, igraph_spmatrix_destroy */
    igraph_spmatrix_init(&m, 10, 10);
    igraph_spmatrix_destroy(&m);

    igraph_spmatrix_init(&m, 0, 0);
    igraph_spmatrix_destroy(&m);

    /* igraph_spmatrix_ncol, igraph_spmatrix_nrow */
    igraph_spmatrix_init(&m, 10, 5);
    if (igraph_spmatrix_nrow(&m) != 10) {
        return 1;
    }
    if (igraph_spmatrix_ncol(&m) != 5) {
        return 2;
    }

    /* igraph_spmatrix_size, igraph_spmatrix_resize */
    igraph_spmatrix_resize(&m, 6, 5);
    if (igraph_spmatrix_size(&m) != 30) {
        return 3;
    }
    if (igraph_spmatrix_nrow(&m) != 6) {
        return 4;
    }
    if (igraph_spmatrix_ncol(&m) != 5) {
        return 5;
    }
    igraph_spmatrix_resize(&m, 2, 4);
    if (igraph_spmatrix_nrow(&m) != 2) {
        return 6;
    }
    if (igraph_spmatrix_ncol(&m) != 4) {
        return 7;
    }
    igraph_spmatrix_destroy(&m);

    /* igraph_spmatrix_get, igraph_spmatrix_set, igraph_spmatrix_null */
    igraph_spmatrix_init(&m, 3, 4);
    for (i = 0; i < igraph_spmatrix_nrow(&m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(&m); j++) {
            igraph_spmatrix_set(&m, i, j, (i + j) % 3);
        }
    }
    print_result(&m, stdout);
    igraph_spmatrix_null(&m);
    print_result(&m, stdout);
    /* now fill it in shuffled order */
    for (i = 0; i < 12; i++) {
        igraph_spmatrix_set(&m, order[i] / 4, order[i] % 4, (order[i] / 4 + order[i] % 4) % 3);
    }
    print_result(&m, stdout);
    /* now decrease all elements by two in shuffled order */
    for (i = 0; i < 12; i++) {
        igraph_spmatrix_add_e(&m, order[i] / 4, order[i] % 4, -2);
    }
    print_result(&m, stdout);
    /* now increase all elements by one in shuffled order */
    for (i = 0; i < 12; i++) {
        igraph_spmatrix_add_e(&m, order[i] / 4, order[i] % 4, 1);
    }
    print_result(&m, stdout);

    igraph_spmatrix_destroy(&m);

    /* igraph_matrix_add_cols, igraph_matrix_add_rows */
    igraph_spmatrix_init(&m, 4, 3);
    for (i = 0; i < igraph_spmatrix_nrow(&m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(&m); j++) {
            igraph_spmatrix_set(&m, i, j, (i + 1) * (j + 1));
        }
    }
    igraph_spmatrix_add_cols(&m, 2);
    igraph_spmatrix_add_rows(&m, 2);
    if (igraph_spmatrix_ncol(&m) != 5) {
        return 8;
    }
    if (igraph_spmatrix_nrow(&m) != 6) {
        return 9;
    }
    print_result(&m, stdout);
    igraph_spmatrix_destroy(&m);

    /* igraph_spmatrix_count_nonzero */
    igraph_spmatrix_init(&m, 5, 3);
    for (i = 0; i < igraph_spmatrix_nrow(&m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(&m); j++) {
            igraph_spmatrix_set(&m, i, j, i * j);
        }
    }
    print_result(&m, stdout);
    if (igraph_spmatrix_count_nonzero(&m) != 8) {
        return 10;
    }
    igraph_spmatrix_destroy(&m);

    /* igraph_spmatrix_copy */
    igraph_spmatrix_init(&m, 3, 4);
    for (i = 0; i < igraph_spmatrix_nrow(&m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(&m); j++) {
            igraph_spmatrix_set(&m, i, j, i * j);
        }
    }
    igraph_spmatrix_copy(&m1, &m);
    print_result(&m1, stdout);
    igraph_spmatrix_destroy(&m);
    igraph_spmatrix_destroy(&m1);

    /* igraph_spmatrix_copy_to */
    igraph_spmatrix_init(&m, 3, 4);
    for (i = 0; i < igraph_spmatrix_nrow(&m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(&m); j++) {
            igraph_spmatrix_set(&m, i, j, i * j);
        }
    }
    igraph_spmatrix_copy_to(&m, arr);
    for (i = 0; i < 12; i++) {
        printf(" %ld", (long)arr[i]);
    }
    printf("\n=========================\n");

    /* igraph_spmatrix_max */
    arr[0] = igraph_spmatrix_max(&m, arr + 1, arr + 2);
    for (i = 0; i < 3; i++) {
        printf(" %ld", (long)arr[i]);
    }
    printf("\n=========================\n");

    igraph_spmatrix_destroy(&m);

    /* igraph_spmatrix_colsums */
    igraph_spmatrix_init(&m, 3, 5);
    for (i = 0; i < igraph_spmatrix_nrow(&m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(&m); j++) {
            igraph_spmatrix_set(&m, i, j, i + j - 4);
        }
    }
    igraph_vector_init(&v, 0);
    igraph_spmatrix_colsums(&m, &v);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);
    igraph_spmatrix_destroy(&m);

    /* igraph_spmatrix_iter_t */
    igraph_spmatrix_init(&m, 5, 5);
    for (i = 0; i < igraph_spmatrix_nrow(&m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(&m); j++) {
            if (labs(i - j) == 1) {
                igraph_spmatrix_set(&m, i, j, (i + 1) * (j + 1));
            }
        }
    }
    igraph_spmatrix_iter_create(&mit, &m);
    while (!igraph_spmatrix_iter_end(&mit)) {
        printf("%ld %ld %ld\n", mit.ri, mit.ci, (long int)mit.value);
        igraph_spmatrix_iter_next(&mit);
    }
    igraph_spmatrix_iter_destroy(&mit);
    igraph_spmatrix_destroy(&m);
    printf("=========================\n");

    VERIFY_FINALLY_STACK();

    return 0;
}
