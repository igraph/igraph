/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include <stdio.h>

#include "test_utilities.inc"

void byrow(igraph_matrix_t *m) {
    long int r = igraph_matrix_nrow(m), c = igraph_matrix_ncol(m);
    long int n = 0, i, j;
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            MATRIX(*m, i, j) = n++;
        }
    }
}

#define apply(m,a,b) \
    for (i=0; i<igraph_matrix_nrow(&(m)); i++) { \
        for (j=0; j<igraph_matrix_ncol(&(m)); j++) { \
            (a); \
        } \
        (b); \
    }


int main() {
    igraph_matrix_t m, m2;
    igraph_vector_t v;
    long int i, j, i2, j2;
    igraph_real_t r1, r2;

    igraph_matrix_init(&m, 4, 3);
    byrow(&m);

    /* igraph_matrix_e */
    printf("igraph_matrix_e\n");
    apply(m, printf("%i ", (int)igraph_matrix_e(&m, i, j)), printf("\n"));

    /* igraph_matrix_e_ptr */
    printf("igraph_matrix_e_ptr\n");
    apply(m, printf("%i ", (int)igraph_matrix_e_ptr(&m, i, j)[0]), printf("\n"));

    /* igraph_matrix_set */
    printf("igraph_matrix_set\n");
    apply(m, igraph_matrix_set(&m, i, j, i), (void) 0 );
    print_matrix(&m);
    apply(m, igraph_matrix_set(&m, i, j, j), (void) 0 );
    print_matrix(&m);

    /* igraph_matrix_fill */
    printf("igraph_matrix_fill\n");
    igraph_matrix_fill(&m, 42);
    print_matrix(&m);
    igraph_matrix_fill(&m, -42.1);
    print_matrix(&m);

    /* igraph_matrix_update */
    printf("igraph_matrix_update\n");
    igraph_matrix_init(&m2, 0, 0);
    byrow(&m);
    igraph_matrix_update(&m2, &m);
    print_matrix(&m2);

    /* igraph_matrix_rbind */
    printf("igraph_matrix_rbind\n");
    igraph_matrix_rbind(&m2, &m);
    print_matrix(&m2);
    printf("\n");
    igraph_matrix_resize(&m, 0, igraph_matrix_ncol(&m2));
    igraph_matrix_rbind(&m2, &m);
    print_matrix(&m2);
    printf("\n");
    igraph_matrix_rbind(&m, &m2);
    print_matrix(&m);

    /* igraph_matrix_cbind */
    printf("igraph_matrix_cbind\n");
    igraph_matrix_resize(&m, 4, 3);
    igraph_matrix_resize(&m2, 4, 2);
    byrow(&m);
    byrow(&m2);
    igraph_matrix_cbind(&m, &m2);
    print_matrix(&m);

    /* igraph_matrix_swap */
    printf("igraph_matrix_swap\n");
    igraph_matrix_update(&m, &m2);
    igraph_matrix_null(&m);
    igraph_matrix_swap(&m, &m2);
    print_matrix(&m);
    print_matrix(&m2);

    /* igraph_matrix_get_row */
    /* igraph_matrix_set_row */
    printf("igraph_matrix_get_row\n");
    printf("igraph_matrix_set_row\n");
    igraph_vector_init(&v, 0);
    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        igraph_matrix_get_row(&m, &v, i);
        igraph_matrix_set_row(&m2, &v, i);
    }
    print_matrix(&m2);

    /* igraph_matrix_set_col */
    printf("igraph_matrix_set_col\n");
    igraph_matrix_null(&m2);
    for (i = 0; i < igraph_matrix_ncol(&m); i++) {
        igraph_matrix_get_col(&m, &v, i);
        igraph_matrix_set_col(&m2, &v, i);
    }
    print_matrix(&m2);

    /* igraph_matrix_swap_rows */
    printf("igraph_matrix_swap_rows\n");
    igraph_matrix_swap_rows(&m2, 0, 0);
    igraph_matrix_swap_rows(&m2, 0, 2);
    print_matrix(&m2);

    /* igraph_matrix_swap_cols */
    printf("igraph_matrix_swap_cols\n");
    igraph_matrix_swap_cols(&m2, 0, 0);
    igraph_matrix_swap_cols(&m2, 0, 1);
    print_matrix(&m2);

    /* igraph_matrix_add_constant */
    printf("igraph_matrix_add_constant\n");
    igraph_matrix_add_constant(&m2, 0);
    print_matrix(&m2);
    igraph_matrix_add_constant(&m2, -1);
    print_matrix(&m2);

    /* igraph_matrix_add */
    printf("igraph_matrix_add\n");
    byrow(&m2);
    byrow(&m);
    igraph_matrix_add(&m2, &m);
    print_matrix(&m2);

    /* igraph_matrix_sub */
    printf("igraph_matrix_sub\n");
    igraph_matrix_sub(&m2, &m);
    print_matrix(&m2);

    /* igraph_matrix_mul_elements */
    printf("igraph_matrix_mul_elements\n");
    igraph_matrix_mul_elements(&m2, &m);
    print_matrix(&m2);

    /* igraph_matrix_div_elements */
    printf("igraph_matrix_div_elements\n");
    igraph_matrix_fill(&m, 2);
    igraph_matrix_div_elements(&m2, &m);
    print_matrix(&m2);

    /* igraph_matrix_min */
    printf("igraph_matrix_min\n");
    if (igraph_matrix_min(&m2) != 0) {
        return 1;
    }
    if (igraph_matrix_min(&m) != 2) {
        return 1;
    }

    /* igraph_matrix_which_min */
    printf("igraph_matrix_which_min\n");
    igraph_matrix_which_min(&m2, &i, &j);
    if (i != 0 || j != 0) {
        return 2;
    }
    MATRIX(m2, 0, 1) = -1;
    igraph_matrix_which_min(&m2, &i, &j);
    if (i != 0 || j != 1) {
        return 2;
    }
    MATRIX(m2, 3, 1) = -2;
    igraph_matrix_which_min(&m2, &i, &j);
    if (i != 3 || j != 1) {
        return 2;
    }

    /* igraph_matrix_which_max */
    printf("igraph_matrix_which_max\n");
    MATRIX(m2, 3, 0) = 100;
    igraph_matrix_which_max(&m2, &i, &j);
    if (i != 3 || j != 0) {
        return 3;
    }

    /* igraph_matrix_minmax */
    printf("igraph_matrix_minmax\n");
    igraph_matrix_minmax(&m2, &r1, &r2);
    printf("%g %g\n", r1, r2);

    /* igraph_matrix_which_minmax */
    printf("igraph_matrix_which_minmax\n");
    igraph_matrix_which_minmax(&m2, &i, &j, &i2, &j2);
    if (i != 3 || j != 1 || i2 != 3 || j2 != 0) {
        return 4;
    }

    /* igraph_matrix_isnull */
    printf("igraph_matrix_isnull\n");
    if (igraph_matrix_isnull(&m2)) {
        return 5;
    }
    igraph_matrix_null(&m);
    if (!igraph_matrix_isnull(&m)) {
        return 5;
    }
    igraph_matrix_resize(&m2, 5, 0);
    if (!igraph_matrix_isnull(&m2)) {
        return 5;
    }

    /* igraph_matrix_empty */
    printf("igraph_matrix_empty\n");
    if (!igraph_matrix_empty(&m2)) {
        return 6;
    }
    igraph_matrix_resize(&m2, 5, 5);
    if (igraph_matrix_empty(&m2)) {
        return 6;
    }

    /* igraph_matrix_is_symmetric */
    printf("igraph_matrix_is_symmetric\n");
    byrow(&m2);
    if (igraph_matrix_is_symmetric(&m2)) {
        return 7;
    }
    igraph_matrix_update(&m, &m2);
    igraph_matrix_transpose(&m);
    igraph_matrix_add(&m, &m2);
    if (!igraph_matrix_is_symmetric(&m)) {
        return 7;
    }

    /* igraph_matrix_prod */
    printf("igraph_matrix_prod\n");
    igraph_matrix_resize(&m, 3, 2);
    byrow(&m);
    igraph_matrix_add_constant(&m, 1);
    print_matrix(&m);
    printf("product: %g\n", igraph_matrix_prod(&m));

    /* igraph_matrix_rowsum */
    printf("igraph_matrix_rowsum\n");
    igraph_matrix_rowsum(&m, &v);
    print_vector(&v);

    /* igraph_matrix_colsum */
    printf("igraph_matrix_colsum\n");
    igraph_matrix_colsum(&m, &v);
    print_vector(&v);

    /* igraph_matrix_contains */
    printf("igraph_matrix_contains\n");
    if (igraph_matrix_contains(&m, 0)) {
        return 8;
    }
    if (igraph_matrix_contains(&m, 6.0001)) {
        return 8;
    }
    if (igraph_matrix_contains(&m, 7)) {
        return 8;
    }
    if (!igraph_matrix_contains(&m, 1)) {
        return 8;
    }
    if (!igraph_matrix_contains(&m, 6)) {
        return 8;
    }

    /* igraph_matrix_search */
    printf("igraph_matrix_search\n");
    if (!igraph_matrix_search(&m, 0, 6.0, &i2, &i, &j)) {
        return 9;
    }
    if (i2 != 5 || i != 2 || j != 1) {
        return 9;
    }

    /* igraph_matrix_remove_row */
    printf("igraph_matrix_remove_row\n");
    igraph_matrix_remove_row(&m, 1);
    print_matrix(&m);
    igraph_matrix_resize(&m, 5, 4);
    byrow(&m);
    igraph_matrix_remove_row(&m, 4);
    print_matrix(&m);
    igraph_matrix_remove_row(&m, 0);
    print_matrix(&m);

    /* igraph_matrix_select_cols */
    printf("igraph_matrix_select_cols\n");
    igraph_matrix_resize(&m, 6, 5);
    apply(m, igraph_matrix_set(&m, i, j, j), (void) 0 );
    igraph_vector_resize(&v, 3);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 4;
    VECTOR(v)[2] = 2;
    igraph_matrix_select_cols(&m, &m2, &v);
    print_matrix(&m2);
    igraph_vector_resize(&v, 1);
    igraph_matrix_select_cols(&m, &m2, &v);
    print_matrix(&m2);
    igraph_vector_clear(&v);
    igraph_matrix_select_cols(&m, &m2, &v);
    if (!igraph_matrix_empty(&m2)) {
        return 9;
    }

    igraph_vector_destroy(&v);
    igraph_matrix_destroy(&m2);
    igraph_matrix_destroy(&m);

    VERIFY_FINALLY_STACK();

    return 0;
}
