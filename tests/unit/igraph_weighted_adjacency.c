/* -*- mode: C -*-  */
/*
    IGraph library.
    Copyright (C) 2006-2022  The igraph development team

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
#include <stdarg.h>

#include "test_utilities.h"

void test_single_matrix(const igraph_matrix_t* mat, igraph_adjacency_t mode) {
    igraph_t g;

    printf("No loops\n--------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, 0, IGRAPH_NO_LOOPS);
    print_weighted_graph_attr(&g, "weight");
    igraph_destroy(&g);
    printf("\n");

    printf("Loops once\n----------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, 0, IGRAPH_LOOPS_ONCE);
    print_weighted_graph_attr(&g, "weight");
    igraph_destroy(&g);
    printf("\n");

    printf("Loops twice\n-----------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, 0, IGRAPH_LOOPS_TWICE);
    print_weighted_graph_attr(&g, "weight");
    igraph_destroy(&g);
    printf("\n");

    VERIFY_FINALLY_STACK();
}

int main() {
    igraph_matrix_t mat;
    igraph_t g;
    int m[4][4] = { { 0, 1, 2, 0 }, { 2, 0, 0, 1 }, { 0, 0, 4, 0 }, { 0, 1, 0, 0 } };
    igraph_integer_t i, j;

    igraph_matrix_init(&mat, 4, 4);
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        MATRIX(mat, i, j) = m[i][j];
    }
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* [ 0 1 2 0 ]
       [ 2 0 0 1 ]
       [ 0 0 1 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_DIRECTED\n");
    printf("===================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_DIRECTED);

    /* [ 0 1 2 0 ]
       [ - 0 0 1 ]
       [ - - 1 0 ]
       [ - - - 0 ] */
    printf("IGRAPH_ADJ_UPPER\n");
    printf("================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_UPPER);

    /* [ 0 - - - ]
       [ 2 0 - - ]
       [ 0 0 1 - ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_LOWER\n");
    printf("================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_LOWER);

    /* [ 0 1 0 0 ]
       [ 1 0 0 1 ]
       [ 0 0 1 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_MIN\n");
    printf("==============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_MIN);

    /* [ 0 2 2 0 ]
       [ 2 0 0 1 ]
       [ 2 0 1 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_MAX\n");
    printf("==============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_MAX);

    /* [ 0 3 2 0 ]
       [ 3 0 0 2 ]
       [ 2 0 1 0 ]
       [ 0 2 0 0 ] */
    printf("IGRAPH_ADJ_PLUS\n");
    printf("===============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_PLUS);

    igraph_matrix_destroy(&mat);

    VERIFY_FINALLY_STACK();

    printf("\nCheck handling of non-square matrix error.\n");
    {
        igraph_real_t e[] = {1, 2, 0};
        igraph_matrix_view(&mat, e, 1, 3);
        CHECK_ERROR(igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_DIRECTED, 0, IGRAPH_NO_LOOPS), IGRAPH_NONSQUARE);
    }
    printf("\nCheck handling of invalid adjacency mode.\n");
    {
        igraph_real_t e[] = {0, 2, 0, 3, 0, 4, 0, 5, 6};
        igraph_matrix_view(&mat, e, 3, 3);
        CHECK_ERROR(igraph_weighted_adjacency(&g, &mat, 42, 0, IGRAPH_LOOPS_TWICE), IGRAPH_EINVAL);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
