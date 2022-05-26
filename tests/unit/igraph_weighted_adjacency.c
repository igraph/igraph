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
    igraph_vector_t weights;

    igraph_vector_init(&weights, 0);

    printf("No loops\n--------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, &weights, IGRAPH_NO_LOOPS);
    print_weighted_graph(&g, &weights);
    igraph_destroy(&g);
    printf("\n");

    printf("Loops once\n----------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, &weights, IGRAPH_LOOPS_ONCE);
    print_weighted_graph(&g, &weights);
    igraph_destroy(&g);
    printf("\n");

    printf("Loops twice\n-----------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, &weights, IGRAPH_LOOPS_TWICE);
    print_weighted_graph(&g, &weights);
    igraph_destroy(&g);
    printf("\n");

    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();
}

int main() {
    igraph_matrix_t mat;
    igraph_vector_t weights;
    igraph_t g;
    int m[4][4] = { { 0, 1, 2, 0 }, { 2, 0, 0, 1 }, { 0, 0, 4, 0 }, { 0, 1, 0, 0 } };
    igraph_integer_t i, j;

    printf("0x0 matrix\n");
    printf("==========\n\n");
    igraph_matrix_init(&mat, 0, 0);
    test_single_matrix(&mat, IGRAPH_ADJ_DIRECTED);
    igraph_matrix_destroy(&mat);

    igraph_matrix_init(&mat, 4, 4);
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        MATRIX(mat, i, j) = m[i][j];
    }

    /* [ 0 1 2 0 ]
       [ 2 0 0 1 ]
       [ 0 0 4 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_DIRECTED\n");
    printf("===================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_DIRECTED);

    /* [ 0 1 2 0 ]
       [ - 0 0 1 ]
       [ - - 4 0 ]
       [ - - - 0 ] */
    printf("IGRAPH_ADJ_UPPER\n");
    printf("================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_UPPER);

    /* [ 0 - - - ]
       [ 2 0 - - ]
       [ 0 0 4 - ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_LOWER\n");
    printf("================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_LOWER);

    /* [ 0 1 0 0 ]
       [ 1 0 0 1 ]
       [ 0 0 4 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_MIN\n");
    printf("==============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_MIN);

    /* [ 0 2 2 0 ]
       [ 2 0 0 1 ]
       [ 2 0 4 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_MAX\n");
    printf("==============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_MAX);

    /* [ 0 3 2 0 ]
       [ 3 0 0 2 ]
       [ 2 0 4 0 ]
       [ 0 2 0 0 ] */
    printf("IGRAPH_ADJ_PLUS\n");
    printf("===============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_PLUS);

    igraph_matrix_destroy(&mat);

    VERIFY_FINALLY_STACK();

    igraph_vector_init(&weights, 0);

    {
        printf("Check handling of non-square matrix error.\n");
        igraph_real_t e[] = {1, 2, 0};
        igraph_matrix_view(&mat, e, 1, 3);
        CHECK_ERROR(igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_DIRECTED, &weights, IGRAPH_NO_LOOPS), IGRAPH_NONSQUARE);
    }

    {
        printf("Check handling of invalid adjacency mode.\n");
        igraph_real_t e[] = {0, 2, 0, 3, 0, 4, 0, 5, 6};
        igraph_matrix_view(&mat, e, 3, 3);
        CHECK_ERROR(igraph_weighted_adjacency(&g, &mat, 42, &weights, IGRAPH_LOOPS_TWICE), IGRAPH_EINVAL);
    }

    {
        printf("Check error for 0x1 matrix.\n");
        igraph_matrix_init(&mat, 0, 1);
        CHECK_ERROR(igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_DIRECTED, &weights, 1), IGRAPH_NONSQUARE);
        igraph_matrix_destroy(&mat);
    }

    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();

    return 0;
}
