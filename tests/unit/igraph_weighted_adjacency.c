/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include <stdarg.h>
#include "test_utilities.h"

void print_and_destroy(const igraph_matrix_t *adjmatrix,
        igraph_adjacency_t mode, igraph_bool_t loops) {
    igraph_vector_int_t el;
    igraph_vector_t weights;
    igraph_t g;
    igraph_integer_t i, j, n;
    char ch = (mode == IGRAPH_ADJ_DIRECTED) ? '>' : '-';

    igraph_vector_int_init(&el, 0);
    igraph_vector_init(&weights, 0);

    igraph_weighted_adjacency(&g, adjmatrix, mode, &weights, loops);

    igraph_get_edgelist(&g, &el, 0);
    n = igraph_ecount(&g);

    for (i = 0, j = 0; i < n; i++, j += 2) {
        printf("%" IGRAPH_PRId " --%c %" IGRAPH_PRId ": %g\n",
               VECTOR(el)[j], ch, VECTOR(el)[j + 1], VECTOR(weights)[i]);
    }
    printf("\n");

    igraph_vector_int_destroy(&el);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);
}

int main() {
    igraph_matrix_t mat;
    int m[4][4] = { { 0, 1, 2, 0 }, { 2, 0, 0, 1 }, { 0, 0, 1, 0 }, { 0, 1, 0, 0 } };
    igraph_integer_t i, j;

    printf("Use 0x0 matrix:\n");
    igraph_matrix_init(&mat, 0, 0);
    print_and_destroy(&mat, IGRAPH_ADJ_DIRECTED, /*loops=*/ 1);
    igraph_matrix_destroy(&mat);

    igraph_matrix_init(&mat, 4, 4);
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        MATRIX(mat, i, j) = m[i][j];
    }

    printf("Use all modes of 4x4 matrix:\n");
    /* [ 0 1 2 0 ]
       [ 2 0 0 1 ]
       [ 0 0 1 0 ]
       [ 0 1 0 0 ] */
    print_and_destroy(&mat, IGRAPH_ADJ_DIRECTED, /*loops=*/ 1);

    /* [ 0 1 2 0 ]
       [ - 0 0 1 ]
       [ - - 1 0 ]
       [ - - - 0 ] */
    print_and_destroy(&mat, IGRAPH_ADJ_UPPER, /*loops=*/ 1);

    /* [ 0 - - - ]
       [ 2 0 - - ]
       [ 0 0 1 - ]
       [ 0 1 0 0 ] */
    print_and_destroy(&mat, IGRAPH_ADJ_LOWER, /*loops=*/ 1);

    /* [ 0 1 0 0 ]
       [ 1 0 0 1 ]
       [ 0 0 1 0 ]
       [ 0 1 0 0 ] */
    print_and_destroy(&mat, IGRAPH_ADJ_MIN, /*loops=*/ 1);

    /* [ 0 2 2 0 ]
       [ 2 0 0 1 ]
       [ 2 0 1 0 ]
       [ 0 1 0 0 ] */
    print_and_destroy(&mat, IGRAPH_ADJ_MAX, /*loops=*/ 1);

    /* [ 0 3 2 0 ]
       [ 3 0 0 2 ]
       [ 2 0 1 0 ]
       [ 0 2 0 0 ] */
    print_and_destroy(&mat, IGRAPH_ADJ_PLUS, /*loops=*/ 1);

    igraph_matrix_destroy(&mat);

    VERIFY_FINALLY_STACK();

    {
        printf("Check error for 0x1 matrix.\n");
        igraph_t g;
        igraph_matrix_init(&mat, 0, 1);
        CHECK_ERROR(igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_DIRECTED, NULL, 1), IGRAPH_NONSQUARE);
        igraph_matrix_destroy(&mat);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
