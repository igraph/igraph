/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.inc"

void check_and_destroy(igraph_matrix_t *matrix, igraph_bool_t pivot) {
    igraph_vector_int_t pivot_indices;
    int info;
    igraph_vector_int_init(&pivot_indices, 0);
    printf("Starting matrix:\n");
    igraph_matrix_print(matrix);
    if (pivot) {
        IGRAPH_ASSERT(igraph_lapack_dgetrf(matrix, &pivot_indices, &info) == IGRAPH_SUCCESS);
    } else {
        IGRAPH_ASSERT(igraph_lapack_dgetrf(matrix, NULL, &info) == IGRAPH_SUCCESS);
    }
    printf("Returned matrix:\n");
    igraph_matrix_print(matrix);
    if (pivot) {
        printf("Returned pivot indices:\n");
        igraph_vector_int_print(&pivot_indices);
    }
    printf("info: %d\n", info);
    igraph_vector_int_destroy(&pivot_indices);
    igraph_matrix_destroy(matrix);
    printf("\n");
}

int main() {
    igraph_matrix_t matrix;

    printf("Empty matrix:\n");
    igraph_matrix_init(&matrix, 0, 0);
    check_and_destroy(&matrix, 1);

    int elements_1[9] = {7, 8, 9, 2, 2, 3, 1, 1, 1};
    matrix_init_int_row_major(&matrix, 3, 3, elements_1);
    check_and_destroy(&matrix, 1);

    int elements_2[9] = {1, 1, 1, 2, 2, 3, 7, 8, 9};
    matrix_init_int_row_major(&matrix, 3, 3, elements_2);
    check_and_destroy(&matrix, 1);

    int elements_3[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    matrix_init_int_row_major(&matrix, 3, 3, elements_3);
    check_and_destroy(&matrix, 0);

    int elements_4[6] = {1, 2, 3, 4, 5, 6};
    matrix_init_int_row_major(&matrix, 2, 3, elements_4);
    check_and_destroy(&matrix, 1);

    int elements_5[6] = {1, 2, 3, 4, 5, 6};
    matrix_init_int_row_major(&matrix, 3, 2, elements_5);
    check_and_destroy(&matrix, 1);

    VERIFY_FINALLY_STACK();
    return 0;
}
