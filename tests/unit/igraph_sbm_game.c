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

void call_and_print(igraph_integer_t n, igraph_matrix_t *pref_matrix, igraph_vector_int_t *block_sizes, igraph_bool_t directed, igraph_bool_t loops) {
    igraph_t result;
    IGRAPH_ASSERT(igraph_sbm_game(&result, n, pref_matrix, block_sizes, directed, loops) == IGRAPH_SUCCESS);
    print_graph_canon(&result);
    printf("\n");
    igraph_destroy(&result);
}


int main() {
    igraph_t result;
    igraph_matrix_t pref_matrix_0, pref_matrix_1, pref_matrix_3, pref_matrix_3u, pref_matrix_nonsq, pref_matrix_oor, pref_matrix_nsym;
    igraph_vector_int_t block_sizes_0, block_sizes_1, block_sizes_3, block_sizes_neg;

    igraph_matrix_init(&pref_matrix_0, 0, 0);

    igraph_matrix_init(&pref_matrix_1, 1, 1);
    MATRIX(pref_matrix_1, 0, 0) = 1;

    igraph_matrix_init(&pref_matrix_3, 3, 3);
    igraph_matrix_null(&pref_matrix_3);
    MATRIX(pref_matrix_3, 0, 1) = 1;
    MATRIX(pref_matrix_3, 2, 2) = 1;

    igraph_matrix_init(&pref_matrix_3u, 3, 3);
    igraph_matrix_null(&pref_matrix_3u);
    MATRIX(pref_matrix_3u, 0, 1) = 1;
    MATRIX(pref_matrix_3u, 1, 0) = 1;
    MATRIX(pref_matrix_3u, 2, 2) = 1;

    igraph_matrix_init(&pref_matrix_nonsq, 3, 2);

    igraph_matrix_init(&pref_matrix_oor, 3, 3);
    igraph_matrix_null(&pref_matrix_oor);
    MATRIX(pref_matrix_oor, 0, 1) = 10;

    igraph_matrix_init(&pref_matrix_nsym, 3, 3);
    igraph_matrix_null(&pref_matrix_nsym);
    MATRIX(pref_matrix_nsym, 0, 1) = 1;

    igraph_vector_int_init_int(&block_sizes_0, 0);
    igraph_vector_int_init_int(&block_sizes_1, 1, 1);
    igraph_vector_int_init_int(&block_sizes_3, 3, 2, 2, 2);
    igraph_vector_int_init_int(&block_sizes_neg, 3, 2, 2, -2);

    printf("No vertices.\n");
    call_and_print(0, &pref_matrix_0, &block_sizes_0, 0, 0);

    printf("One vertex, directed, with loops.\n");
    call_and_print(1, &pref_matrix_1, &block_sizes_1, 1, 1);

    printf("Six vertices, directed, only edges from block 0 to 1 and 2 to 2.\n");
    call_and_print(6, &pref_matrix_3, &block_sizes_3, 1, 1);

    printf("Six vertices, directed, only edges from block 0 to 1 and 2 to 2, no loops.\n");
    call_and_print(6, &pref_matrix_3, &block_sizes_3, 1, 0);

    printf("Six vertices, undirected, only edges between block 0 and 1, and inside block 2.\n");
    call_and_print(6, &pref_matrix_3u, &block_sizes_3, 0, 1);

    printf("Six vertices, undirected, only edges between block 0 and 1, and inside block 2, no loops.\n");
    call_and_print(6, &pref_matrix_3u, &block_sizes_3, 0, 0);

    VERIFY_FINALLY_STACK();

    printf("Check for nonsquare matrix error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, 6, &pref_matrix_nonsq, &block_sizes_3, 0, 0), IGRAPH_NONSQUARE);

    printf("Check for preference matrix probability out of range error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, 6, &pref_matrix_oor, &block_sizes_3, 0, 0), IGRAPH_EINVAL);

    printf("Check for nonsymmetric preference matrix for undirected graph error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, 6, &pref_matrix_nsym, &block_sizes_3, 0, 0), IGRAPH_EINVAL);

    printf("Check for incorrect block size vector error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, 6, &pref_matrix_3, &block_sizes_1, 1, 0), IGRAPH_EINVAL);

    printf("Check for negative block size error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, 6, &pref_matrix_3, &block_sizes_neg, 1, 0), IGRAPH_EINVAL);

    printf("Check for sum of block sizes not equal to number of vertices error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, 3, &pref_matrix_3, &block_sizes_3, 1, 0), IGRAPH_EINVAL);

    igraph_matrix_destroy(&pref_matrix_0);
    igraph_matrix_destroy(&pref_matrix_1);
    igraph_matrix_destroy(&pref_matrix_3);
    igraph_matrix_destroy(&pref_matrix_3u);
    igraph_matrix_destroy(&pref_matrix_oor);
    igraph_matrix_destroy(&pref_matrix_nsym);
    igraph_matrix_destroy(&pref_matrix_nonsq);
    igraph_vector_int_destroy(&block_sizes_0);
    igraph_vector_int_destroy(&block_sizes_1);
    igraph_vector_int_destroy(&block_sizes_3);
    igraph_vector_int_destroy(&block_sizes_neg);

    VERIFY_FINALLY_STACK();
    return 0;
}
