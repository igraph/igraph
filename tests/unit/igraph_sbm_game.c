/*
   igraph library.
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
#include "test_utilities.h"

void call_and_print(
        igraph_matrix_t *pref_matrix,
        igraph_vector_int_t *block_sizes,
        igraph_bool_t directed,
        igraph_bool_t loops, igraph_bool_t multiple) {

    igraph_t result;
    igraph_edge_type_sw_t allowed_edge_types = IGRAPH_SIMPLE_SW;
    if (loops) allowed_edge_types |= IGRAPH_LOOPS_SW;
    if (multiple) allowed_edge_types |= IGRAPH_MULTI_SW;

    IGRAPH_ASSERT(igraph_sbm_game(&result, pref_matrix, block_sizes, directed, allowed_edge_types) == IGRAPH_SUCCESS);
    print_graph_canon(&result);
    printf("\n");

    if (!loops) {
        igraph_bool_t has_loops;
        igraph_has_loop(&result, &has_loops);
        IGRAPH_ASSERT(! has_loops);
    }

    if (!multiple) {
        igraph_bool_t has_multi;
        igraph_has_multiple(&result, &has_multi);
        IGRAPH_ASSERT(! has_multi);
    }

    igraph_destroy(&result);
}


int main(void) {
    igraph_t result;
    igraph_matrix_t pref_matrix_0, pref_matrix_1, pref_matrix_2, pref_matrix_3, pref_matrix_3u, pref_matrix_nonsq, pref_matrix_oor, pref_matrix_nsym;
    igraph_vector_int_t block_sizes_0, block_sizes_1, block_sizes_2, block_sizes_3, block_sizes_neg;

    igraph_matrix_init(&pref_matrix_0, 0, 0);

    igraph_matrix_init(&pref_matrix_1, 1, 1);
    MATRIX(pref_matrix_1, 0, 0) = 1;

    igraph_matrix_init(&pref_matrix_2, 2, 2);

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
    igraph_vector_int_init_int(&block_sizes_2, 2, 1, 1);
    igraph_vector_int_init_int(&block_sizes_3, 3, 2, 2, 2);
    igraph_vector_int_init_int(&block_sizes_neg, 3, 2, 2, -2);

    /* Tests without multi-edges */

    printf("No vertices.\n");
    call_and_print(&pref_matrix_0, &block_sizes_0, false, false, false);

    printf("One vertex, directed, with loops.\n");
    call_and_print(&pref_matrix_1, &block_sizes_1, true, true, false);

    printf("Six vertices, directed, only edges from block 0 to 1 and 2 to 2.\n");
    call_and_print(&pref_matrix_3, &block_sizes_3, true, true, false);

    printf("Six vertices, directed, only edges from block 0 to 1 and 2 to 2, no loops.\n");
    call_and_print(&pref_matrix_3, &block_sizes_3, true, false, false);

    printf("Six vertices, undirected, only edges between block 0 and 1, and inside block 2.\n");
    call_and_print(&pref_matrix_3u, &block_sizes_3, false, true, false);

    printf("Six vertices, undirected, only edges between block 0 and 1, and inside block 2, no loops.\n");
    call_and_print(&pref_matrix_3u, &block_sizes_3, false, false, false);

    /* Smoke test for multi-edges */

    {
        const int trials = 100;
        igraph_real_t mean_ecount = 0;

        igraph_matrix_fill(&pref_matrix_2, 100);

        for (int i=0; i < trials; i++) {
            igraph_sbm_game(&result, &pref_matrix_2, &block_sizes_2, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW);
            mean_ecount += igraph_ecount(&result);
            igraph_destroy(&result);
        }
        mean_ecount /= trials;

        /* Since we disallowed loops, edges are present only between the two
         * distinct vertices of this graph. For 100 trials and a rate of p=100,
         * the probability that the edge count is below 80 is ~1.8%
         * and that it's above 120 is ~2.8%. */
        IGRAPH_ASSERT(80 < mean_ecount && mean_ecount < 120);
    }

    VERIFY_FINALLY_STACK();

    /* Verify validation */

    printf("Check for nonsquare matrix error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, &pref_matrix_nonsq, &block_sizes_3, false, IGRAPH_SIMPLE_SW), IGRAPH_EINVAL);

    printf("Check for preference matrix probability out of range error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, &pref_matrix_oor, &block_sizes_3, false, IGRAPH_SIMPLE_SW), IGRAPH_EINVAL);

    printf("Check for nonsymmetric preference matrix for undirected graph error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, &pref_matrix_nsym, &block_sizes_3, false, IGRAPH_SIMPLE_SW), IGRAPH_EINVAL);

    printf("Check for incorrect block size vector error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, &pref_matrix_3, &block_sizes_1, true, IGRAPH_SIMPLE_SW), IGRAPH_EINVAL);

    printf("Check for negative block size error handling.\n");
    CHECK_ERROR(igraph_sbm_game(&result, &pref_matrix_3, &block_sizes_neg, true, IGRAPH_SIMPLE_SW), IGRAPH_EINVAL);

    igraph_matrix_destroy(&pref_matrix_0);
    igraph_matrix_destroy(&pref_matrix_1);
    igraph_matrix_destroy(&pref_matrix_2);
    igraph_matrix_destroy(&pref_matrix_3);
    igraph_matrix_destroy(&pref_matrix_3u);
    igraph_matrix_destroy(&pref_matrix_oor);
    igraph_matrix_destroy(&pref_matrix_nsym);
    igraph_matrix_destroy(&pref_matrix_nonsq);
    igraph_vector_int_destroy(&block_sizes_0);
    igraph_vector_int_destroy(&block_sizes_1);
    igraph_vector_int_destroy(&block_sizes_2);
    igraph_vector_int_destroy(&block_sizes_3);
    igraph_vector_int_destroy(&block_sizes_neg);

    VERIFY_FINALLY_STACK();
    return 0;
}
