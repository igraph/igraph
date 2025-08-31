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

int main(void) {
    igraph_t g_empty, g_lm;
    igraph_matrix_t result;
    igraph_matrix_t bf_result;
    igraph_vs_t vids;
    igraph_vector_t weights_empty, weights_lm, weights_lm_neg_loop;

    igraph_matrix_init(&result, 0, 0);
    igraph_matrix_init(&bf_result, 0, 0);
    igraph_vs_all(&vids);
    igraph_vector_init(&weights_empty, 0);
    igraph_vector_init_int(&weights_lm_neg_loop, 9, -4, -3, -2, -1, 0, 1, 2, 3, 4);
    igraph_vector_init_int(&weights_lm, 9, -1, 0, 1, -2, 2, 3, 4, 5, 6);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,2, 1,3, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices, undirected:\n");
    igraph_empty(&g_empty, 0, IGRAPH_UNDIRECTED);
    igraph_vector_resize(&weights_empty, igraph_ecount(&g_empty));
    igraph_distances_johnson(&g_empty, &result, vids, vids, &weights_empty, IGRAPH_OUT);
    igraph_destroy(&g_empty);
    print_matrix(&result);

    printf("No vertices, directed:\n");
    igraph_empty(&g_empty, 0, IGRAPH_DIRECTED);
    igraph_vector_resize(&weights_empty, igraph_ecount(&g_empty));
    igraph_distances_johnson(&g_empty, &result, vids, vids, &weights_empty, IGRAPH_OUT);
    igraph_destroy(&g_empty);
    print_matrix(&result);

    printf("No edges, undirected:\n");
    igraph_empty(&g_empty, 3, IGRAPH_UNDIRECTED);
    igraph_vector_resize(&weights_empty, igraph_ecount(&g_empty));
    igraph_distances_johnson(&g_empty, &result, vids, vids, &weights_empty, IGRAPH_OUT);
    igraph_destroy(&g_empty);
    print_matrix(&result);

    printf("No edges, directed:\n");
    igraph_empty(&g_empty, 4, IGRAPH_DIRECTED);
    igraph_vector_resize(&weights_empty, igraph_ecount(&g_empty));
    igraph_distances_johnson(&g_empty, &result, vids, vids, &weights_empty, IGRAPH_OUT);
    igraph_destroy(&g_empty);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges:\n");
    igraph_distances_johnson(&g_lm, &result, vids, vids, &weights_lm, IGRAPH_OUT);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges, select vertices 1 and 2:\n");
    igraph_distances_johnson(&g_lm, &result, igraph_vss_range(1, 3), igraph_vss_range(1, 3), &weights_lm, IGRAPH_OUT);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges, select 0 -> 2:\n");
    igraph_distances_johnson(&g_lm, &result, igraph_vss_1(0), igraph_vss_1(2), &weights_lm, IGRAPH_OUT);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges, select none:\n");
    igraph_distances_johnson(&g_lm, &result, igraph_vss_none(), igraph_vss_none(), &weights_lm, IGRAPH_OUT);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges, IGRAPH_IN:\n");
    igraph_distances_johnson(&g_lm, &result, vids, vids, &weights_lm, IGRAPH_IN);
    igraph_distances_bellman_ford(&g_lm, &bf_result, vids, vids, &weights_lm, IGRAPH_IN);
    print_matrix(&result);
    IGRAPH_ASSERT(igraph_matrix_all_e(&result, &bf_result));

    VERIFY_FINALLY_STACK();

    printf("Checking error for directed graph with loops and multi-edges with negative loop.\n");
    CHECK_ERROR(igraph_distances_johnson(&g_lm, &result, vids, vids, &weights_lm_neg_loop, IGRAPH_OUT), IGRAPH_ENEGCYCLE);

    printf("Directed graph with loops and multi-edges, IGRAPH_ALL:\n");
    CHECK_ERROR(igraph_distances_johnson(&g_lm, &result, vids, vids, &weights_lm, IGRAPH_ALL), IGRAPH_ENEGCYCLE);

    igraph_matrix_destroy(&result);
    igraph_matrix_destroy(&bf_result);
    igraph_destroy(&g_lm);
    igraph_vector_destroy(&weights_empty);
    igraph_vector_destroy(&weights_lm);
    igraph_vector_destroy(&weights_lm_neg_loop);

    VERIFY_FINALLY_STACK();
    return 0;
}
