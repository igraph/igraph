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

int main() {
    igraph_t g_empty, g_empty_dir, g_lm;
    igraph_matrix_t result;
    igraph_vs_t vids;
    igraph_vector_t weights_empty, weights_lm, weights_lm_neg_loop;

    igraph_matrix_init(&result, 0, 0);
    igraph_vs_all(&vids);
    igraph_vector_init(&weights_empty, 0);
    igraph_vector_init_int(&weights_lm_neg_loop, 9, -4, -3, -2, -1, 0, 1, 2, 3, 4);
    igraph_vector_init_int(&weights_lm, 9, -1, 0, 1, -2, 2, 3, 4, 5, 6);
    igraph_small(&g_empty, 0, 0, -1);
    igraph_small(&g_empty_dir, 0, 1, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,2, 1,3, 2,0, 2,3, 3,4, 3,4, -1);

    igraph_set_error_handler(igraph_error_handler_printignore);

    printf("No vertices, not directed:\n");
    IGRAPH_ASSERT(igraph_shortest_paths_johnson(&g_empty, &result, vids, vids, &weights_empty) == IGRAPH_SUCCESS);
    print_matrix(&result);

    printf("No vertices, directed:\n");
    IGRAPH_ASSERT(igraph_shortest_paths_johnson(&g_empty_dir, &result, vids, vids, &weights_empty) == IGRAPH_SUCCESS);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges with negative loop:\n");
    IGRAPH_ASSERT(igraph_shortest_paths_johnson(&g_lm, &result, vids, vids, &weights_lm_neg_loop) == IGRAPH_ENEGLOOP);

    printf("Directed graph with loops and multi-edges:\n");
    IGRAPH_ASSERT(igraph_shortest_paths_johnson(&g_lm, &result, vids, vids, &weights_lm) == IGRAPH_SUCCESS);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges, select vertices 1 and 2:\n");
    IGRAPH_ASSERT(igraph_shortest_paths_johnson(&g_lm, &result, igraph_vss_seq(1, 2), igraph_vss_seq(1, 2), &weights_lm) == IGRAPH_SUCCESS);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges, select 0 -> 2:\n");
    IGRAPH_ASSERT(igraph_shortest_paths_johnson(&g_lm, &result, igraph_vss_1(0), igraph_vss_1(2), &weights_lm) == IGRAPH_SUCCESS);
    print_matrix(&result);

    printf("Directed graph with loops and multi-edges, select none:\n");
    IGRAPH_ASSERT(igraph_shortest_paths_johnson(&g_lm, &result, igraph_vss_none(), igraph_vss_none(), &weights_lm) == IGRAPH_SUCCESS);
    print_matrix(&result);

    igraph_matrix_destroy(&result);
    igraph_destroy(&g_empty);
    igraph_destroy(&g_empty_dir);
    igraph_destroy(&g_lm);
    igraph_vector_destroy(&weights_empty);
    igraph_vector_destroy(&weights_lm);
    igraph_vector_destroy(&weights_lm_neg_loop);

    VERIFY_FINALLY_STACK();
    return 0;
}
