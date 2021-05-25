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

void compute_and_print(igraph_t *graph, igraph_vector_t *weights, igraph_bool_t directed, igraph_bool_t unconn) {
    igraph_real_t result;
    igraph_real_t unconn_pairs;

    IGRAPH_ASSERT(igraph_average_path_length_dijkstra(graph, &result, &unconn_pairs, weights,
                  directed, unconn) == IGRAPH_SUCCESS);

    printf("Result: ");
    print_real(stdout, result, "%8g");
    printf("\nUnconnected pairs: ");
    print_real(stdout, unconn_pairs, "%8g");
    printf("\n\n");
}

int main() {
    igraph_t g_0, g_1, g_2, g_3, g_lm;
    igraph_vector_t weights_0, weights_3, weights_lm, weights_lm_neg;
    igraph_real_t result;

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_2, 2, 0, -1);
    igraph_small(&g_3, 2, 1, 0,1, 0,2, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);

    igraph_vector_init(&weights_0, 0);
    igraph_vector_init_int(&weights_3, 2, 1, 1);
    igraph_vector_init_int(&weights_lm, 8, 0, 1, 2, 3, 4, 5, 6, 7);
    igraph_vector_init_int(&weights_lm_neg, 8, -10, 1, 2, 3, 4, 5, 6, 7);

    printf("No vertices:\n");
    compute_and_print(&g_0, &weights_0, 1, 1);

    printf("One vertex:\n");
    compute_and_print(&g_1, &weights_0, 1, 1);

    printf("Two vertices:\n");
    compute_and_print(&g_2, &weights_0, 1, 1);

    printf("Two vertices, inf for unconnected pairs:\n");
    compute_and_print(&g_2, &weights_0, 1, 0);

    printf("Smallest bifurcating directed tree:\n");
    compute_and_print(&g_3, &weights_3, 1, 1);

    printf("Smallest bifurcating directed tree, inf for unconnected pairs:\n");
    compute_and_print(&g_3, &weights_3, 1, 0);

    printf("Graph with loops and multiple edges:\n");
    compute_and_print(&g_lm, &weights_lm, 1, 1);

    printf("Graph with loops and multiple edges, ignoring direction:\n");
    compute_and_print(&g_lm, &weights_lm, 0, 1);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Checking incorrect weight length error handling.\n");
    IGRAPH_ASSERT(igraph_average_path_length_dijkstra(&g_lm, &result, NULL, &weights_0,
                  1, 1) == IGRAPH_EINVAL);

    printf("Checking negative weight error handling.\n");
    IGRAPH_ASSERT(igraph_average_path_length_dijkstra(&g_lm, &result, NULL, &weights_lm_neg,
                  1, 1) == IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_2);
    igraph_destroy(&g_3);
    igraph_destroy(&g_lm);
    igraph_vector_destroy(&weights_0);
    igraph_vector_destroy(&weights_3);
    igraph_vector_destroy(&weights_lm);
    igraph_vector_destroy(&weights_lm_neg);

    VERIFY_FINALLY_STACK();
    return 0;
}
