/* igraph library.
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

void call_and_print(igraph_t *graph, igraph_vector_t *weights, igraph_vector_int_list_t *subsets) {
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    IGRAPH_ASSERT(igraph_local_scan_subset_ecount(graph, &result, weights, subsets) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&result);
    printf("\n");
}


int main(void) {
    igraph_t g_0, g_1, g_lmu, g_lm;
    igraph_vector_t weights, result;
    igraph_vector_int_list_t subsets;
    igraph_vector_int_t n1;

    igraph_vector_init_real(&weights, 8, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
    igraph_vector_init(&result, 0);
    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1); //undirected

    igraph_vector_int_list_init(&subsets, 0);


    printf("No vertices:\n");
    call_and_print(&g_0, NULL, &subsets);

    printf("one vertex, empty subset:\n");
    igraph_vector_int_init(&n1, 0);
    igraph_vector_int_list_push_back(&subsets, &n1);
    call_and_print(&g_1, NULL, &subsets);

    printf("Graph with unconnected vertices, loops and multiple edges, empty subsets:\n");
    for (int i = 0; i < 5; i ++) {
        igraph_vector_int_init_int(&n1, 0);
        igraph_vector_int_list_push_back(&subsets, &n1);
    }
    call_and_print(&g_lm, NULL, &subsets);
    igraph_vector_int_list_clear(&subsets);

    printf("Same graph, every vertex as a subset:\n");
    for (int i = 0; i < 6; i ++) {
        igraph_vector_int_init_int(&n1, 1, i);
        igraph_vector_int_list_push_back(&subsets, &n1);
    }
    call_and_print(&g_lm, NULL, &subsets);
    igraph_vector_int_list_clear(&subsets);

    printf("Same graph, every vertex and next one as a subset:\n");
    for (int i = 0; i < 6; i ++) {
        igraph_vector_int_init_int(&n1, 2, i, (i + 1) % 6);
        igraph_vector_int_list_push_back(&subsets, &n1);
    }
    call_and_print(&g_lm, NULL, &subsets);

    printf("Same situation, but with weights:\n");
    call_and_print(&g_lm, &weights, &subsets);

    printf("Same situation, no weights, but undirected graph:\n");
    call_and_print(&g_lmu, NULL, &subsets);
    igraph_vector_int_list_clear(&subsets);

    printf("Same graph, whole graph as subset:\n");
    igraph_vector_int_init_int(&n1, 6, 0, 1, 2, 3, 4, 5);
    igraph_vector_int_list_push_back(&subsets, &n1);
    call_and_print(&g_lmu, NULL, &subsets);
    igraph_vector_int_list_destroy(&subsets);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lmu);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&result);
    igraph_vector_int_list_destroy(&subsets);

    VERIFY_FINALLY_STACK();
    return 0;
}
