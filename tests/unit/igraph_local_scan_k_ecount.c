/* IGraph library.  Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

void call_and_print(igraph_t *graph, int k, igraph_vector_t *weights, igraph_neimode_t mode) {
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    IGRAPH_ASSERT(igraph_local_scan_k_ecount(graph, k, &result, weights, mode) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&result);
    printf("\n");
}


int main() {
    igraph_t g_0, g_1, g_lmu, g_lm, g_lm_nl;
    igraph_vector_t weights, result;

    igraph_vector_init_real(&weights, 8, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
    igraph_vector_init(&result, 0);
    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1); //undirected
    igraph_small(&g_lm_nl, 6, 1, 0,1, 0,2, 1,3, 2,0, 2,3, 3,4, 3,4, -1); // no loop

    printf("No vertices:\n");
    call_and_print(&g_0, 2, NULL, IGRAPH_ALL);

    printf("One vertex:\n");
    call_and_print(&g_1, 2, NULL, IGRAPH_ALL);

    printf("Directed disconnected graph with loops and multiple edges, no weights, k = 0, IGRAPH_IN:\n");
    call_and_print(&g_lm, 0, NULL, IGRAPH_IN);

    printf("Same graph, k=1:\n");
    call_and_print(&g_lm, 1, NULL, IGRAPH_IN);

    printf("Same graph, without loops, k=1:\n");
    call_and_print(&g_lm_nl, 1, NULL, IGRAPH_IN);

    printf("Same graph with loop, k=1, undirected:\n");
    call_and_print(&g_lmu, 1, NULL, IGRAPH_IN);

    printf("Checking if calling igraph_local_scan_1_ecount properly redirects:\n");
    igraph_vector_clear(&result);
    IGRAPH_ASSERT(igraph_local_scan_1_ecount(&g_lmu, &result, NULL, IGRAPH_IN) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("Same graph, directed, k=2:\n");
    call_and_print(&g_lm, 2, NULL, IGRAPH_IN);

    printf("Same graph, undirected, k=2:\n");
    call_and_print(&g_lmu, 2, NULL, IGRAPH_IN);

    printf("Same graph, weighted:\n");
    call_and_print(&g_lmu, 2, &weights, IGRAPH_IN);

    VERIFY_FINALLY_STACK();

    printf("Wrong size weights.\n");
    igraph_vector_clear(&weights);
    CHECK_ERROR(igraph_local_scan_k_ecount(&g_lmu, 3, &result, &weights, IGRAPH_ALL), IGRAPH_EINVAL);

    printf("Negative k.\n");
    CHECK_ERROR(igraph_local_scan_k_ecount(&g_lmu, -3, &result, NULL, IGRAPH_ALL), IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_lmu);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lm_nl);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
