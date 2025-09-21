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

void call_and_print(igraph_t *us, igraph_t *them, int k, igraph_vector_t *weights, igraph_neimode_t mode) {
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    IGRAPH_ASSERT(igraph_local_scan_k_ecount_them(us, them, k, &result, weights, mode) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&result);
    printf("\n");
}


int main(void) {
    igraph_t g_0, g_1, g_lmu, g_lm, g_lm_nl, g_6, g_6_1, g_6_full;
    igraph_vector_t weights, result;

    igraph_vector_init_real(&weights, 8, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
    igraph_vector_init(&result, 0);
    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_6, 6, 1, -1);
    igraph_small(&g_6_1, 6, 1, 0,1, -1);
    igraph_full(&g_6_full, 6, 1, 0);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1); //undirected
    igraph_small(&g_lm_nl, 6, 1, 0,1, 0,2, 1,3, 2,0, 2,3, 3,4, 3,4, -1); // no loop

    printf("First some tests where us and them are equal:\n");
    printf("No vertices:\n");
    call_and_print(&g_0, &g_0, 2, NULL, IGRAPH_ALL);

    printf("One vertex:\n");
    call_and_print(&g_1, &g_1, 2, NULL, IGRAPH_ALL);

    printf("Directed disconnected graph with loops and multiple edges, no weights, k = 0, IGRAPH_IN:\n");
    call_and_print(&g_lm, &g_lm, 0, NULL, IGRAPH_IN);

    printf("Same graph, with weights:\n");
    call_and_print(&g_lm, &g_lm, 0, &weights, IGRAPH_IN);

    printf("Same graph, k=1:\n");
    call_and_print(&g_lm, &g_lm, 1, NULL, IGRAPH_IN);

    printf("Same graph, k=1, IGRAPH_ALL:\n");
    call_and_print(&g_lm, &g_lm, 1, NULL, IGRAPH_ALL);

    printf("Same graph, without loops, k=1:\n");
    call_and_print(&g_lm_nl, &g_lm_nl, 1, NULL, IGRAPH_IN);

    printf("Same graph with loop, k=1, undirected:\n");
    call_and_print(&g_lmu, &g_lmu, 1, NULL, IGRAPH_IN);

    printf("Same graph, directed, k=2:\n");
    call_and_print(&g_lm, &g_lm, 2, NULL, IGRAPH_IN);

    printf("Same graph, undirected, k=2:\n");
    call_and_print(&g_lmu, &g_lmu, 2, NULL, IGRAPH_IN);

    printf("Same graph, weighted:\n");
    call_and_print(&g_lmu, &g_lmu, 2, &weights, IGRAPH_IN);

    printf("Now some tests where us and them are not equal:\n");
    printf("Us = same graph, them = edgless, directed, k=1, "
           "should show 0, because there are no edges:\n");
    call_and_print(&g_lm, &g_6, 1, NULL, IGRAPH_IN);

    printf("Switched us and them, "
           "should show only 1 at the loop:\n");
    call_and_print(&g_6, &g_lm, 1, NULL, IGRAPH_IN);

    printf("Us = same graph, them = only edge from 0 to 1, "
           "directed, k=1, IGRAPH_IN, "
           "should show edge from 0 to 1:\n");
    call_and_print(&g_lm, &g_6_1, 1, NULL, IGRAPH_IN);

    printf("Switched us and them, "
           "should show edge from 0 to 1 and loop:\n");
    call_and_print(&g_6_1, &g_lm, 1, NULL, IGRAPH_IN);

    printf("Us = same graph, them = full graph, "
           "directed, k=3, IGRAPH_ALL, "
           "should show 4*5=20 edges for the connected part:\n");
    call_and_print(&g_lm, &g_6_full, 3, NULL, IGRAPH_ALL);

    printf("Switched us and them, "
           "should show 8 edges for the whole graph:\n");
    call_and_print(&g_6_full, &g_lm, 3, NULL, IGRAPH_ALL);
    VERIFY_FINALLY_STACK();

    printf("Check error handling:\n");
    printf("Wrong size weights.\n");
    igraph_vector_clear(&weights);
    CHECK_ERROR(igraph_local_scan_k_ecount_them(&g_lmu, &g_lmu, 3, &result, &weights, IGRAPH_ALL), IGRAPH_EINVAL);

    printf("Negative k.\n");
    CHECK_ERROR(igraph_local_scan_k_ecount_them(&g_lmu, &g_lmu, -3, &result, NULL, IGRAPH_ALL), IGRAPH_EINVAL);

    printf("Number of vertices in us and them not equal.\n");
    CHECK_ERROR(igraph_local_scan_k_ecount_them(&g_lmu, &g_1, 3, &result, NULL, IGRAPH_ALL), IGRAPH_EINVAL);

    printf("Directedness in us and them not equal.\n");
    CHECK_ERROR(igraph_local_scan_k_ecount_them(&g_lmu, &g_lm, 3, &result, NULL, IGRAPH_ALL), IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_lmu);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lm_nl);
    igraph_destroy(&g_6);
    igraph_destroy(&g_6_1);
    igraph_destroy(&g_6_full);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
