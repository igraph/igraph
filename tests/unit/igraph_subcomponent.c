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

void call_and_print(igraph_t *graph, igraph_real_t vertex, igraph_neimode_t mode) {
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    IGRAPH_ASSERT(igraph_subcomponent(graph, &result, vertex, mode) == IGRAPH_SUCCESS);
    igraph_vector_sort(&result);
    igraph_vector_print(&result);
    igraph_vector_destroy(&result);
    printf("\n");
}


int main() {
    igraph_t g_0, g_1, g_lm, g_lmu;
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    int i;

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices, should give error for impossible starting vertex.\n");
    CHECK_ERROR(igraph_subcomponent(&g_0, &result, 0, IGRAPH_ALL), IGRAPH_EINVVID);

    printf("One vertex.\n");
    call_and_print(&g_1, 0, IGRAPH_ALL);

    printf("All vertices of a graph, IGRAPH_OUT:\n");
    for (i = 0; i < 6; i++) {
        call_and_print(&g_lm, i, IGRAPH_OUT);
    }

    printf("All vertices of a graph, IGRAPH_IN:\n");
    for (i = 0; i < 6; i++) {
        call_and_print(&g_lm, i, IGRAPH_IN);
    }

    printf("All vertices of a graph, IGRAPH_ALL:\n");
    for (i = 0; i < 6; i++) {
        call_and_print(&g_lm, i, IGRAPH_ALL);
    }

    printf("All vertices of a graph, undirected:\n");
    for (i = 0; i < 6; i++) {
        call_and_print(&g_lmu, i, IGRAPH_OUT);
    }

    printf("Check for invalid mode error handling.\n");
    CHECK_ERROR(igraph_subcomponent(&g_1, &result, 0, 100), IGRAPH_EINVMODE);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lmu);
    igraph_vector_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
