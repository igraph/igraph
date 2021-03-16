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

void call_and_print(igraph_t *graph) {
    igraph_vector_int_t result;
    igraph_vector_int_init(&result, 0);
    IGRAPH_ASSERT(igraph_list_triangles(graph, &result) == IGRAPH_SUCCESS);
    print_vector_int(&result);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) % 3 == 0);
    igraph_vector_int_destroy(&result);
    printf("\n");
}


int main() {
    igraph_t g_0, g_1, g_5_full, g_lm;

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_full(&g_5_full, 5, 0, IGRAPH_NO_LOOPS);
    igraph_small(&g_lm, 6, 1, 0,1, 0,1, 0,2, 0,2, 1,1, 1,2, 1,2, 1,3, 2,0, 2,3, 3,4, 3,4, -1);


    printf("No vertices:\n");
    call_and_print(&g_0);

    printf("One vertex:\n");
    call_and_print(&g_1);

    printf("Full graph of 5 vertices:\n");
    call_and_print(&g_5_full);

    printf("Graph with loops and multiple edges:\n");
    call_and_print(&g_lm);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_5_full);
    igraph_destroy(&g_lm);

    VERIFY_FINALLY_STACK();
    return 0;
}
