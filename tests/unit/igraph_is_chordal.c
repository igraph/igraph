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

void call_and_print(igraph_t *graph, igraph_vector_t *alpha, igraph_vector_t *alpham1,
                    igraph_bool_t fill, igraph_bool_t ng) {
    igraph_bool_t chordal;
    igraph_vector_t fill_in;
    igraph_t newgraph;
    igraph_vector_init(&fill_in, 0);
    IGRAPH_ASSERT(igraph_is_chordal(graph, alpha, alpham1, &chordal,
                  fill ? &fill_in : NULL, ng? &newgraph : NULL) == IGRAPH_SUCCESS);
    printf("Is chordal: %d\nFill in:\n", chordal);
    print_vector(&fill_in);
    if (ng) {
        printf("New graph:\n");
        print_graph_canon(&newgraph);
        igraph_destroy(&newgraph);
    }
    printf("\n");
    igraph_vector_destroy(&fill_in);
}


int main() {
    igraph_t g_0, g_1, g_lmu;
    igraph_bool_t chordal;
    igraph_vector_t alpha, alpham1;

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices:\n");
    call_and_print(&g_0, NULL, NULL, 1, 1);

    printf("One vertex:\n");
    call_and_print(&g_1, NULL, NULL, 1, 1);

    printf("One vertex, don't calculate anything.\n\n");
    IGRAPH_ASSERT(igraph_is_chordal(&g_1, NULL, NULL, NULL, NULL, NULL) == IGRAPH_SUCCESS);

    printf("Disconnected graph with loops and multiple edges:\n");
    call_and_print(&g_lmu, NULL, NULL, 1, 1);

    printf("Same graph, don't ask for fill_in vector:\n");
    call_and_print(&g_lmu, NULL, NULL, 0, 1);

    printf("Same graph, don't ask for fill_in vector or newgraph:\n");
    call_and_print(&g_lmu, NULL, NULL, 0, 0);

    printf("Same graph, own calculation of alpha and its inverse:\n");
    igraph_vector_init(&alpha, 0);
    igraph_vector_init(&alpham1, 0);
    igraph_maximum_cardinality_search(&g_lmu, &alpha, &alpham1);
    call_and_print(&g_lmu, &alpha, &alpham1, 1, 1);

    printf("Same graph, own calculation of alpha:\n");
    call_and_print(&g_lmu, &alpha, NULL, 1, 1);

    printf("Same graph, own calculation of inverse alpha:\n");
    call_and_print(&g_lmu, NULL, &alpham1, 1, 1);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Wrong size alpha.\n");
    igraph_vector_clear(&alpha);
    IGRAPH_ASSERT(igraph_is_chordal(&g_lmu, &alpha, NULL, &chordal, NULL, NULL) == IGRAPH_EINVAL);

    printf("Wrong size alpham1.\n");
    IGRAPH_ASSERT(igraph_is_chordal(&g_lmu, NULL, &alpha, &chordal, NULL, NULL) == IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_lmu);
    igraph_vector_destroy(&alpha);
    igraph_vector_destroy(&alpham1);

    VERIFY_FINALLY_STACK();
    return 0;
}
