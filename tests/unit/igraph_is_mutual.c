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

void call_and_print(igraph_t *graph, igraph_es_t es) {
    igraph_vector_bool_t result;
    igraph_vector_bool_init(&result, 0);
    IGRAPH_ASSERT(igraph_is_mutual(graph, &result, es) == IGRAPH_SUCCESS);
    igraph_vector_bool_print(&result);
    printf("\n");
    igraph_vector_bool_destroy(&result);
}


int main() {
    igraph_t g_0, g_lm, g_lmu;
    igraph_vector_bool_t result;

    igraph_vector_bool_init(&result, 0);

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices:\n");
    call_and_print(&g_0, igraph_ess_all(IGRAPH_EDGEORDER_ID));

    printf("Graph with loops and multiple edges:\n");
    call_and_print(&g_lm, igraph_ess_all(IGRAPH_EDGEORDER_ID));

    printf("Same graph, selecting edge 4:\n");
    call_and_print(&g_lm, igraph_ess_1(4));

    printf("Same graph, but undirected:\n");
    call_and_print(&g_lmu, igraph_ess_all(IGRAPH_EDGEORDER_ID));

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Edge out of range.\n");
    IGRAPH_ASSERT(igraph_is_mutual(&g_lm, &result, igraph_ess_1(100)) == IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lmu);
    igraph_vector_bool_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
