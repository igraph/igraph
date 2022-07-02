/* IGraph library.
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

void call_and_print(igraph_t *graph, igraph_vector_t *weights, igraph_vector_int_list_t *neighborhoods) {
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    IGRAPH_ASSERT(igraph_local_scan_neighborhood_ecount(graph, &result, weights, neighborhoods) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&result);
    printf("\n");
}


int main() {
    igraph_t g_0, g_1, g_lmu, g_lm;
    igraph_vector_t weights, result;
    igraph_vector_int_list_t neighborhoods;
    igraph_vector_int_t n1;

    igraph_vector_init_real(&weights, 8, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
    igraph_vector_init(&result, 0);
    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1); //undirected

    igraph_vector_int_list_init(&neighborhoods, 0);


    printf("No vertices:\n");
    call_and_print(&g_0, NULL, &neighborhoods);

    printf("one vertex, empty neighborhood:\n");
    igraph_vector_int_init(&n1, 0);
    igraph_vector_int_list_push_back(&neighborhoods, &n1);
    call_and_print(&g_1, NULL, &neighborhoods);

    printf("Graph with unconnected vertices, loops and multiple edges, empty neighborhoods:\n");
    for (int i = 0; i < 5; i ++) {
        igraph_vector_int_init_int(&n1, 0);
        igraph_vector_int_list_push_back(&neighborhoods, &n1);
    }
    call_and_print(&g_lm, NULL, &neighborhoods);
    igraph_vector_int_list_clear(&neighborhoods);

    printf("Graph with unconnected vertices, loops and multiple edges, all the same neighborhood:\n");
    for (int i = 0; i < 6; i ++) {
        igraph_vector_int_init_int(&n1, 3, 0, 1, 2);
        igraph_vector_int_list_push_back(&neighborhoods, &n1);
    }
    call_and_print(&g_lm, NULL, &neighborhoods);
    igraph_vector_int_list_destroy(&neighborhoods);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lmu);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&result);
    igraph_vector_int_list_destroy(&neighborhoods);

    VERIFY_FINALLY_STACK();
    return 0;
}
