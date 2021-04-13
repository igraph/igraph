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

void call_and_print(igraph_t *graph, igraph_vector_t *weights, igraph_integer_t k, igraph_integer_t from, igraph_integer_t to, igraph_neimode_t mode) {
    igraph_vector_ptr_t paths;
    igraph_vector_ptr_init(&paths, 0);
    IGRAPH_ASSERT(igraph_k_shortest_paths(graph, weights, &paths, k, from, to, mode) == IGRAPH_SUCCESS);
    printf("result: \n");
    for (int i = 0; i < igraph_vector_ptr_size(&paths); i++) {
        print_vector(VECTOR(paths)[i]);
    }
    igraph_vector_ptr_destroy_all(&paths);
    printf("\n");
}


int main() {
    igraph_t g_0, g_1, g_2, g_wiki;
    igraph_vector_t weights, weights_wiki;
    igraph_vector_ptr_t paths;

    igraph_vector_ptr_init(&paths, 0);

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_2, 2, 0, -1);
    igraph_small(&g_wiki, 6, 1, 0,1, 0,2, 1,3, 2,1, 2,3, 2,4, 3,4, 3,5, 4,5, -1);

    igraph_vector_init(&weights, 0);
    igraph_vector_init_int(&weights_wiki, 9, 3, 2, 4, 1, 2, 3, 2, 1, 2);

    printf("One vertex:\n");
    call_and_print(&g_1, &weights, 2, 0, 0, IGRAPH_ALL);

    printf("Two disconnected vertices:\n");
    call_and_print(&g_2, &weights, 2, 0, 1, IGRAPH_ALL);

    printf("wiki example:\n");
    call_and_print(&g_wiki, &weights_wiki, 10, 0, 5, IGRAPH_OUT);

    printf("wiki example, 0 shortest paths:\n");
    call_and_print(&g_wiki, &weights_wiki, 0, 0, 5, IGRAPH_OUT);

    //printf("wiki example, 2 shortest paths:\n"); //memory leak
    //call_and_print(&g_wiki, &weights_wiki, 2, 0, 5, IGRAPH_OUT);

    printf("wiki example, other direction:\n");
    call_and_print(&g_wiki, &weights_wiki, 10, 5, 0, IGRAPH_IN);

    printf("wiki example, direction ignored:\n");
    call_and_print(&g_wiki, &weights_wiki, 20, 5, 0, IGRAPH_ALL);

    printf("wiki example, no weights:\n");
    call_and_print(&g_wiki, NULL, 10, 0, 5, IGRAPH_OUT);

    igraph_set_error_handler(igraph_error_handler_printignore);

    printf("Zero vertices, from and to don't exist:\n");
    IGRAPH_ASSERT(igraph_k_shortest_paths(&g_0, &weights, &paths, 4, 0, 0, IGRAPH_ALL) == IGRAPH_EINVAL);

    printf("Wrong weights length:\n");
    IGRAPH_ASSERT(igraph_k_shortest_paths(&g_wiki, &weights, &paths, 4, 0, 5, IGRAPH_ALL) == IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_2);
    igraph_destroy(&g_wiki);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&weights_wiki);
    igraph_vector_ptr_destroy_all(&paths);

    VERIFY_FINALLY_STACK();
    return 0;
}
