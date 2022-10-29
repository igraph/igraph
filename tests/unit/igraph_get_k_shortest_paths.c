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
#include "test_utilities.h"

void call_and_print(
    igraph_t *graph, igraph_vector_t *weights, igraph_integer_t k,
    igraph_integer_t from, igraph_integer_t to, igraph_neimode_t mode
) {
    igraph_vector_int_list_t vertex_paths;
    igraph_vector_int_list_t edge_paths;
    igraph_vector_int_list_t paths_for_verification;
    igraph_integer_t i, n;

    igraph_vector_int_list_init(&vertex_paths, 0);
    igraph_vector_int_list_init(&edge_paths, 0);
    igraph_vector_int_list_init(&paths_for_verification, 0);

    IGRAPH_ASSERT(igraph_get_k_shortest_paths(graph, weights, &vertex_paths, &edge_paths, k, from, to, mode) == IGRAPH_SUCCESS);

    printf("result (vertex IDs): \n");
    print_vector_int_list(&vertex_paths);

    printf("result (edge IDs): \n");
    print_vector_int_list(&edge_paths);

    /* now we execute the same test but with only one of vertex_paths or edge_paths,
     * and we check that we get the same result */
    IGRAPH_ASSERT(igraph_get_k_shortest_paths(graph, weights, NULL, &paths_for_verification, k, from, to, mode) == IGRAPH_SUCCESS);
    n = igraph_vector_int_list_size(&paths_for_verification);
    IGRAPH_ASSERT(n == igraph_vector_int_list_size(&edge_paths));
    for (i = 0; i < n; i++) {
        IGRAPH_ASSERT(igraph_vector_int_is_equal(
            igraph_vector_int_list_get_ptr(&edge_paths, i),
            igraph_vector_int_list_get_ptr(&paths_for_verification, i)
        ));
    }

    IGRAPH_ASSERT(igraph_get_k_shortest_paths(graph, weights, &paths_for_verification, NULL, k, from, to, mode) == IGRAPH_SUCCESS);
    n = igraph_vector_int_list_size(&paths_for_verification);
    IGRAPH_ASSERT(n == igraph_vector_int_list_size(&vertex_paths));
    for (i = 0; i < n; i++) {
        IGRAPH_ASSERT(igraph_vector_int_is_equal(
            igraph_vector_int_list_get_ptr(&vertex_paths, i),
            igraph_vector_int_list_get_ptr(&paths_for_verification, i)
        ));
    }

    igraph_vector_int_list_destroy(&paths_for_verification);
    igraph_vector_int_list_destroy(&vertex_paths);
    igraph_vector_int_list_destroy(&edge_paths);
    printf("\n");

    VERIFY_FINALLY_STACK();
}


int main(void) {
    /* Wiki example taken from https://en.wikipedia.org/wiki/Yen's_algorithm */
    igraph_t g_0, g_1, g_2, g_2c, g_wiki, g_wiki_u, g_sz, g_lattice;
    igraph_vector_t weights, weights_wiki, weights_inf;
    igraph_vector_int_t dims;
    igraph_vector_int_list_t paths;

    igraph_vector_int_list_init(&paths, 0);

    igraph_small(&g_0, 0, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_1, 1, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_2, 2, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_2c, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    igraph_small(&g_wiki, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,3, 2,1, 2,3, 2,4, 3,4, 3,5, 4,5, -1);
    igraph_small(&g_wiki_u, 6, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,3, 2,1, 2,3, 2,4, 3,4, 3,5, 4,5, -1);
    igraph_small(&g_sz, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 1, 2, 2, 0, 0, 3, 3, 4, 4, 2, 2, 2, 2, 4, 0, 4, -1);

    igraph_vector_int_init_int(&dims, 2, 3, 3);
    igraph_square_lattice(&g_lattice, &dims, 1, IGRAPH_UNDIRECTED, /* mutual = */ 0, /* periodic = */ NULL);
    igraph_vector_int_destroy(&dims);

    igraph_vector_init(&weights, 0);
    igraph_vector_init_int(&weights_wiki, 9, 3, 2, 4, 1, 2, 3, 2, 1, 2);
    igraph_vector_init_real(&weights_inf, 1, IGRAPH_INFINITY);

    printf("One vertex:\n");
    call_and_print(&g_1, &weights, 2, 0, 0, IGRAPH_ALL);

    printf("Two disconnected vertices:\n");
    call_and_print(&g_2, &weights, 2, 0, 1, IGRAPH_ALL);

    printf("Two vertices with infinite weight edge in between:\n");
    call_and_print(&g_2c, &weights_inf, 2, 0, 1, IGRAPH_ALL);

    printf("Wiki example:\n");
    call_and_print(&g_wiki, &weights_wiki, 10, 0, 5, IGRAPH_OUT);

    printf("Wiki example, 0 shortest paths:\n");
    call_and_print(&g_wiki, &weights_wiki, 0, 0, 5, IGRAPH_OUT);

    printf("Wiki example, 2 shortest paths:\n");
    call_and_print(&g_wiki, &weights_wiki, 2, 0, 5, IGRAPH_OUT);

    printf("Wiki example, other direction:\n");
    call_and_print(&g_wiki, &weights_wiki, 10, 5, 0, IGRAPH_IN);

    printf("Wiki example, direction ignored:\n");
    call_and_print(&g_wiki, &weights_wiki, 20, 5, 0, IGRAPH_ALL);

    printf("Wiki example, undirected:\n");
    call_and_print(&g_wiki_u, &weights_wiki, 20, 5, 0, IGRAPH_ALL);

    printf("Wiki example, no weights:\n");
    call_and_print(&g_wiki, NULL, 10, 0, 5, IGRAPH_OUT);

    printf("Directed unweighted graph:\n");
    call_and_print(&g_sz, NULL, 4, 0, 4, IGRAPH_OUT);

    printf("3x3 square lattice:\n");
    call_and_print(&g_lattice, NULL, 6, 0, 8, IGRAPH_OUT);

    printf("Zero vertices, from and to don't exist:\n");
    CHECK_ERROR(igraph_get_k_shortest_paths(&g_0, &weights, NULL, &paths, 4, 0, 0, IGRAPH_ALL), IGRAPH_EINVVID);

    printf("Wrong weights length:\n");
    CHECK_ERROR(igraph_get_k_shortest_paths(&g_wiki, &weights, NULL, &paths, 4, 0, 5, IGRAPH_ALL), IGRAPH_EINVAL);

    printf("Non-existent mode:\n");
    CHECK_ERROR(igraph_get_k_shortest_paths(&g_1, &weights, NULL, &paths, 4, 0, 0, (igraph_neimode_t) 100), IGRAPH_EINVMODE);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_2);
    igraph_destroy(&g_2c);
    igraph_destroy(&g_wiki);
    igraph_destroy(&g_wiki_u);
    igraph_destroy(&g_sz);
    igraph_destroy(&g_lattice);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&weights_wiki);
    igraph_vector_destroy(&weights_inf);
    igraph_vector_int_list_destroy(&paths);

    VERIFY_FINALLY_STACK();
    return 0;
}
