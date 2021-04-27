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

void chop_print_destroy(igraph_matrix_t *result) {
    matrix_chop(result, 1e-10);
    print_matrix(result);
    igraph_matrix_destroy(result);
}

int main() {
    igraph_t g;
    igraph_matrix_t result;
    igraph_vector_t roots, rootlevel;

    printf("Empty graph check:\n");
    igraph_small(&g, 0, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_ALL, /*roots*/ NULL, /*rootlevel*/ NULL) == IGRAPH_SUCCESS);
    chop_print_destroy(&result);
    igraph_destroy(&g);

    printf("Singleton graph check:\n");
    igraph_small(&g, 1, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_ALL, /*roots*/ NULL, /*rootlevel*/ NULL) == IGRAPH_SUCCESS);
    chop_print_destroy(&result);
    igraph_destroy(&g);

    printf("Star graph check with given root:\n");
    igraph_small(&g, 5, 1, 0,1, 0,2, 0,3, 0,4, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init_int(&roots, 1, 1);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_OUT, &roots, /*rootlevel*/ NULL) == IGRAPH_SUCCESS);
    chop_print_destroy(&result);
    igraph_destroy(&g);
    igraph_vector_destroy(&roots);

    printf("Star graph check with root found by topological sort:\n");
    igraph_small(&g, 5, 1, 1,0, 2,0, 3,0, 4,0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_IN, NULL, /*rootlevel*/ NULL) == IGRAPH_SUCCESS);
    chop_print_destroy(&result);
    igraph_destroy(&g);

    printf("Two minitrees without rootlevel:\n");
    igraph_small(&g, 6, 1, 0,1, 0,2, 3,4, 3,5, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init_int(&roots, 2, 0, 3);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_OUT, &roots, NULL) == IGRAPH_SUCCESS);
    chop_print_destroy(&result);
    igraph_destroy(&g);
    igraph_vector_destroy(&roots);

    printf("Two minitrees with rootlevel 10 and 20:\n");
    igraph_small(&g, 6, 1, 0,1, 0,2, 3,4, 3,5, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init_int(&roots, 2, 0, 3);
    igraph_vector_init_int(&rootlevel, 2, 10, 20);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_OUT, &roots, &rootlevel) == IGRAPH_SUCCESS);
    chop_print_destroy(&result);
    igraph_destroy(&g);
    igraph_vector_destroy(&roots);
    igraph_vector_destroy(&rootlevel);

    printf("Graph with just loops, triple edges and disconnected vertices:\n");
    igraph_small(&g, 5, 1, 0,0, 0,0, 0,0, 1,2, 1,2, 1,2, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_ALL, NULL, /*rootlevel*/ NULL) == IGRAPH_SUCCESS);
    chop_print_destroy(&result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Checking proper error handling:\n");
    printf("Giving negative root.\n");
    igraph_small(&g, 5, 1, 0,1, 0,2, 0,3, 0,4, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init_int(&roots, 1, -1);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_OUT, &roots, /*rootlevel*/ NULL) == IGRAPH_EINVVID);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);
    igraph_vector_destroy(&roots);

    printf("Giving negative rootlevel.\n");
    igraph_small(&g, 6, 1, 0,1, 0,2, 3,4, 3,5, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init_int(&roots, 2, 0, 3);
    igraph_vector_init_int(&rootlevel, 2, -10, -20);
    IGRAPH_ASSERT(igraph_layout_reingold_tilford_circular(&g, &result, IGRAPH_OUT, &roots, &rootlevel) == IGRAPH_EINVAL);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);
    igraph_vector_destroy(&roots);
    igraph_vector_destroy(&rootlevel);

    VERIFY_FINALLY_STACK();
    return 0;
}
