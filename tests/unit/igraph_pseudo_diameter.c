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

void print_and_destroy(igraph_t *g) {
    igraph_real_t result;
    IGRAPH_ASSERT(igraph_pseudo_diameter(g, &result) == IGRAPH_SUCCESS);
    printf("result:");
    print_real(stdout, result, "%g");
    printf("\n\n");
    igraph_destroy(g);
}

int main() {
    igraph_t graph;
    igraph_real_t result;
    int i;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No vertices:\n");
    igraph_small(&graph, 0, 0, -1);
    print_and_destroy(&graph);

    printf("1 vertex:\n");
    igraph_small(&graph, 1, 0, -1);
    print_and_destroy(&graph);

    printf("2 vertices:\n");
    igraph_small(&graph, 2, 0, -1);
    print_and_destroy(&graph);

    printf("Disconnected graph with loops and multiple edges.\n");
    igraph_small(&graph, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    for (i = 0; i < 10; i ++) {
        IGRAPH_ASSERT(igraph_pseudo_diameter(&graph, &result) == IGRAPH_SUCCESS);
        IGRAPH_ASSERT(result == 0 || result == 2 || result == 3);
    }
    igraph_destroy(&graph);

    printf("Tutte graph.\n");
    igraph_famous(&graph, "tutte");
    for (i = 0; i < 10; i ++) {
        IGRAPH_ASSERT(igraph_pseudo_diameter(&graph, &result) == IGRAPH_SUCCESS);
        IGRAPH_ASSERT(result == 8);
    }
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
    return 0;
}
