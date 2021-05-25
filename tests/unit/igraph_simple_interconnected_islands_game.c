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

int main() {
    igraph_t g;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No islands:\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/0, /*size of islands*/ 1,
                  /*islands_pin*/ 1, /*number of edges between two islands*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("One island, no edges:\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/1, /*size of islands*/ 4,
                  /*islands_pin*/ 0, /*number of edges between two islands*/ 2) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("One island, full graph:\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/1, /*size of islands*/ 4,
                  /*islands_pin*/ 1, /*number of edges between two islands*/ 2) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("Two islands, full graphs, no connections between islands:\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/2, /*size of islands*/ 4,
                  /*islands_pin*/ 1, /*number of edges between two islands*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("Three islands, full graphs, 20 connections between islands.\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/3, /*size of islands*/ 4,
                  /*islands_pin*/ 1, /*number of edges between two islands*/ 20) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 18 + 60);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Negative number of islands.\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/-2, /*size of islands*/ 4,
                  /*islands_pin*/ 1, /*number of edges between two islands*/ 0) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("Negative island size.\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/2, /*size of islands*/ -4,
                  /*islands_pin*/ 1, /*number of edges between two islands*/ 0) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("Probability out of range.\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/2, /*size of islands*/ 4,
                  /*islands_pin*/ 2, /*number of edges between two islands*/ 0) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("Negative number of edges between islands.\n");
    IGRAPH_ASSERT(igraph_simple_interconnected_islands_game(&g, /*number of islands*/2, /*size of islands*/ 4,
                  /*islands_pin*/ 1, /*number of edges between two islands*/ -3) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
