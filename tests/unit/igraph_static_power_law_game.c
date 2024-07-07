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

int main(void) {
    igraph_t g;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No vertices:\n");
    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 0, /*number of edges*/ 0,
                  /*exponent_out*/ 2.0, /*exponent in*/ 2.0, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE,
                  /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("No edges, undirected:\n");
    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 10, /*number of edges*/ 0,
                  /*exponent_out*/ 2.0, /*exponent in*/ -2.0, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE,
                  /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("Checking some basic outputs.\n");
    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 100, /*number of edges*/ 30,
                  /*exponent_out*/ 2.0, /*exponent in*/ -2.0, IGRAPH_LOOPS, IGRAPH_MULTIPLE,
                  /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) == 30);
    igraph_destroy(&g);

    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 90, /*number of edges*/ 40,
                                               /*exponent_out*/ 2.0, /*exponent in*/ -2.0, IGRAPH_LOOPS, IGRAPH_NO_MULTIPLE,
                                               /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 90);
    IGRAPH_ASSERT(igraph_ecount(&g) == 40);
    igraph_destroy(&g);

    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 110, /*number of edges*/ 50,
                                               /*exponent_out*/ 2.0, /*exponent in*/ -2.0, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE,
                                               /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 110);
    IGRAPH_ASSERT(igraph_ecount(&g) == 50);
    igraph_destroy(&g);

    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 100, /*number of edges*/ 30,
                                               /*exponent_out*/ 2.0, /*exponent in*/ 2.0, IGRAPH_LOOPS, IGRAPH_MULTIPLE,
                                               /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) == 30);
    igraph_destroy(&g);

    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 90, /*number of edges*/ 40,
                                               /*exponent_out*/ 2.0, /*exponent in*/ 2.2, IGRAPH_LOOPS, IGRAPH_NO_MULTIPLE,
                                               /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 90);
    IGRAPH_ASSERT(igraph_ecount(&g) == 40);
    igraph_destroy(&g);

    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 110, /*number of edges*/ 50,
                                               /*exponent_out*/ 2.0, /*exponent in*/ 2.5, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE,
                                               /*finite_size_correction*/ true) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 110);
    IGRAPH_ASSERT(igraph_ecount(&g) == 50);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Negative number of vertices.\n");
    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ -100, /*number of edges*/ 30,
                  /*exponent_out*/ 2.0, /*exponent in*/ -2.0, IGRAPH_LOOPS, IGRAPH_MULTIPLE,
                  /*finite_size_correction*/ true) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("Negative number of edges.\n");
    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 100, /*number of edges*/ -30,
                  /*exponent_out*/ 2.0, /*exponent in*/ -2.0, IGRAPH_LOOPS, IGRAPH_MULTIPLE,
                  /*finite_size_correction*/ true) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("Exponent out too low.\n");
    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 100, /*number of edges*/ 30,
                  /*exponent_out*/ 1.0, /*exponent in*/ -2.0, IGRAPH_LOOPS, IGRAPH_MULTIPLE,
                  /*finite_size_correction*/ true) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("Exponent in too low but not negative.\n");
    IGRAPH_ASSERT(igraph_static_power_law_game(&g, /*number of vertices*/ 100, /*number of edges*/ 30,
                  /*exponent_out*/ 2.0, /*exponent in*/ 0.5, IGRAPH_LOOPS, IGRAPH_MULTIPLE,
                  /*finite_size_correction*/ true) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
