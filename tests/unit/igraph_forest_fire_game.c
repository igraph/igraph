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

    printf("No vertices:\n");
    IGRAPH_ASSERT(igraph_forest_fire_game(&g, /*number of vertices*/0, /*fw_prob*/ 0.0,
                  /*bw_factor*/ 0.0, /*pambs*/1, /*directed*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("No ambassadors:\n");
    IGRAPH_ASSERT(igraph_forest_fire_game(&g, /*number of vertices*/10, /*fw_prob*/ 0.0,
                  /*bw_factor*/ 0.0, /*pambs*/0, /*directed*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("More ambassadors than nodes:\n");
    IGRAPH_ASSERT(igraph_forest_fire_game(&g, /*number of vertices*/5, /*fw_prob*/ 0.0,
                  /*bw_factor*/ 0.0, /*pambs*/100, /*directed*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("Some normal inputs, just to check for memory problems, no output checking.\n");
    IGRAPH_ASSERT(igraph_forest_fire_game(&g, /*number of vertices*/50, /*fw_prob*/ 0.5,
                  /*bw_factor*/ 0.5, /*pambs*/3, /*directed*/ 1) == IGRAPH_SUCCESS);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Negative fw_prob.\n");
    IGRAPH_ASSERT(igraph_forest_fire_game(&g, /*number of vertices*/5, /*fw_prob*/ -0.5,
                  /*bw_factor*/ 0.0, /*pambs*/100, /*directed*/ 1) == IGRAPH_EINVAL);

    printf("Negative bw_factor.\n");
    IGRAPH_ASSERT(igraph_forest_fire_game(&g, /*number of vertices*/5, /*fw_prob*/ 0.5,
                  /*bw_factor*/ -0.5, /*pambs*/100, /*directed*/ 0) == IGRAPH_EINVAL);

    printf("Negative number of ambassadors.\n");
    IGRAPH_ASSERT(igraph_forest_fire_game(&g, /*number of vertices*/5, /*fw_prob*/ 0.5,
                  /*bw_factor*/ 0.5, /*pambs*/-100, /*directed*/ 1) == IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();
    return 0;
}
