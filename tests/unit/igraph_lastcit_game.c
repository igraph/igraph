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
    igraph_vector_t preference;

    igraph_rng_seed(igraph_rng_default(), 42);

    /*No nodes*/
    igraph_vector_init_int_end(&preference, -1, 1, 1, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 0, /*edges_per_node*/ 5, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);
    igraph_vector_destroy(&preference);

    /*No edges*/
    igraph_vector_init_int_end(&preference, -1, 1, 1, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 9, /*edges_per_node*/ 0, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 9);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);
    igraph_vector_destroy(&preference);

    /*Only cite un-cited to make a line*/
    igraph_vector_init_int_end(&preference, -1, 0, 1, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 9, /*edges_per_node*/ 1, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);
    igraph_vector_destroy(&preference);

    /*Hugely prefer cited to make a star*/
    igraph_vector_init_real(&preference, 2, 1e30, 1e-30);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 9, /*edges_per_node*/ 1, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);
    igraph_vector_destroy(&preference);

    igraph_set_error_handler(igraph_error_handler_ignore);

    /*Negative number of nodes*/
    igraph_vector_init_int_end(&preference, -1, 1, 1, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ -9, /*edges_per_node*/ 1, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&preference);

    /*Too few agebins*/
    igraph_vector_init_int_end(&preference, -1, 1, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 9, /*edges_per_node*/ 1, /*agebins*/ 0, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&preference);

    /*Wrong vector size*/
    igraph_vector_init_int_end(&preference, -1, 1, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 9, /*edges_per_node*/ 1, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&preference);

    /*No uncited preference*/
    igraph_vector_init_int_end(&preference, -1, 1, 0, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 9, /*edges_per_node*/ 1, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&preference);

    /*Negative preference*/
    igraph_vector_init_int_end(&preference, -1, -1, 1, -1);
    IGRAPH_ASSERT(igraph_lastcit_game(&g, /*nodes*/ 9, /*edges_per_node*/ 1, /*agebins*/ 1, /*preference*/ &preference, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&preference);

    VERIFY_FINALLY_STACK();
    return 0;
}
