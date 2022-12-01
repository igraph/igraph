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
    igraph_t graph1, graph2;
    igraph_bool_t same;

    igraph_rng_seed(igraph_rng_default(), 137); // Makes test deterministic

    // singleton with loop

    igraph_erdos_renyi_game_gnm_multi(&graph2, 1, 2, 0, 1, 1); // undirected with loops with mutiedges
    IGRAPH_ASSERT(igraph_vcount(&graph2) == 1);
    print_graph(&graph2);
    igraph_destroy(&graph2);

    // directed with loops

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 99, 1, 1, 0);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 99);
    IGRAPH_ASSERT(igraph_is_directed(&graph1));
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 560, 1, 1, 1);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 560);
    IGRAPH_ASSERT(igraph_is_directed(&graph1));
    igraph_simplify(&graph1, true, false, NULL);
    IGRAPH_ASSERT(igraph_ecount(&graph1) <= 100);
    igraph_destroy(&graph1);

    // directed without loops

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 89, 1, 0, 0);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 89);
    IGRAPH_ASSERT(igraph_is_directed(&graph1));
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 890, 1, 0, 1);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 890);
    IGRAPH_ASSERT(igraph_is_directed(&graph1));
    igraph_simplify(&graph1, true, false, NULL);
    IGRAPH_ASSERT(igraph_ecount(&graph1) <= 90);
    igraph_destroy(&graph1);

    // undirected with loops

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 54, 0, 1, 0);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 54);
    IGRAPH_ASSERT(!igraph_is_directed(&graph1));
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 540, 0, 1, 1);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 540);
    IGRAPH_ASSERT(!igraph_is_directed(&graph1));
    igraph_simplify(&graph1, true, false, NULL);
    IGRAPH_ASSERT(igraph_ecount(&graph1) <= 55);
    igraph_destroy(&graph1);

    // undirected without loops

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 44, 0, 0, 0);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 44);
    IGRAPH_ASSERT(!igraph_is_directed(&graph1));
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm_multi(&graph1, 10, 440, 0, 0, 1);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 440);
    IGRAPH_ASSERT(!igraph_is_directed(&graph1));
    igraph_copy(&graph2, &graph1);
    igraph_simplify(&graph1, true, true, NULL);
    igraph_simplify(&graph2, true, false, NULL);
    igraph_isomorphic(&graph1, &graph2, &same);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph1);
    igraph_destroy(&graph2);

    // large graphs

    igraph_erdos_renyi_game_gnm_multi(&graph1, 100000, 2.0 * 100000, 1, 0, 1);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 2.0 * 100000);
    IGRAPH_ASSERT(igraph_is_directed(&graph1));
    igraph_destroy(&graph1);


    igraph_erdos_renyi_game_gnm_multi(&graph1, 100000, 2.0 * 100000, 0, 0, 1);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 2.0 * 100000);
    IGRAPH_ASSERT(!igraph_is_directed(&graph1));
    igraph_destroy(&graph1);

    
    VERIFY_FINALLY_STACK();
    return 0;
}
