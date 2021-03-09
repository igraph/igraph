/* -*- mode: C -*-  */
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
    igraph_t graph, comp_graph;
    igraph_bool_t is_perfect;

    igraph_rng_seed(igraph_rng_default(), 0);

    //bipartite 
    //==========================================================
    igraph_bipartite_game(&graph, NULL, IGRAPH_ERDOS_RENYI_GNM, 10, 10, 0, 20, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect);
    igraph_destroy(&graph);

    //complement to star graph size 10 - chordal
    //==========================================================
    igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_complementer(&comp_graph, &graph, 0);
    igraph_is_perfect(&comp_graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect);
    igraph_destroy(&graph);
    igraph_destroy(&comp_graph);

    
    //A cycle of size 5
    //==========================================================
    igraph_ring(&graph, 5, IGRAPH_UNDIRECTED, 0, 1);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(!is_perfect);
    igraph_destroy(&graph);

    //Paley graph of order 9
    //==========================================================
    igraph_small(&graph, 9, IGRAPH_UNDIRECTED,
                0, 1, 0, 3, 0, 6, 0, 2, 1, 2, 1, 4, 1, 7, 2, 5, 2, 8,
                3, 4, 3, 5, 3, 6, 4, 5, 4, 7, 5, 8, 6, 7, 7, 8, 6, 8, -1);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect);
    igraph_destroy(&graph);

    //Null graph
    //==========================================================
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect);
    igraph_destroy(&graph);

    //singleton graph
    //==========================================================
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect);
    igraph_destroy(&graph);

    //Empty graph
    //==========================================================
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect);
    igraph_destroy(&graph);

    //Test directed paths
    igraph_bipartite_game(&graph, NULL, IGRAPH_ERDOS_RENYI_GNM, 10, 10, 0, 20, IGRAPH_DIRECTED, IGRAPH_ALL);
    igraph_set_error_handler(igraph_error_handler_ignore);
    IGRAPH_ASSERT(igraph_is_perfect(&graph, &is_perfect) == IGRAPH_EINVAL);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
