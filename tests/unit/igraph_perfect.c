/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2021  Gabor Csardi <csardi.gabor@gmail.com>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include "igraph_structural.h"

#include "test_utilities.inc"

int main() {
    igraph_t graph, comp_graph;
    igraph_bool_t is_perfect;
    int x;

    printf("bipartite \n");
    printf("==========================================================\n");
    igraph_bipartite_game(&graph, NULL, IGRAPH_ERDOS_RENYI_GNM, 10, 10, 0, 20, 0, IGRAPH_ALL);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect == 1);
    igraph_destroy(&graph);
    x = IGRAPH_FINALLY_STACK_SIZE();

    printf("complanetry to star graph size 10 - chordal\n");
    printf("==========================================================\n");
    igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_complementer(&comp_graph, &graph, 0);
    igraph_is_perfect(&comp_graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect == 1);
    igraph_destroy(&graph);
    igraph_destroy(&comp_graph);
    x = IGRAPH_FINALLY_STACK_SIZE();

    
    printf("A cycle of size 5\n");
    printf("==========================================================\n");
    igraph_ring(&graph, 5, 0, 0, 1);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect == 0);
    igraph_destroy(&graph);
    x = IGRAPH_FINALLY_STACK_SIZE();

    printf("Paley graph of order 9\n");
    printf("==========================================================\n");
    igraph_small(&graph, 9, 0,
                0, 1, 0, 3, 0, 6, 0, 2, 1, 2, 1, 4, 1, 7, 2, 5, 2, 8,
                3, 4, 3, 5, 3, 6, 4, 5, 4, 7, 5, 8, 6, 7, 7, 8, 6, 8, -1);
    igraph_is_perfect(&graph, &is_perfect);
    IGRAPH_ASSERT(is_perfect == 1);
    igraph_destroy(&graph);
    x = IGRAPH_FINALLY_STACK_SIZE();

    VERIFY_FINALLY_STACK();

    return 0;
}
