/*
   igraph library.
   Copyright (C) 2017-2025  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_bool_t valid_coloring;

    /* Initialize the library. */
    igraph_setup();

    /* Setting a seed makes the result of erdos_renyi_game_gnm deterministic. */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* IGRAPH_UNDIRECTED and IGRAPH_NO_LOOPS are both equivalent to 0/FALSE, but
       communicate intent better in this context. */
    igraph_erdos_renyi_game_gnm(&graph, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    /* As with all igraph functions, the vector in which the result is returned must
       be initialized in advance. */
    igraph_vector_int_init(&colors, 0);
    igraph_vertex_coloring_greedy(&graph, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);

    /* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
    igraph_is_vertex_coloring(&graph, &colors, &valid_coloring);
    IGRAPH_ASSERT(valid_coloring);

    /* Destroy data structure when we are done. */
    igraph_vector_int_destroy(&colors);
    igraph_destroy(&graph);

    return 0;
}
