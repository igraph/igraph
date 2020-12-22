/*
   IGraph library.
   Copyright (C) 2017-2020  The igraph development team <igraph@igraph.org>

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
#include <assert.h>

int main() {
    igraph_t graph;
    igraph_vector_int_t colors;

    /* Setting a seed makes the result of erdos_renyi_game deterministic. */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* IGRAPH_UNDIRECTED and IGRAPH_NO_LOOPS are both equivalent to 0/FALSE, but
       communicate intent better in this context. */
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    /* As with all igraph functions, the vector in which the result is returned must
       be initialized in advance. */
    igraph_vector_int_init(&colors, 0);
    igraph_vertex_coloring_greedy(&graph, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);

    /* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
    {
        long int i;
        /* Store the edge count to avoid the overhead from igraph_ecount in the for loop. */
        long int no_of_edges = igraph_ecount(&graph);
        for (i = 0; i < no_of_edges; ++i) {
            assert( VECTOR(colors)[ IGRAPH_FROM(&graph, i) ] != VECTOR(colors)[ IGRAPH_TO(&graph, i) ]  );
        }
    }

    /* Destroy data structure when we are done. */
    igraph_vector_int_destroy(&colors);
    igraph_destroy(&graph);

    return 0;
}
