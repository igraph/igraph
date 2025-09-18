/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include "bench.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_t colors;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_vector_int_init(&colors, 0);

    igraph_erdos_renyi_game_gnm(&g, 50000, 1000000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    BENCH(" 1 COLORED_NEIGHBORS, random graph with 50,000 vertices and 2,000,000 edges",
          igraph_vertex_coloring_greedy(&g, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS)
         );
    BENCH(" 2 DSATUR,            random graph with 50,000 vertices and 2,000,000 edges",
          igraph_vertex_coloring_greedy(&g, &colors, IGRAPH_COLORING_GREEDY_DSATUR)
         );
    igraph_destroy(&g);

    printf("\n");

    igraph_barabasi_game(&g, 100000, 1, 15, NULL, 0, 0, 0, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 1 COLORED_NEIGHBORS, pref. attach. graph n=100,000 m=15",
          igraph_vertex_coloring_greedy(&g, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS)
         );
    BENCH(" 2 DSATUR,            pref. attach. graph n=100,000 m=15",
          igraph_vertex_coloring_greedy(&g, &colors, IGRAPH_COLORING_GREEDY_DSATUR)
         );
    igraph_destroy(&g);

    igraph_vector_int_destroy(&colors);


    return 0;
}
