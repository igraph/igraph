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
    igraph_t graph;
    igraph_vector_t betweenness;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&betweenness, 0);

    /* Kautz and De Bruijn graphs are connected, therefore there should not be a dramatic difference
     * in the performance of the directed and undirected calculations. */

    igraph_kautz(&graph, 4, 5);
    BENCH(" 1 Betweenness, Kautz(4,5), directed",
          igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false));
    BENCH(" 2 Betweenness, Kautz(4,5), undirected",
          igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false));
    igraph_destroy(&graph);

    igraph_de_bruijn(&graph, 6, 5);
    BENCH(" 3 Betweenness, DeBruijn(6,5), directed",
          igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false));
    igraph_destroy(&graph);

    {
        igraph_vector_int_t dims;

        igraph_vector_int_init_int_end(&dims, -1, 8, 8, 8, 8, -1);
        igraph_square_lattice(&graph, &dims, 1, IGRAPH_UNDIRECTED, /* mutual */ 0, /* periodic */ 0);
        BENCH(" 4 Betweenness, Grid(8,8,8,8), undirected",
              igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false));
        igraph_destroy(&graph);
        igraph_vector_int_destroy(&dims);

        igraph_vector_int_init_int_end(&dims, -1, 10, 10, 10, 10, -1);
        igraph_square_lattice(&graph, &dims, 1, IGRAPH_UNDIRECTED, /* mutual */ 0, /* periodic */ 0);
        BENCH(" 5 Betweenness, Grid(10,10,10,10), cutoff 5",
              igraph_betweenness_cutoff(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false, 5));
        BENCH(" 6 Betweenness, Grid(10,10,10,10), cutoff 8",
              igraph_betweenness_cutoff(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false, 8));
        igraph_destroy(&graph);
        igraph_vector_int_destroy(&dims);
    }

    igraph_barabasi_game(&graph, 8000, 1, 1, NULL, 1, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 7 Betweenness, Barabasi n=8000 m=1, undirected",
          igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false));
    BENCH(" 8 Betweenness, Barabasi n=8000 m=1, undirected, cutoff 6",
          igraph_betweenness_cutoff(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false, 6));
    igraph_destroy(&graph);

    igraph_barabasi_game(&graph, 30000, 1, 5, NULL, 1, 0, IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 9 Betweenness, Barabasi n=30000 m=5, directed",
          igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false));
    BENCH("10 Betweenness, Barabasi n=30000 m=5, directed, cutoff 5",
          igraph_betweenness_cutoff(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false, 5));
    igraph_destroy(&graph);

    igraph_vector_destroy(&betweenness);

    return 0;
}
