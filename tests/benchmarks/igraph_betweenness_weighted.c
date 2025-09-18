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

void rand_weight_vec(igraph_vector_t *vec, const igraph_t *graph) {
    const igraph_int_t n = igraph_ecount(graph);
    igraph_vector_resize(vec, n);
    for (igraph_int_t i=0; i < n; ++i) {
        VECTOR(*vec)[i] = RNG_UNIF(1, 10);
    }
}

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

int main(void) {
    igraph_t graph;
    igraph_vector_t betweenness, weight;

    /* These betweenness benchmarks are identical to the weighted closeness ones. */

    /* This benchmark compares directed/undirected and weighted/unweighted calculations
     * on the same graphs. */

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&betweenness, 0);
    igraph_vector_init(&weight, 0);

    igraph_kautz(&graph, 4, 3);

    /* Kautz and De Bruijn graphs are connected, therefore there should not be a dramatic difference
     * in the performance of the directed and undirected calculations. */

#define NAME "Kautz(4,3)"
#define REP 100

    BENCH(" 1 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH(" 2 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH(" 3 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH(" 4 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

#define NAME "DeBruijn(5,5)"
#define REP 1

    igraph_de_bruijn(&graph, 5, 5);

    BENCH(" 5 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH(" 6 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH(" 7 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH(" 8 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

    /* Choose the parameters of the Erdős-Rényi model so that the graph will have a large strongly connected
     * giant component. With 3000 vertices and 10000 edges, it is likely to contain over 90% of vertices.
     * If the graph does not have a giant component, then the directed betweenness calculation will be very
     * fast, and not directly comparable to the undirected one. */

#define NAME "GNM(3000,10000)"
#define REP 1

    igraph_erdos_renyi_game_gnm(&graph, 3000, 10000, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);

    BENCH(" 9 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH("10 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH("11 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH("12 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

    /* Benchmark a much denser Erdős-Rényi graph as well. */

#define NAME "GNM(3000,30000)"
#define REP 1

    igraph_erdos_renyi_game_gnm(&graph, 3000, 30000, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);

    BENCH("13 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH("14 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, NULL, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH("15 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, false), REP)
    );
    BENCH("16 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &weight, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, false), REP)
    );

    igraph_destroy(&graph);

    igraph_vector_destroy(&weight);
    igraph_vector_destroy(&betweenness);

    return 0;
}
