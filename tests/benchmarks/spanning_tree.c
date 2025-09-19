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

void run_bench(const igraph_t *graph, const igraph_vector_t *weights, int rep, const char *name) {
    igraph_vector_int_t edges;
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t ecount = igraph_ecount(graph);
    char msg[128];

    igraph_vector_int_init(&edges, vcount - 1);

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%s, vcount=%" IGRAPH_PRId ", ecount=%" IGRAPH_PRId ", unweigthed, %dx",
             name, vcount, ecount, rep);
    BENCH(msg, REPEAT(igraph_minimum_spanning_tree(graph, &edges, NULL, IGRAPH_MST_UNWEIGHTED), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%s, vcount=%" IGRAPH_PRId ", ecount=%" IGRAPH_PRId ", Prim, %dx",
             name, vcount, ecount, rep);
    BENCH(msg, REPEAT(igraph_minimum_spanning_tree(graph, &edges, weights, IGRAPH_MST_UNWEIGHTED), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%s, vcount=%" IGRAPH_PRId ", ecount=%" IGRAPH_PRId ", Kruskal, %dx",
             name, vcount, ecount, rep);
    BENCH(msg, REPEAT(igraph_minimum_spanning_tree(graph, &edges, weights, IGRAPH_MST_UNWEIGHTED), rep));

    printf("\n");

    igraph_vector_int_destroy(&edges);
}

void rand_weights(const igraph_t *graph, igraph_vector_t *weights) {
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_vector_resize(weights, ecount);
    for (igraph_int_t i=0; i < ecount; i++) {
        VECTOR(*weights)[i] = RNG_UNIF01();
    }
}

int main(void) {
    igraph_t g;
    igraph_vector_t weights;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_init(&weights, 0);

    igraph_barabasi_game(&g, 100, 1, 1, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 50000, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 100, 1, 5, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 10000, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 1000, 1, 5, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 1000, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 10000, 1, 5, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 100, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 100, 1, 50, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 1000, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 1000, 1, 50, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 100, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 10000, 1, 50, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 10, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 1000, 1, 500, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 10, "PA");
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 10000, 1, 500, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 1, "PA");
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 1000, 500, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 60000, "G(n,m)");
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 1000, 1000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 30000, "G(n,m)");
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 1000, 3000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 10000, "G(n,m)");
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 3000, "G(n,m)");
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 1000, 30000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 1000, "G(n,m)");
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 1000, 100000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 300, "G(n,m)");
    igraph_destroy(&g);

    igraph_full(&g, 100, false, false);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 10000, "K");
    igraph_destroy(&g);

    igraph_full(&g, 1000, false, false);
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 100, "K");
    igraph_destroy(&g);

    {
        igraph_t K;
        igraph_full(&K, 100, false, false);
        igraph_disjoint_union(&g, &K, &K);
        igraph_destroy(&K);
    }
    rand_weights(&g, &weights);
    run_bench(&g, &weights, 10000, "K + K");
    igraph_destroy(&g);

    return 0;
}
