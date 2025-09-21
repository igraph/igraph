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

void run_bench(const igraph_t *graph, const igraph_vector_t *weights,
               const char *name, int rep) {
    igraph_vector_int_t membership;
    igraph_vector_t vertex_weight;
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t ecount = igraph_ecount(graph);
    char msg[256], msg2[128];

    igraph_vector_int_init(&membership, vcount);
    igraph_vector_init(&vertex_weight, vcount);

    igraph_strength(graph, &vertex_weight, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, weights);

    snprintf(msg2, sizeof(msg2) / sizeof(msg2[0]),
             "%s, vcount=%" IGRAPH_PRId ", ecount=%" IGRAPH_PRId ", %s, %dx",
             name, vcount, ecount, weights == NULL ? "unweighted" : "weighted",
             rep);


    snprintf(msg, sizeof(msg) / sizeof(msg[0]), "1 Louvain, %s", msg2);
    BENCH(msg, REPEAT(igraph_community_multilevel(graph, weights, 1.0, &membership, NULL, NULL), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]), "2 Leiden , %s", msg2);
    BENCH(msg, REPEAT(igraph_community_leiden(graph, weights, &vertex_weight, NULL, 1.0 / igraph_vector_sum(weights), 0.01, false, 1, &membership, NULL, NULL), rep));

    printf("\n");

    igraph_vector_destroy(&vertex_weight);
    igraph_vector_int_destroy(&membership);
}

void rand_weights(const igraph_t *graph, igraph_vector_t *weights) {
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_vector_resize(weights, ecount);
    for (igraph_int_t i=0; i < ecount; i++) {
        VECTOR(*weights)[i] = RNG_UNIF01();
    }
}

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_init(&weights, 0);

    igraph_erdos_renyi_game_gnm(&graph, 100, 500, false, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 1000);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 1000, 5000, false, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 100);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 1000, 50000, false, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 10);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 10000, 50000, false, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 10);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 100000, 500000, false, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 1);
    igraph_destroy(&graph);

    igraph_forest_fire_game(&graph, 1000, 0.2, 1, 2, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "forest fire", 100);
    igraph_destroy(&graph);

    igraph_barabasi_game(&graph, 1000, 1, 5, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "PA", 100);
    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);

    return 0;
}
