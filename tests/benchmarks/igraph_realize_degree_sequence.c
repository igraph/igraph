/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

void bench(const igraph_t *graph, const char *what, int rep) {
    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    igraph_vector_int_t degrees;
    igraph_t new_graph;

    igraph_vector_int_init(&degrees, 0);
    igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    char name[200];
    sprintf(name, "%s, n=%5" IGRAPH_PRId ", m=%5" IGRAPH_PRId ", LARGEST, %dx", what, vcount, ecount, rep);
    BENCH(name, REPEAT(igraph_realize_degree_sequence(&new_graph, &degrees, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST), rep));
    igraph_destroy(&new_graph);

    sprintf(name, "%s, n=%5" IGRAPH_PRId ", m=%5" IGRAPH_PRId ", SMALLEST, %dx", what, vcount, ecount, rep);
    BENCH(name, REPEAT(igraph_realize_degree_sequence(&new_graph, &degrees, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST), rep));
    igraph_destroy(&new_graph);

    sprintf(name, "%s, n=%5" IGRAPH_PRId ", m=%5" IGRAPH_PRId ", INDEX, %dx", what, vcount, ecount, rep);
    BENCH(name, REPEAT(igraph_realize_degree_sequence(&new_graph, &degrees, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_INDEX), rep));
    igraph_destroy(&new_graph);

    printf("\n");

    igraph_vector_int_destroy(&degrees);
}

void bench_gnm(igraph_int_t n, igraph_int_t m, int rep) {
    igraph_t graph;
    igraph_erdos_renyi_game_gnm(&graph, n, m, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    bench(&graph, "G(n,m)", rep);
    igraph_destroy(&graph);
}

void bench_ba(igraph_int_t n, igraph_int_t m, int rep) {
    igraph_t graph;
    igraph_barabasi_game(&graph, n, 1, m, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    bench(&graph, "BA", rep);
    igraph_destroy(&graph);
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 789);
    BENCH_INIT();

    bench_gnm(100, 200, 10000);
    bench_gnm(100, 1000, 10000);

    bench_gnm(10000, 20000, 1);
    bench_gnm(10000, 100000, 1);

    bench_ba(100, 1, 10000);
    bench_ba(100, 10, 10000);

    bench_ba(10000, 1, 1);
    bench_ba(10000, 10, 1);

    return 0;
}
