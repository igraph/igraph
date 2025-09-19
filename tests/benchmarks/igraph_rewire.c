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

void bench(const char *name, igraph_t *graph, int rep) {
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_int_t trials = 10*ecount;
    char msg[200];

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%s, |V|=%6" IGRAPH_PRId ", |E|=%6" IGRAPH_PRId ", %6" IGRAPH_PRId " swaps, %dx",
             name, vcount, ecount, trials, rep);

    BENCH(msg, REPEAT(igraph_rewire(graph, trials, IGRAPH_SIMPLE_SW, NULL), rep));
}

void bench_er(igraph_int_t n, igraph_int_t m, int rep) {
    igraph_t graph;

    igraph_erdos_renyi_game_gnm(&graph, n, m, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    bench("G(n,m)", &graph, rep);

    igraph_destroy(&graph);
}

void bench_ba(igraph_int_t n, igraph_int_t m, int rep) {
    igraph_t graph;

    igraph_barabasi_game(&graph, n, 1, m ,NULL, true, 1, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_BAG, NULL);

    bench("BA", &graph, rep);

    igraph_destroy(&graph);
}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    bench_er(100, 3000, 100);
    bench_er(300, 3000, 100);
    bench_er(1000, 3000, 100);
    bench_er(3000, 3000, 100);
    bench_er(10000, 3000, 100);

    printf("\n");
    bench_er(300, 30000, 10);
    bench_er(1000, 30000, 10);
    bench_er(3000, 30000, 10);
    bench_er(10000, 30000, 10);

    printf("\n");
    bench_ba(1000, 3, 100);
    bench_ba(10000, 3, 10);
    bench_ba(100000, 3, 1);

    return 0;
}
