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

void run_bench(int n) {
    igraph_t graph;
    igraph_matrix_t modmat;
    int rep;
    char msg[255];

    rep = 20000 / n;
    rep = rep*rep;

    igraph_matrix_init(&modmat, n, n);

    igraph_erdos_renyi_game_gnm(&graph, n, 10 * n, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "G(n,m),   directed, n = %5d, m = %7" IGRAPH_PRId ", %dx", n, igraph_ecount(&graph), rep);
    BENCH(msg, REPEAT(igraph_modularity_matrix(&graph, NULL, 1.0, &modmat, IGRAPH_DIRECTED), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "G(n,m), undirected, n = %5d, m = %7" IGRAPH_PRId ", %dx", n, igraph_ecount(&graph), rep);
    BENCH(msg, REPEAT(igraph_modularity_matrix(&graph, NULL, 1.0, &modmat, IGRAPH_UNDIRECTED), rep));

    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnp(&graph, n, 0.1, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "G(n,p),   directed, n = %5d, m = %7" IGRAPH_PRId ", %dx", n, igraph_ecount(&graph), rep);
    BENCH(msg, REPEAT(igraph_modularity_matrix(&graph, NULL, 1.0, &modmat, IGRAPH_DIRECTED), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "G(n,p), undirected, n = %5d, m = %7" IGRAPH_PRId ", %dx", n, igraph_ecount(&graph), rep);
    BENCH(msg, REPEAT(igraph_modularity_matrix(&graph, NULL, 1.0, &modmat, IGRAPH_UNDIRECTED), rep));
    igraph_destroy(&graph);

    igraph_matrix_destroy(&modmat);

}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    run_bench(100);
    run_bench(1000);
    run_bench(10000);

    return 0;
}
