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

void bench_gnp(igraph_int_t n, igraph_real_t p, igraph_bool_t directed, igraph_real_t negloop, int rep) {
    igraph_t g;
    igraph_matrix_t res;
    igraph_vector_t weights;
    char msg[200];

    igraph_matrix_init(&res, 0, 0);
    igraph_vector_init(&weights, 0);

    /* IGRAPH_NO_LOOPS, as loops are not allowed with negative weights. */
    igraph_erdos_renyi_game_gnp(&g, n, p, directed, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    igraph_matrix_resize(&res, igraph_vcount(&g), igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    for (igraph_int_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 1 n=%d, p=%g, %s, Unweighted, %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances(&g, NULL, &res, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT), rep);
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 2 n=%d, p=%g, %s, Dijkstra, %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_dijkstra(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), rep)
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 3 n=%d, p=%g, %s, Unweighted Floyd-Warshall, %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), rep)
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 4 n=%d, p=%g, %s, Unweighted Floyd-Warshall-tree-speedup, %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), rep)
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 5 n=%d, p=%g, %s, Floyd-Warshall, %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), rep)
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 6 n=%d, p=%g, %s, Floyd-Warshall-tree-speedup, %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), rep)
    );

    /* Negative weights, only for directed graphs. */

    if (! directed) {
        return;
    }

    /* Add a small number of negative weights, rely on luck, as well as tuning 'negloop',
     * to avoid negative loops. */
    for (igraph_int_t i=0; i < trunc(negloop * igraph_ecount(&g)); i++) {
        /* For reproducibility, do not write two RNG_...() calls within the same statement,
         * as the C language does not guarantee any evaluation order between them. */
        igraph_real_t w = RNG_UNIF(-negloop, 0);
        VECTOR(weights)[RNG_INTEGER(0, igraph_ecount(&g)-1)] = w;
    }

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 7 n=%d, p=%g, %s, Bellman-Ford (negative), %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_bellman_ford(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), rep)
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 8 n=%d, p=%g, %s, Johnson (negative), %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_johnson(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), rep)
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 9 n=%d, p=%g, %s, Floyd-Warshall (negative), %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), rep)
    );

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "10 n=%d, p=%g, %s, Floyd-Warshall-tree-speedup (negative), %dx",
             (int) n, p, directed ? "  directed" : "undirected", rep);
    BENCH(msg,
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), rep)
    );
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    /*
    bench_gnp(100, 0.5, IGRAPH_DIRECTED, 0.01, 200); printf("\n");
    bench_gnp(100, 0.5, IGRAPH_UNDIRECTED, 0.01, 200); printf("\n");

    bench_gnp(30, 0.5, IGRAPH_DIRECTED, 0.01, 2000); printf("\n");
    bench_gnp(30, 0.5, IGRAPH_UNDIRECTED, 0.01, 2000); printf("\n");

    bench_gnp(100, 0.1, IGRAPH_DIRECTED, 0.01, 1000); printf("\n");
    bench_gnp(100, 0.1, IGRAPH_UNDIRECTED, 0.01, 1000); printf("\n");

    bench_gnp(500, 0.1, IGRAPH_DIRECTED, 0.005, 10); printf("\n");
    bench_gnp(500, 0.1, IGRAPH_UNDIRECTED, 0.005, 10); printf("\n");

    bench_gnp(1500, 0.02, IGRAPH_DIRECTED, 0.01, 1); printf("\n");
    bench_gnp(1500, 0.02, IGRAPH_UNDIRECTED, 0.01, 1); printf("\n");
    */

    bench_gnp(1000, 0.1, IGRAPH_DIRECTED, 0.002, 1); printf("\n");

    bench_gnp(300, 0.1, IGRAPH_DIRECTED, 0.005, 50); printf("\n");

    bench_gnp(100, 0.1, IGRAPH_DIRECTED, 0.01, 500); printf("\n");

    bench_gnp(50, 0.1, IGRAPH_DIRECTED, 0.01, 10000); printf("\n");

    bench_gnp(30, 0.1, IGRAPH_DIRECTED, 0.01, 20000); printf("\n");

    bench_gnp(20, 0.2, IGRAPH_DIRECTED, 0.01, 40000); printf("\n");

    bench_gnp(10, 0.2, IGRAPH_DIRECTED, 0.01, 100000); printf("\n");

    bench_gnp(1500, 0.01, IGRAPH_DIRECTED, 0.01, 1); printf("\n");

    bench_gnp(100, 0.5, IGRAPH_DIRECTED, 0.01, 1000); printf("\n");

    return 0;
}
