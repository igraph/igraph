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

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

#define BENCH_UNDIR(N, REP) \
    do { \
        igraph_barabasi_game(&graph, N, 1, 3, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL); \
        igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); \
        igraph_destroy(&graph); \
        \
        BENCH("Undirected simple,       PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_SIMPLE_SW, &graphical), REP); \
              ); \
        \
        BENCH("Undirected loops,        PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_LOOPS_SW, &graphical), REP); \
              ); \
        \
        BENCH("Undirected multi,        PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_MULTI_SW, &graphical), REP); \
              ); \
        \
        BENCH("Undirected multi, loops, PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW, &graphical), 100); \
              ); \
        \
        igraph_erdos_renyi_game_gnp(&graph, N, 12.0/N, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED); \
        igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); \
        igraph_destroy(&graph); \
        \
        BENCH("Undirected simple,       G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_SIMPLE_SW, &graphical), REP); \
              ); \
            \
        BENCH("Undirected loops,        G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_LOOPS_SW, &graphical), REP); \
              ); \
            \
        BENCH("Undirected multi,        G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_MULTI_SW, &graphical), REP); \
              ); \
            \
        BENCH("Undirected multi, loops, G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&deg, NULL, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW, &graphical), 100); \
              ); \
            \
    } while (0)


#define BENCH_DIR(N, REP) \
do { \
        igraph_barabasi_game(&graph, N, 1, 3, NULL, true, 0, IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL); \
        igraph_degree(&graph, &outdeg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS); \
        igraph_degree(&graph, &indeg, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS); \
        igraph_destroy(&graph); \
        \
        BENCH("Directed simple,       PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_SIMPLE_SW, &graphical), REP); \
              ); \
        \
        BENCH("Directed loops,        PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_LOOPS_SW, &graphical), REP); \
              ); \
        \
        BENCH("Directed multi,        PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_MULTI_SW, &graphical), REP); \
              ); \
        \
        BENCH("Directed multi, loops, PA,     n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW, &graphical), REP); \
              ); \
        \
        igraph_erdos_renyi_game_gnp(&graph, N, 12.0/N, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED); \
        igraph_degree(&graph, &outdeg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS); \
        igraph_degree(&graph, &indeg, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS); \
        igraph_destroy(&graph); \
        \
        BENCH("Directed simple,       G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_SIMPLE_SW, &graphical), REP); \
              ); \
            \
        BENCH("Directed loops,        G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_LOOPS_SW, &graphical), REP); \
              ); \
            \
        BENCH("Directed multi,        G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_MULTI_SW, &graphical), REP); \
              ); \
            \
        BENCH("Directed multi, loops, G(n,p), n = " TOSTR(N) ", " TOSTR(REP) "x", \
              REPEAT(igraph_is_graphical(&outdeg, &indeg, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW, &graphical), REP); \
              ); \
        \
} while (0)


int main(void) {
    igraph_t graph;
    igraph_vector_int_t deg, outdeg, indeg;
    igraph_bool_t graphical;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_int_init(&deg, 0);
    igraph_vector_int_init(&outdeg, 0);
    igraph_vector_int_init(&indeg, 0);

    BENCH_UNDIR(50, 2000000);
    BENCH_DIR(50, 2000000);
    printf("\n");
    BENCH_UNDIR(1000, 100000);
    BENCH_DIR(1000, 100000);
    printf("\n");
    BENCH_UNDIR(1000000, 100);
    BENCH_DIR(1000000, 100);

    igraph_vector_int_destroy(&deg);
    igraph_vector_int_destroy(&outdeg);
    igraph_vector_int_destroy(&indeg);

    return 0;
}
