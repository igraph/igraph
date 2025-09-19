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

void gnp(igraph_int_t n, igraph_real_t p, igraph_bool_t directed, igraph_edge_type_sw_t allowed_edge_types) {
    igraph_t g;
    igraph_erdos_renyi_game_gnp(&g, n, p, directed, allowed_edge_types, IGRAPH_EDGE_UNLABELED);
    igraph_destroy(&g);
}

void gnm(igraph_int_t n, igraph_real_t meandeg, igraph_bool_t directed, igraph_edge_type_sw_t allowed_edge_types) {
    igraph_t g;
    igraph_erdos_renyi_game_gnm(&g, n, round(directed ? n * meandeg : 0.5 * n * meandeg), directed, allowed_edge_types, IGRAPH_EDGE_UNLABELED);
    igraph_destroy(&g);
}

void chung_lu(const igraph_vector_t *outdeg, const igraph_vector_t *indeg, igraph_bool_t loops) {
    igraph_t g;
    igraph_chung_lu_game(&g, outdeg, indeg, loops, IGRAPH_CHUNG_LU_ORIGINAL);
    igraph_destroy(&g);
}

void run_bench(igraph_int_t vcount, igraph_real_t meandeg, igraph_int_t rep) {
    igraph_t g;
    igraph_real_t p = meandeg / vcount;
    igraph_vector_t outdeg, indeg;
    char msg[128], msg2[80];

    snprintf(msg2, sizeof(msg2) / sizeof(msg2[0]),
             "vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
             vcount, meandeg, rep);

    igraph_vector_init(&outdeg, 0);
    igraph_vector_init(&indeg, 0);

    igraph_erdos_renyi_game_gnp(&g, vcount, p, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_strength(&g, &outdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, NULL);
    igraph_destroy(&g);

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 1 G(n,p)   undirected, no loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnp(vcount, p, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 2 G(n,m)   undirected, no loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnm(vcount, meandeg, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 3 Chung-Lu undirected, no loops, %s", msg2);
    BENCH(msg, REPEAT(chung_lu(&outdeg, NULL, IGRAPH_NO_LOOPS), rep));

    printf("\n");

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 4 G(n,p)   undirected,    loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnp(vcount, p, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 5 G(n,m)   undirected,    loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnm(vcount, meandeg, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 6 Chung-Lu undirected,    loops, %s", msg2);
    BENCH(msg, REPEAT(chung_lu(&outdeg, NULL, IGRAPH_LOOPS), rep));

    printf("\n");

    igraph_erdos_renyi_game_gnp(&g, vcount, p, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_strength(&g, &outdeg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS, NULL);
    igraph_strength(&g, &indeg, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS, NULL);
    igraph_destroy(&g);

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 7 G(n,p)     directed, no loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnp(vcount, p, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 8 G(n,m)     directed, no loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnm(vcount, meandeg, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 9 Chung-Lu   directed, no loops, %s", msg2);
    BENCH(msg, REPEAT(chung_lu(&outdeg, &indeg, IGRAPH_NO_LOOPS), rep));

    printf("\n");

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "10 G(n,p)     directed,    loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnp(vcount, p, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "11 G(n,m)     directed,    loops / no multi, %s", msg2);
    BENCH(msg, REPEAT(gnm(vcount, meandeg, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "12 Chung-Lu   directed,    loops, %s", msg2);
    BENCH(msg, REPEAT(chung_lu(&outdeg, &indeg, IGRAPH_LOOPS), rep));

    printf("\n\n");

    igraph_vector_destroy(&indeg);
    igraph_vector_destroy(&outdeg);
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    run_bench(100, 3, 10000);
    run_bench(10000, 3, 100);
    run_bench(1000000, 3, 1);

    run_bench(10000, 1, 1000);
    run_bench(10000, 100, 10);

    return 0;
}
