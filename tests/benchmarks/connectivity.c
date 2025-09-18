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

/* Benchmarks related to biconnected components. */

void run_bench(igraph_int_t vcount, igraph_real_t meandeg, igraph_int_t rep) {
    igraph_t g;
    igraph_bool_t b;
    igraph_vector_int_t ivec;
    char msg[128];

    igraph_erdos_renyi_game_gnm(&g, vcount, round(meandeg * vcount / 2), IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_int_init(&ivec, igraph_vcount(&g));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 1 Weakly connected?   vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
              vcount, meandeg, 10*rep);

    igraph_invalidate_cache(&g);
    BENCH(msg, REPEAT(igraph_is_connected(&g, &b, IGRAPH_WEAK), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 2 Weak components.    vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
             vcount, meandeg, rep);

    igraph_invalidate_cache(&g);
    BENCH(msg, REPEAT(igraph_connected_components(&g, &ivec, NULL, NULL, IGRAPH_WEAK), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 3 Strongly connected? vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
             vcount, meandeg, 10*rep);

    igraph_invalidate_cache(&g);
    BENCH(msg, REPEAT(igraph_is_connected(&g, &b, IGRAPH_STRONG), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 4 Strong components.  vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
             vcount, meandeg, rep);

    igraph_invalidate_cache(&g);
    BENCH(msg, REPEAT(igraph_connected_components(&g, &ivec, NULL, NULL, IGRAPH_STRONG), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 5 Biconnected?        vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
             vcount, meandeg, rep);

    igraph_invalidate_cache(&g);
    BENCH(msg, REPEAT(igraph_is_biconnected(&g, &b), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 6 Bridges.            vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
             vcount, meandeg, rep);

    igraph_invalidate_cache(&g);
    BENCH(msg, REPEAT(igraph_bridges(&g, &ivec), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             " 7 Articulation pts.   vcount=%" IGRAPH_PRId ", meandeg=%g, %" IGRAPH_PRId "x",
             vcount, meandeg, rep);

    igraph_invalidate_cache(&g);
    BENCH(msg, REPEAT(igraph_articulation_points(&g, &ivec), rep));

    igraph_vector_int_destroy(&ivec);
    igraph_destroy(&g);

    printf("\n");

}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    // Note that whether these random graphs end up being connected
    // with high probability depends not only on their mean degree,
    // but also their vertex count.

    run_bench(10000, 0.5, 100); // no giant component
    run_bench(10000, 3, 100);   // not weakly connected
    run_bench(10000, 10, 100);  // not strongly connected
    run_bench(10000, 20, 100);  // strongly connnected

    run_bench(100000, 0.5, 100); // no giant component
    run_bench(100000, 3, 100);   // not weakly connected
    run_bench(100000, 10, 100);  // not strongly connected
    run_bench(100000, 20, 100);  // strongly connnected

    return 0;
}
