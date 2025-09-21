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

void rand_weight_vec(igraph_vector_t *vec, const igraph_t *graph) {
    igraph_int_t i, n = igraph_ecount(graph);
    igraph_vector_resize(vec, n);
    for (i=0; i < n; ++i) {
        VECTOR(*vec)[i] = RNG_UNIF(1, 10);
    }
}

int main(void) {

    igraph_t g;
    igraph_vector_t strength, weights;

    igraph_rng_seed(igraph_rng_default(), 54);
    BENCH_INIT();

    igraph_vector_init(&strength, 0);
    igraph_vector_init(&weights, 0);

    /* igraph_strength() uses an optimized, cache friendly code path when given
     * a vertex selector for which igraph_vs_is_all() returns true (this is
     * currently igraph_vs_all()). We pass both igraph_vss_all() and
     * igraph_vss_range(0, igraph_vcount(graph)) as vertex selectors
     * to compare the performance of the two code paths.
     *
     * NOTE: While currently igraph_vs_is_all() does not return true for
     * a range-type vertex selector, this may change in the future.
     * An altrenative is to to use igraph_vs_vector(), with the vector
     * initialized to a range.
     */

#define N 1000
#define M 10
#define REP 10000

    igraph_barabasi_game(&g, N, 1, M, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    rand_weight_vec(&weights, &g);

    BENCH(" 1a igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
          );
    BENCH(" 1b igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_range(0, igraph_vcount(&g)), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
          );
    printf("\n");

    igraph_destroy(&g);

#undef N
#undef M
#undef REP

#define N 10000
#define M 10
#define REP 1000

    igraph_barabasi_game(&g, N, 1, M, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    rand_weight_vec(&weights, &g);

    BENCH(" 2a igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );
    BENCH(" 2b igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_range(0, igraph_vcount(&g)), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );
    printf("\n");

    igraph_destroy(&g);

#undef N
#undef M
#undef REP

#define N 100000
#define M 10
#define REP 100

    igraph_barabasi_game(&g, N, 1, M, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    rand_weight_vec(&weights, &g);

    BENCH(" 3a igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );
    BENCH(" 3b igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_range(0, igraph_vcount(&g)), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );
    printf("\n");

    igraph_destroy(&g);

#undef N
#undef M
#undef REP

#define N 100000
#define M 100
#define REP 10

    igraph_barabasi_game(&g, N, 1, M, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    rand_weight_vec(&weights, &g);

    BENCH(" 4a igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );
    BENCH(" 4b igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_range(0, igraph_vcount(&g)), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );
    printf("\n");

    igraph_destroy(&g);

#undef N
#undef M
#undef REP

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&strength);

    return 0;
}
