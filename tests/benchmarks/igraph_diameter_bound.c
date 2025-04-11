/*
   IGraph library.
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

int main(void) {
    igraph_t g;
    igraph_real_t res, ref;
    igraph_vector_t weights;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_init(&weights, 0);

    printf("Sparse G(n, p)\n");

#define VCOUNT 250
#define DENS 3.0/VCOUNT
#define REP 1000

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Bound, " TOSTR(REP) "x",
          REPEAT(igraph_diameter_bound(&g, NULL, &res, IGRAPH_UNDIRECTED, true), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Original, " TOSTR(REP) "x",
          REPEAT(igraph_diameter(&g, &ref, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, true), REP);
    );
    IGRAPH_ASSERT(res == ref);
    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Weighted Bound, " TOSTR(REP) "x",
        REPEAT(igraph_diameter_bound(&g, &weights, &res, IGRAPH_UNDIRECTED, true), REP);
    );
    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Weighted Original, " TOSTR(REP) "x",
        REPEAT(igraph_diameter_dijkstra(&g, &weights, &ref, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, true), REP);
    );
    IGRAPH_ASSERT(igraph_almost_equals(res, ref, 1e-10));

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("Dense G(n, p)\n");

#define VCOUNT 250
#define DENS 0.3
#define REP 100

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

    BENCH(" 5 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Bound, " TOSTR(REP) "x",
          REPEAT(igraph_diameter_bound(&g, NULL, &res, IGRAPH_UNDIRECTED, true), REP);
    );
    BENCH(" 6 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Original, " TOSTR(REP) "x",
          REPEAT(igraph_diameter(&g, &ref, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, true), REP);
    );
    IGRAPH_ASSERT(res == ref);
    BENCH(" 7 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Weighted Bound, " TOSTR(REP) "x",
        REPEAT(igraph_diameter_bound(&g, &weights, &res, IGRAPH_UNDIRECTED, true), REP);
    );
    BENCH(" 8 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Weighted Original, " TOSTR(REP) "x",
        REPEAT(igraph_diameter_dijkstra(&g, &weights, &ref, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, true), REP);
    );
    IGRAPH_ASSERT(igraph_almost_equals(res, ref, 1e-10));

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("Barabasi-Albert G(n, p)\n");

#define VCOUNT 250
#define ECOUNT 15
#define REP 1000

    igraph_barabasi_game(&g, VCOUNT, 1, ECOUNT, 0, 0, 1, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_BAG, NULL);
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

    BENCH(" 9 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(ECOUNT) ", Unweighted Bound, " TOSTR(REP) "x",
          REPEAT(igraph_diameter_bound(&g, NULL, &res, IGRAPH_UNDIRECTED, true), REP);
    );
    BENCH("10 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(ECOUNT) ", Unweighted Original, " TOSTR(REP) "x",
          REPEAT(igraph_diameter(&g, &ref, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, true), REP);
    );
    IGRAPH_ASSERT(res == ref);
    BENCH("11 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(ECOUNT) ", Weighted Bound, " TOSTR(REP) "x",
        REPEAT(igraph_diameter_bound(&g, &weights, &res, IGRAPH_UNDIRECTED, true), REP);
    );
    BENCH("12 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(ECOUNT) ", Weighted Original, " TOSTR(REP) "x",
        REPEAT(igraph_diameter_dijkstra(&g, &weights, &ref, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, true), REP);
    );
    IGRAPH_ASSERT(igraph_almost_equals(res, ref, 1e-10));

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

    igraph_vector_destroy(&weights);

    return 0;
}
