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

int main(void) {
    igraph_t g;
    igraph_vector_int_t membership;
    igraph_vector_t distances;
    igraph_vector_t weights;
    igraph_vector_int_t generators;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_int_init(&membership, 0);
    igraph_vector_init(&distances, 0);
    igraph_vector_init(&weights, 0);
    igraph_vector_int_init(&generators, 0);

#define VCOUNT 100
#define DENS 0.5
#define REP 1000

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    for (igraph_int_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }

    igraph_random_sample(&generators, 0, igraph_vcount(&g)-1, 20);

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Dijkstra, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, &weights, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 500
#define DENS 0.1
#define REP 1000

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    for (igraph_int_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }

    igraph_random_sample(&generators, 0, igraph_vcount(&g)-1, 20);

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Dijkstra, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, &weights, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 5000
#define DENS 0.01
#define REP 100

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    for (igraph_int_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }

    igraph_random_sample(&generators, 0, igraph_vcount(&g)-1, 20);

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Dijkstra, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, &weights, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 100000
#define DENS 0.0001
#define REP 1

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    for (igraph_int_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }

    igraph_random_sample(&generators, 0, igraph_vcount(&g)-1, 100);

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Dijkstra, " TOSTR(REP) "x",
          REPEAT(igraph_voronoi(&g, &membership, &distances, &generators, &weights, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM), REP);
    );

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    igraph_vector_int_destroy(&generators);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&distances);
    igraph_vector_int_destroy(&membership);

    return 0;
}
