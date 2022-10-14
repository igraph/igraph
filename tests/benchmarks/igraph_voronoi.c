
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

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

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

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

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

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

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

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_int_resize(&membership, igraph_vcount(&g));
    igraph_vector_resize(&distances, igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

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
