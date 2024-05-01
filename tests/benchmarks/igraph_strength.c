#include <igraph.h>

#include "bench.h"

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

void rand_weight_vec(igraph_vector_t *vec, const igraph_t *graph) {
    igraph_integer_t i, n = igraph_ecount(graph);
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

#define N 10000
#define M 10
#define REP 100

    igraph_barabasi_game(&g, N, 1, M, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    rand_weight_vec(&weights, &g);

    BENCH(" 1 igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );

    igraph_destroy(&g);

#undef N
#undef M
#undef REP

#define N 100000
#define M 10
#define REP 100

    igraph_barabasi_game(&g, N, 1, M, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    rand_weight_vec(&weights, &g);

    BENCH(" 2 igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );

    igraph_destroy(&g);

#undef N
#undef M
#undef REP

#define N 100000
#define M 100
#define REP 10

    igraph_barabasi_game(&g, N, 1, M, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    rand_weight_vec(&weights, &g);

    BENCH(" 3 igraph_strength(), preferential attachment n=" TOSTR(N) ", m=" TOSTR(M) ", " TOSTR(REP) "x",
          REPEAT(igraph_strength(&g, &strength, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), REP)
    );

    igraph_destroy(&g);

#undef N
#undef M
#undef REP

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&strength);

    return 0;
}
