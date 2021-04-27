#include <igraph.h>

#include "bench.h"

void rand_weight_vec(igraph_vector_t *vec, const igraph_t *graph) {
    long i, n = igraph_ecount(graph);
    igraph_vector_resize(vec, n);
    for (i=0; i < n; ++i) {
        VECTOR(*vec)[i] = RNG_UNIF(1, 10);
    }
}

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

int main() {
    igraph_t graph;
    igraph_vector_t closeness, weight;

    /* These closeness benchmarks are identical to the weighted betweenness ones. */

    /* This benchmark compares directed/undirected and weighted/unweighted calculations
     * on the same graphs. */

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&closeness, 0);
    igraph_vector_init(&weight, 0);

    igraph_kautz(&graph, 4, 3);

    /* Kautz and De Bruijn graphs are connected, therefore there should not be a dramatic difference
     * in the performance of the directed and undirected calculations. */

#define NAME "Kautz(4,3)"
#define REP 100

    BENCH(" 1 Closeness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, NULL, 1), REP)
    );
    BENCH(" 2 Closeness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, NULL, 1), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH(" 3 Closeness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, &weight, 1), REP)
    );
    BENCH(" 4 Closeness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, &weight, 1), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

#define NAME "DeBruijn(5,5)"
#define REP 1

    igraph_de_bruijn(&graph, 5, 5);

    BENCH(" 5 Closeness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, NULL, 1), REP)
    );
    BENCH(" 6 Closeness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, NULL, 1), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH(" 7 Closeness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, &weight, 1), REP)
    );
    BENCH(" 8 Closeness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, &weight, 1), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

    /* Choose the parameters of the Erdős-Rényi model so that the graph will have a large strongly connected
     * giant component. With 3000 vertices and 10000 edges, it is likely to contain over 90% of vertices.
     * If the graph does not have a giant component, then the directed betweenness calculation will be very
     * fast, and not directly comparable to the undirected one. */

#define NAME "GNM(3000,10000)"
#define REP 1

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 3000, 10000, IGRAPH_DIRECTED, IGRAPH_LOOPS);

    BENCH(" 9 Closeness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, NULL, 1), REP)
    );
    BENCH("10 Closeness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, NULL, 1), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH("11 Closeness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, &weight, 1), REP)
    );
    BENCH("12 Closeness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, &weight, 1), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

    /* Benchmark a much denser Erdős-Rényi graph as well. */

#define NAME "GNM(3000,30000)"
#define REP 1

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 3000, 30000, IGRAPH_DIRECTED, IGRAPH_LOOPS);

    BENCH("13 Closeness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, NULL, 1), REP)
    );
    BENCH("14 Closeness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, NULL, 1), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH("15 Closeness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_OUT, &weight, 1), REP)
    );
    BENCH("16 Closeness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_closeness(&graph, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, &weight, 1), REP)
    );

    igraph_destroy(&graph);

    igraph_vector_destroy(&weight);
    igraph_vector_destroy(&closeness);

    return 0;
}
