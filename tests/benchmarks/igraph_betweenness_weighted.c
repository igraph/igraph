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
    igraph_vector_t betweenness, weight;

    /* These betweenness benchmarks are identical to the weighted closeness ones. */

    /* This benchmark compares directed/undirected and weighted/unweighted calculations
     * on the same graphs. */

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&betweenness, 0);
    igraph_vector_init(&weight, 0);

    igraph_kautz(&graph, 4, 3);

    /* Kautz and De Bruijn graphs are connected, therefore there should not be a dramatic difference
     * in the performance of the directed and undirected calculations. */

#define NAME "Kautz(4,3)"
#define REP 100

    BENCH(" 1 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL), REP)
    );
    BENCH(" 2 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH(" 3 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, &weight), REP)
    );
    BENCH(" 4 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, &weight), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

#define NAME "DeBruijn(5,5)"
#define REP 1

    igraph_de_bruijn(&graph, 5, 5);

    BENCH(" 5 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL), REP)
    );
    BENCH(" 6 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH(" 7 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, &weight), REP)
    );
    BENCH(" 8 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, &weight), REP)
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

    BENCH(" 9 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL), REP)
    );
    BENCH("10 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH("11 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, &weight), REP)
    );
    BENCH("12 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, &weight), REP)
    );

    igraph_destroy(&graph);

#undef NAME
#undef REP

    /* Benchmark a much denser Erdős-Rényi graph as well. */

#define NAME "GNM(3000,30000)"
#define REP 1

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 3000, 30000, IGRAPH_DIRECTED, IGRAPH_LOOPS);

    BENCH("13 Betweenness, unweighted, " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL), REP)
    );
    BENCH("14 Betweenness, unweighted, " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL), REP)
    );

    rand_weight_vec(&weight, &graph);

    BENCH("15 Betweenness, weighted,   " NAME ", directed, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, &weight), REP)
    );
    BENCH("16 Betweenness, weighted,   " NAME ", undirected, " TOSTR(REP) "x",
          REPEAT(igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, &weight), REP)
    );

    igraph_destroy(&graph);

    igraph_vector_destroy(&weight);
    igraph_vector_destroy(&betweenness);

    return 0;
}
