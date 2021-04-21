#include <igraph.h>

#include "bench.h"

int main() {
    igraph_t graph;
    igraph_vector_t res;
    igraph_arpack_options_t arpack_opts;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_arpack_options_init(&arpack_opts);

    igraph_vector_init(&res, 0);

    igraph_barabasi_game(&graph, 100000, 1, 4, NULL, 1, 0, IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 1 PageRank, Barabasi n=100000 m=4, PRPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 10)
    );
    BENCH(" 2 PageRank, Barabasi n=100000 m=4, ARPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 10)
    );
    igraph_destroy(&graph);

    igraph_barabasi_game(&graph, 100000, 1, 10, NULL, 1, 0, IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 3 PageRank, Barabasi n=100000 m=10, PRPACK, 5x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 5)
    );
    BENCH(" 4 PageRank, Barabasi n=100000 m=10, ARPACK, 5x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 5)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100, 1000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    BENCH(" 5 PageRank, GNM(100,1000), PRPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 1000)
    );
    BENCH(" 6 PageRank, GNM(100,1000), ARPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 1000)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 200, 4000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    BENCH(" 7 PageRank, GNM(200,4000), PRPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 1000)
    );
    BENCH(" 8 PageRank, GNM(200,4000), ARPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 1000)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 10000, 20000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    BENCH(" 9 PageRank, GNM(10000,20000), PRPACK, 100x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 100)
    );
    BENCH("10 PageRank, GNM(10000,20000), ARPACK, 100x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 100)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100000, 100000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    BENCH("11 PageRank, GNM(100000,100000), PRPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 10)
    );
    BENCH("12 PageRank, GNM(100000,100000), ARPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 10)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100000, 500000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    BENCH("13 PageRank, GNM(100000,500000), PRPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 10)
    );
    BENCH("14 PageRank, GNM(100000,500000), ARPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 10)
    );
    igraph_destroy(&graph);

    igraph_kautz(&graph, 6, 6);
    BENCH("13 PageRank, Kautz(6,6), PRPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 1)
    );
    BENCH("14 PageRank, Kautz(6,6), ARPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 1)
    );
    igraph_destroy(&graph);

    igraph_de_bruijn(&graph, 7, 7);
    BENCH("13 PageRank, DeBruijn(7,7), PRPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, NULL), 1)
    );
    BENCH("14 PageRank, DeBruijn(7,7), ARPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, NULL, &arpack_opts), 1)
    );
    igraph_destroy(&graph);

    igraph_vector_destroy(&res);

    return 0;
}
