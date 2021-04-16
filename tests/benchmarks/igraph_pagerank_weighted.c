#include <igraph.h>

#include "bench.h"

void rand_weight_vec(igraph_vector_t *vec, const igraph_t *graph) {
    long i, n = igraph_ecount(graph);
    igraph_vector_resize(vec, n);
    for (i=0; i < n; ++i) {
        VECTOR(*vec)[i] = RNG_UNIF(1, 10);
    }
}

int main() {
    igraph_t graph;
    igraph_vector_t res, weights;
    igraph_arpack_options_t arpack_opts;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_arpack_options_init(&arpack_opts);

    igraph_vector_init(&res, 0);
    igraph_vector_init(&weights, 0);

    igraph_barabasi_game(&graph, 100000, 1, 4, NULL, 1, 0, IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weight_vec(&weights, &graph);
    BENCH(" 1 PageRank weighted, Barabasi n=100000 m=4, PRPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 10)
    );
    BENCH(" 2 PageRank weighted, Barabasi n=100000 m=4, ARPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 10)
    );
    igraph_destroy(&graph);

    igraph_barabasi_game(&graph, 100000, 1, 10, NULL, 1, 0, IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weight_vec(&weights, &graph);
    BENCH(" 3 PageRank weighted, Barabasi n=100000 m=10, PRPACK, 5x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 5)
    );
    BENCH(" 4 PageRank weighted, Barabasi n=100000 m=10, ARPACK, 5x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 5)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100, 1000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    rand_weight_vec(&weights, &graph);
    BENCH(" 5 PageRank weighted, GNM(100,1000), PRPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 1000)
    );
    BENCH(" 6 PageRank weighted, GNM(100,1000), ARPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 1000)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 200, 4000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    rand_weight_vec(&weights, &graph);
    BENCH(" 7 PageRank weighted, GNM(200,4000), PRPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 1000)
    );
    BENCH(" 8 PageRank weighted, GNM(200,4000), ARPACK, 1000x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 1000)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 10000, 20000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    rand_weight_vec(&weights, &graph);
    BENCH(" 9 PageRank weighted, GNM(10000,20000), PRPACK, 100x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 100)
    );
    BENCH("10 PageRank weighted, GNM(10000,20000), ARPACK, 100x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 100)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100000, 100000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    rand_weight_vec(&weights, &graph);
    BENCH("11 PageRank weighted, GNM(100000,100000), PRPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 10)
    );
    BENCH("12 PageRank weighted, GNM(100000,100000), ARPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 10)
    );
    igraph_destroy(&graph);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100000, 500000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    rand_weight_vec(&weights, &graph);
    BENCH("13 PageRank weighted, GNM(100000,500000), PRPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 10)
    );
    BENCH("14 PageRank weighted, GNM(100000,500000), ARPACK, 10x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 10)
    );
    igraph_destroy(&graph);

    igraph_kautz(&graph, 6, 6);
    rand_weight_vec(&weights, &graph);
    BENCH("13 PageRank weighted, Kautz(6,6), PRPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 1)
    );
    BENCH("14 PageRank weighted, Kautz(6,6), ARPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 1)
    );
    igraph_destroy(&graph);

    igraph_de_bruijn(&graph, 7, 7);
    rand_weight_vec(&weights, &graph);
    BENCH("13 PageRank weighted, DeBruijn(7,7), PRPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, NULL), 1)
    );
    BENCH("14 PageRank weighted, DeBruijn(7,7), ARPACK, 1x",
          REPEAT(igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_ARPACK, &res, NULL, igraph_vss_all(), IGRAPH_DIRECTED, 0.85, &weights, &arpack_opts), 1)
    );
    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&res);

    return 0;
}
