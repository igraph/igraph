#include <igraph.h>

#include "bench.h"

int main() {
    igraph_t graph;
    igraph_vector_t betweenness;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&betweenness, 0);

    /* Kautz and De Bruijn graphs are connected, therefore there should not be a dramatic difference
     * in the performance of the directed and undirected calculations. */

    igraph_kautz(&graph, 4, 5);
    BENCH(" 1 Betweenness, Kautz(4,5), directed",
          igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL));
    BENCH(" 2 Betweenness, Kautz(4,5), undirected",
          igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL));
    igraph_destroy(&graph);

    igraph_de_bruijn(&graph, 6, 5);
    BENCH(" 3 Betweenness, DeBruijn(6,5), directed",
          igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL));
    igraph_destroy(&graph);

    {
        igraph_vector_t dims;

        igraph_vector_init_int_end(&dims, -1, 8, 8, 8, 8, -1);
        igraph_lattice(&graph, &dims, 1, IGRAPH_UNDIRECTED, /* mutual */ 0, /* circular */ 0);
        BENCH(" 4 Betweenness, Grid(8,8,8,8), undirected",
              igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL));
        igraph_destroy(&graph);
        igraph_vector_destroy(&dims);

        igraph_vector_init_int_end(&dims, -1, 10, 10, 10, 10, -1);
        igraph_lattice(&graph, &dims, 1, IGRAPH_UNDIRECTED, /* mutual */ 0, /* circular */ 0);
        BENCH(" 5 Betweenness, Grid(10,10,10,10), cutoff 5",
              igraph_betweenness_cutoff(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL, 5));
        BENCH(" 6 Betweenness, Grid(10,10,10,10), cutoff 8",
              igraph_betweenness_cutoff(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL, 8));
        igraph_destroy(&graph);
        igraph_vector_destroy(&dims);
    }

    igraph_barabasi_game(&graph, 8000, 1, 1, NULL, 1, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 7 Betweenness, Barabasi n=8000 m=1, undirected",
          igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL));
    BENCH(" 8 Betweenness, Barabasi n=8000 m=1, undirected, cutoff 6",
          igraph_betweenness_cutoff(&graph, &betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL, 6));
    igraph_destroy(&graph);

    igraph_barabasi_game(&graph, 30000, 1, 5, NULL, 1, 0, IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 9 Betweenness, Barabasi n=30000 m=5, directed",
          igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL));
    BENCH("10 Betweenness, Barabasi n=30000 m=5, directed, cutoff 5",
          igraph_betweenness_cutoff(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED, NULL, 5));
    igraph_destroy(&graph);

    igraph_vector_destroy(&betweenness);

    return 0;
}
