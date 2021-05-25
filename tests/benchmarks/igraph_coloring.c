
#include <igraph.h>

#include "bench.h"

int main() {
    igraph_t g;
    igraph_vector_int_t colors;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_vector_int_init(&colors, 0);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 50000, 1000000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    BENCH(" 1 Vertex coloring random graph with 50,000 vertices and 2,000,000 edges",
          igraph_vertex_coloring_greedy(&g, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS)
         );
    igraph_destroy(&g);

    igraph_barabasi_game(&g, 100000, 1, 15, NULL, 0, 0, 0, IGRAPH_BARABASI_PSUMTREE, NULL);
    BENCH(" 2 Vertex coloring pref. attach. graph n=100,000 m=15",
          igraph_vertex_coloring_greedy(&g, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS)
         );
    igraph_destroy(&g);

    igraph_vector_int_destroy(&colors);


    return 0;
}
