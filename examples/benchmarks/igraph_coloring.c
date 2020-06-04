
#include <igraph.h>

#include "bench.h"

int main() {
    igraph_t g;
    igraph_vector_int_t colors;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 30000, 300000, /* directed = */ 0, /* loops = */ 0);

    igraph_vector_int_init(&colors, 0);

    BENCH("1 Vertex coloring a random graph with 30,000 vertices and 300,000 edges.",
          igraph_vertex_coloring_greedy(&g, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS)
         );

    /* Use the result to prevent optimizing it away. */
    printf("Number of colors used: %d\n", igraph_vector_int_max(&colors));

    igraph_vector_int_destroy(&colors);
    igraph_destroy(&g);

    return 0;
}
