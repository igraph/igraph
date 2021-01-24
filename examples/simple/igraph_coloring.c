
#include <igraph.h>

igraph_long_t main() {
    igraph_t graph;
    igraph_vector_long_t colors;

    /* Setting a seed makes the result of erdos_renyi_game deterministic. */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* IGRAPH_UNDIRECTED and IGRAPH_NO_LOOPS are both equivalent to 0/FALSE, but
       communicate intent better in this context. */
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    /* As with all igraph functions, the vector in which the result is returned must
       be initialized in advance. */
    igraph_vector_long_init(&colors, 0);
    igraph_vertex_coloring_greedy(&graph, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);

    /* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
    {
        igraph_long_t i;
        /* Store the edge count to avoid the overhead from igraph_ecount in the for loop. */
        igraph_long_t no_of_edges = igraph_ecount(&graph);
        for (i = 0; i < no_of_edges; ++i) {
            if ( VECTOR(colors)[ IGRAPH_FROM(&graph, i) ] == VECTOR(colors)[ IGRAPH_TO(&graph, i) ]  ) {
                printf("Inconsistent coloring! Vertices %ld and %ld are adjacent but have the same color.\n",
                       (igraph_long_t) IGRAPH_FROM(&graph, i), (igraph_long_t) IGRAPH_TO(&graph, i));
                abort();
            }
        }
    }

    /* Destroy data structure when we are done. */
    igraph_vector_long_destroy(&colors);
    igraph_destroy(&graph);

    return 0;
}
