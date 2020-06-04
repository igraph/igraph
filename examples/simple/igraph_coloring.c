
#include <igraph.h>
#include <assert.h>

int main() {
    igraph_t graph;
    igraph_vector_int_t colors;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 10000, /* directed = */ 0, /* loops = */ 0);

    igraph_vector_int_init(&colors, 0);

    igraph_vertex_coloring_greedy(&graph, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);

    /* verify that the colouring is valid: */
    {
        long i;
        long no_of_edges = igraph_ecount(&graph);
        for (i = 0; i < no_of_edges; ++i) {
            assert( VECTOR(colors)[ IGRAPH_FROM(&graph, i) ] != VECTOR(colors)[ IGRAPH_TO(&graph, i) ]  );
        }
    }

    igraph_vector_int_destroy(&colors);
    igraph_destroy(&graph);

    return 0;
}
