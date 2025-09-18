#include <igraph.h>

int main(void) {
    igraph_int_t num_vertices = 1000;
    igraph_int_t num_edges = 1000;
    igraph_real_t diameter, mean_degree;
    igraph_t graph;

    /* Initialize the library. */
    igraph_setup();

    /* Ensure identical results across runs. */
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_erdos_renyi_game_gnm(
            &graph, num_vertices, num_edges,
            IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    igraph_diameter(
        &graph, /* weights = */ NULL,
        &diameter,
        /* from = */ NULL, /* to = */ NULL,
        /* vertex_path = */ NULL, /* edge_path = */ NULL,
        IGRAPH_UNDIRECTED, /* unconn= */ true);

    igraph_mean_degree(&graph, &mean_degree, IGRAPH_LOOPS);
    printf("Diameter of a random graph with average degree %g: %g\n",
           mean_degree, diameter);

    igraph_destroy(&graph);

    return 0;
}
