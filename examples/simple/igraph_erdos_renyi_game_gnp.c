
#include <igraph.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_t component_sizes;

    igraph_rng_seed(igraph_rng_default(), 42); /* make program deterministic */

    /* Sample a graph from the Erdős-Rényi G(n,p) model */

    igraph_erdos_renyi_game_gnp(
        &graph, /* n= */ 100, /* p= */ 0.01,
        IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS
    );

    /* Compute the fraction of vertices contained within the largest connected component */

    igraph_vector_int_init(&component_sizes, 0);
    igraph_connected_components(&graph, NULL, &component_sizes, NULL, IGRAPH_STRONG);

    printf(
        "Fraction of vertices in giant component: %g\n",
        ((double) igraph_vector_int_max(&component_sizes)) / igraph_vcount(&graph)
    );

    /* Clean up data structures when no longer needed */

    igraph_vector_int_destroy(&component_sizes);
    igraph_destroy(&graph);

    return 0;
}
