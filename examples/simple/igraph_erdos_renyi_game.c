
#include <igraph.h>

int main() {
    igraph_t graph;
    igraph_vector_t component_sizes;

    igraph_rng_seed(igraph_rng_default(), 42); /* make program deterministic */

    /* Sample a graph from the Erdős-Rényi G(n,m) model */

    igraph_erdos_renyi_game(
                &graph, IGRAPH_ERDOS_RENYI_GNM,
                /* n= */ 100, /* m= */ 100,
                IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    /* Compute the fraction of vertices contained within the largest connected component */

    igraph_vector_init(&component_sizes, 0);
    igraph_clusters(&graph, NULL, &component_sizes, NULL, IGRAPH_STRONG);

    printf("Fraction of vertices in giant component: %g\n", (double) igraph_vector_max(&component_sizes) / igraph_vcount(&graph));

    /* Clean up data structures when no longer needed */

    igraph_vector_destroy(&component_sizes);
    igraph_destroy(&graph);

    return 0;
}
