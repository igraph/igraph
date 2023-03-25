#include <igraph.h>
#include <stdio.h>

int main(void){
    igraph_t g;
    igraph_integer_t vcount = 1000;
    igraph_real_t pf = 0.2;

    /* Seed random number generator to ensure reproducibility. */
    igraph_rng_seed(igraph_rng_default(), 42);

    printf("Forest fire model network with %" IGRAPH_PRId " vertices and %g forward burning probability.\n\n",
           vcount, pf);

    for (int i = 0; i < 5; i++) {
        igraph_real_t assortativity;

        /* Generate graph from the forest fire model. */
        igraph_forest_fire_game(&g, vcount, pf, 1.0, 1, IGRAPH_UNDIRECTED);

        /* Compute assortativity. */
        igraph_assortativity_degree(&g, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
        printf("Assortativity before rewiring = %g\n", assortativity);

        /* Randomize the graph while preserving the degrees. */
        igraph_rewire(&g, 20 * igraph_ecount(&g), IGRAPH_REWIRING_SIMPLE);

        /* Re-compute assortativity. Did it change? */
        igraph_assortativity_degree(&g, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
        printf("Assortativity after rewiring = %g\n\n", assortativity);

        igraph_destroy(&g);
    }
}
