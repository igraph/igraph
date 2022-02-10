#include <igraph.h>
#include <stdio.h>

int main(){
    igraph_t g;
    igraph_integer_t nodes = 100000, A = 0, power = 1, m = 1;
    
    igraph_rng_seed(igraph_rng_default(), 42);
    printf("Barabasi-Albert network with 100000 nodes\n\n");

    for (int i = 0; i < 5; i++) {
        igraph_real_t assortativity;

        /* Generate undirected graph with 100000 nodes */
        igraph_barabasi_game(&g, nodes, power, m, NULL, /* outpref */ 0, A, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, /* start from */ NULL);

        igraph_assortativity_degree(&g, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
        printf("Assortativity before rewiring = %g\n", assortativity);
        
        /* Rewire graph */
        igraph_rewire(&g, 10 * igraph_ecount(&g), IGRAPH_REWIRING_SIMPLE);
        
        igraph_assortativity_degree(&g, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
        printf("Assortativity after rewiring = %g\n\n", assortativity);

        igraph_destroy(&g);
    }
}
