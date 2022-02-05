#include <igraph.h>
#include <stdio.h>

int main(){
    igraph_t g;
    igraph_integer_t nodes = 100000, A = 0, power = 1, m = 1;
    igraph_vector_t degree, weight;
    igraph_vector_init(&degree, 0);
    
    igraph_rng_seed(igraph_rng_default(), 42);
    printf("Barabasi-Albert network with 100000 nodes\n\n");

    for (int i = 0; i < 5; i++) {
        /* Generate undirected graph with 100000 nodes */
        igraph_barabasi_game(&g, nodes, power, m, NULL, IGRAPH_UNDIRECTED, A, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, /* start from */ NULL);

        igraph_integer_t edge_count = igraph_ecount(&g);
        igraph_vector_init(&weight, edge_count);
        igraph_vector_fill(&weight, /* weight of each edge = 1 */ 1);

        /* Get degree sequence */
        igraph_strength(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /* consider self-loops */ IGRAPH_LOOPS, &weight);
       
        igraph_real_t assortativity;
        igraph_assortativity(&g, &degree, NULL, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
        printf("Assortativity before rewiring = %g\n", assortativity);
        
        /* Rewire graph */
        igraph_rewire(&g, 10 * edge_count, IGRAPH_REWIRING_SIMPLE);
        
        igraph_assortativity(&g, &degree, NULL, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
        printf("Assortativity after rewiring = %g\n\n", assortativity);
    }

    igraph_vector_destroy(&weight);
    igraph_vector_destroy(&degree);
    igraph_destroy(&g);
}
