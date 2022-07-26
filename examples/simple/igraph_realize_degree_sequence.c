#include <igraph.h>
#include <stdio.h>

int main(){
    igraph_t g1, g2, g3;
    igraph_integer_t nodes = 500, A = 0, power = 1, m = 1;
    igraph_real_t assortativity;
    
    igraph_rng_seed(igraph_rng_default(), 42);
    printf("Demonstration of difference in assortativities of graphs with the same degree sequence but different linkages:\n\nInitial graph based on the Barabasi-Albert model with %" IGRAPH_PRId " nodes.\n", nodes);
    
    /* Graph 1 generated by a randomized graph generator */
    igraph_barabasi_game(&g1, nodes, power, m, NULL, /* outpref */ 0, A, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, /* start from */ NULL);

    igraph_vector_int_t degree;
    igraph_vector_int_init(&degree, nodes);
    igraph_degree(&g1, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    /* Measuring assortativity of the first graph */
    igraph_assortativity_degree(&g1, &assortativity, IGRAPH_UNDIRECTED);
    printf("Assortativity of initial graph = %g\n\n", assortativity);
    igraph_destroy(&g1);

    /* Graph 2 (with the same degree sequence) generated by selecting vertices with the smallest degree first */
    igraph_realize_degree_sequence(&g2, &degree, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
    igraph_assortativity_degree(&g2, &assortativity, IGRAPH_UNDIRECTED);
    printf("Assortativity after choosing vertices with the smallest degrees first = %g\n\n", assortativity);
    igraph_destroy(&g2);

    /* Graph 3 (with the same degree sequence) generated by selecting vertices with the largest degree first */
    igraph_realize_degree_sequence(&g3, &degree, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);
    igraph_assortativity_degree(&g3, &assortativity, IGRAPH_UNDIRECTED);
    printf("Assortativity after choosing vertices with the largest degrees first = %g\n", assortativity);
    igraph_destroy(&g3);
    igraph_vector_int_destroy(&degree);

    return 0;
}
