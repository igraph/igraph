
#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_vector_t res;
    igraph_vector_t weights;

    igraph_vector_init(&res, 0);    

    /* Path graph */
    igraph_ring(&graph, 7, IGRAPH_DIRECTED, 0, /* circular */ 0);

    printf("Unweighted undirected:\n");
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ NULL, /* normalized= */ 1);
    print_vector(&res);
    printf("Unweighted directed:\n");
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_OUT, /* weights= */ NULL, /* normalized= */ 1);
    print_vector(&res);

    printf("Unweighted undirected, cutoff=0:\n");
    igraph_harmonic_centrality_cutoff(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ NULL, /* normalized= */ 1, 0);
    print_vector(&res);
    printf("Unweighted undirected, cutoff=1:\n");
    igraph_harmonic_centrality_cutoff(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ NULL, /* normalized= */ 1, 1);
    print_vector(&res);

    igraph_vector_init(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1.0);

    printf("Unit-weighted undirected:\n");
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1);    
    print_vector(&res);
    printf("Unit-weighted directed:\n");
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_OUT, /* weights= */ &weights, /* normalized= */ 1);
    print_vector(&res);

    printf("Unit-weighted undirected, cutoff=0:\n");
    igraph_harmonic_centrality_cutoff(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1, 0);
    print_vector(&res);
    printf("Unit-weighted undirected, cutoff=1:\n");
    igraph_harmonic_centrality_cutoff(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1, 1);
    print_vector(&res);

    igraph_vector_destroy(&weights);

    igraph_vector_init_seq(&weights, 1, igraph_ecount(&graph));
    printf("Weighted undirected:\n");
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1);
    print_vector(&res);
    printf("Weighted directed:\n");
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_OUT, /* weights= */ &weights, /* normalized= */ 1);
    print_vector(&res);
    igraph_vector_destroy(&weights);

    igraph_destroy(&graph);

    /* Graphs with no edges */

    igraph_vector_init(&weights, 0);

    /* Null graph */

    printf("Null graph:\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ NULL, /* normalized= */ 1);
    print_vector(&res);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1);
    print_vector(&res);
    igraph_destroy(&graph);

    /* Singleton graph */

    printf("Singleton graph:\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ NULL, /* normalized= */ 1);
    print_vector(&res);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1);
    print_vector(&res);
    igraph_destroy(&graph);

    /* Empty graph with two vertices */

    printf("Empty graph with two vertices:\n");
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ NULL, /* normalized= */ 1);
    print_vector(&res);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1);
    print_vector(&res);
    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);

    /* Graph with multiple connected components and isolated vertices */

    printf("Multiple components, unweighted:\n");
    igraph_small(&graph, 20, IGRAPH_UNDIRECTED,
                 1, 2, 2, 3, 1, 3, 4, 5, 5, 6, 6, 7, 5, 8, 8, 9, 6, 10, 10, 11, 7, 11,
                 9, 13, 9, 14, 9, 15, 9, 16, 13, 14, 13, 15, 13, 16, 14, 15, 14, 16,
                 15, 16, 17, 18, 4, 19, 4, 20,
                 -1);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ NULL, /* normalized= */ 1);
    print_vector(&res);

    printf("Multiple components, constant weight vector:\n");
    igraph_vector_init(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1.0 / 7);
    igraph_harmonic_centrality(&graph, &res, igraph_vss_all(), IGRAPH_ALL, /* weights= */ &weights, /* normalized= */ 1);
    print_vector(&res);
    igraph_vector_destroy(&weights);

    igraph_destroy(&graph);

    igraph_vector_destroy(&res);

    VERIFY_FINALLY_STACK();

    return 0;
}
