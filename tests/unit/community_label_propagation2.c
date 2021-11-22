
#include <igraph.h>

#include "test_utilities.inc"

/* Test case for bug #1852 */

int main() {
    igraph_t graph;
    igraph_vector_t membership, initial_labels;
    igraph_real_t modularity;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Undirected graph with unlabelled components */

    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
                 1,2, -1);

    igraph_vector_init(&membership, 0);
    igraph_vector_init(&initial_labels, igraph_vcount(&graph));
    VECTOR(initial_labels)[0] = 1;
    VECTOR(initial_labels)[1] = -1;
    VECTOR(initial_labels)[2] = -1;
    VECTOR(initial_labels)[3] = -1;

    igraph_community_label_propagation(&graph, &membership, NULL, &initial_labels, NULL, &modularity);
    print_vector(&membership);

    igraph_destroy(&graph);

    /* Directed graph with unlabelled nodes not reachable from any labelled ones. */

    igraph_small(&graph, 8, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 3, 1, 2, 4, 4, 5, 5, 2, 4, 6,
                 -1);

    igraph_vector_resize(&initial_labels, igraph_vcount(&graph));
    igraph_vector_null(&initial_labels);
    VECTOR(initial_labels)[0] = -1;
    VECTOR(initial_labels)[1] = -1;
    VECTOR(initial_labels)[2] = 1;
    VECTOR(initial_labels)[3] = -1;
    VECTOR(initial_labels)[4] = 2;
    VECTOR(initial_labels)[6] = -1;
    VECTOR(initial_labels)[7] = -1;

    igraph_community_label_propagation(&graph, &membership, NULL, &initial_labels, NULL, &modularity);
    print_vector(&membership);

    igraph_destroy(&graph);

    /* None of the nodes are labelled initially */

    igraph_full(&graph, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_resize(&initial_labels, igraph_vcount(&graph));
    igraph_vector_fill(&initial_labels, -1);

    igraph_community_label_propagation(&graph, &membership, NULL, &initial_labels, NULL, &modularity);
    print_vector(&membership);

    igraph_destroy(&graph);

    igraph_vector_destroy(&initial_labels);
    igraph_vector_destroy(&membership);

    return 0;
}
