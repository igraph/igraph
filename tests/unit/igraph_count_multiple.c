
#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_vector_t counts;

    igraph_vector_init(&counts, 0);

    /* undirected case */
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED,
                 0,1, 0,1, 0,1, 0,0, 1,1, 1,1,
                 -1);

    igraph_count_multiple(&graph, &counts, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    print_vector_round(&counts);

    igraph_destroy(&graph);

    /* directed case */
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,1, 0,1, 0,1, 0,0, 1,1, 1,1,
                 -1);

    igraph_count_multiple(&graph, &counts, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    print_vector_round(&counts);

    igraph_destroy(&graph);

    igraph_vector_destroy(&counts);

    VERIFY_FINALLY_STACK();

    return 0;
}
