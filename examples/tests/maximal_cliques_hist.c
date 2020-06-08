
#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_vector_t hist;

    igraph_small(&graph, 6, 0,
                 1, 2, 2, 3, 3, 4, 4, 5, 5, 2, 2, 4,
                 -1);

    igraph_vector_init(&hist, 0);

    igraph_maximal_cliques_hist(&graph, &hist, 0, 0);
    igraph_vector_print(&hist);

    igraph_vector_destroy(&hist);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
