
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"

void sort_and_print_vector(igraph_vector_t *v) {
    igraph_vector_sort(v);
    print_vector_round(v);
}

int main() {
    igraph_t graph;
    igraph_vector_t bridges;

    igraph_vector_init(&bridges, 0);

    igraph_empty(&graph, 0, /* directed = */ 0);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_empty(&graph, 1, /* directed = */ 0);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_small(&graph, /* num_nodes = */ 7, /* directed = */ 0,
                 0, 1, 1, 2, 0, 2, 0, 3, 3, 4, 4, 5, 3, 5, 4, 6, -1);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_small(&graph, /* num_nodes = */ 3, /* directed = */ 0,
                 0, 1, 0, 1, 1, 2, 2, 2, -1);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_vector_destroy(&bridges);

    VERIFY_FINALLY_STACK();

    return 0;
}
