
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.h"

void sort_and_print_vector(igraph_vector_int_t *v) {
    igraph_vector_int_sort(v);
    print_vector_int(v);
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_t bridges;

    igraph_vector_int_init(&bridges, 0);

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_small(&graph, /* num_nodes = */ 7, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 0, 2, 0, 3, 3, 4, 4, 5, 3, 5, 4, 6, -1);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    /* Test with disconnected graph. */
    igraph_small(&graph, /* num_nodes = */ 16, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 1, 3, 4, 5, 5, 6, 4, 6, 4, 7, 7, 8, 4, 8, 9, 10, 10, 11,
                 11, 12, 9, 12, 9, 13, 13, 14, -1);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    /* Test with multi-edges and self-loops. */
    igraph_small(&graph, /* num_nodes = */ 3, IGRAPH_UNDIRECTED,
                 0, 1, 0, 1, 1, 2, 2, 2, -1);
    igraph_bridges(&graph, &bridges);
    sort_and_print_vector(&bridges);
    igraph_destroy(&graph);

    igraph_vector_int_destroy(&bridges);

    VERIFY_FINALLY_STACK();

    return 0;
}
