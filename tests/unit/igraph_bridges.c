
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"

void sort_and_print_vector(igraph_vector_t *v) {
    long int i, n = igraph_vector_size(v);
    igraph_vector_sort(v);
    for (i = 0; i < n; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

int main() {
    igraph_t graph;
    igraph_vector_t bridges;

    igraph_vector_init(&bridges, 0);

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
