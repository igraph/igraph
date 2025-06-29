
#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_t vertices, edges;

    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);

    printf("Null graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_OUT);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    igraph_destroy(&g);

    printf("Several isolated vertices:\n");
    igraph_empty(&g, 3, IGRAPH_UNDIRECTED);
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_OUT);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    igraph_destroy(&g);

    printf("Isolated vertices with self-loops:\n");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED,
                 1,1, 1,1, 2,2,
                 -1);
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_OUT);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    igraph_destroy(&g);

    printf("Small directed graph:\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,4, 2,0,
                 -1);
    printf("OUT\n");
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_OUT);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    printf("IN\n");
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_IN);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    printf("ALL\n");
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_IN);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    igraph_destroy(&g);

    printf("Small undirected multigraph:\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 1,2, 3,4, 3,4, 3,4,
                 -1);
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_ALL);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    igraph_destroy(&g);

    printf("Directed acyclic graph:\n");
    igraph_small(&g, 8, IGRAPH_DIRECTED,
                 0,1, 1,2, 0,3, 3,2, 2,4, 5,6,
                 -1);
    printf("OUT\n");
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_OUT);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    printf("IN\n");
    igraph_find_cycle(&g, &vertices, &edges, IGRAPH_IN);
    print_vector_int(&vertices);
    print_vector_int(&edges);
    igraph_destroy(&g);

    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&vertices);

    VERIFY_FINALLY_STACK();

    return 0;
}
