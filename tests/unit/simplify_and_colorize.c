
#include <igraph.h>

#include "test_utilities.inc"

#define SIMPLIFY_PRINT_DESTROY(name) \
    printf(name "\n"); \
    igraph_simplify_and_colorize(&graph, &res, &vcol, &ecol); \
    print_graph(&res); \
    print_vector_int(&vcol); \
    print_vector_int(&ecol); \
    printf("\n"); \
    igraph_destroy(&res); \
    igraph_destroy(&graph);

int main() {
    igraph_t graph, res;
    igraph_vector_int_t vcol, ecol;

    igraph_vector_int_init(&vcol, 0);
    igraph_vector_int_init(&ecol, 0);

    /* null graph */
    igraph_empty(&graph, 0, 0);
    SIMPLIFY_PRINT_DESTROY("K0");

    /* singleton graph */
    igraph_empty(&graph, 1, 0);
    SIMPLIFY_PRINT_DESTROY("K1");

    /* 4-cycle-graph */
    igraph_ring(&graph, 4, 0, 0, 1);
    SIMPLIFY_PRINT_DESTROY("C4");

    /* both multi-edges and self loops */
    igraph_small(&graph, 2, 0,
                 0, 1, 0, 1, 1, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Undirected graph 1");

    /* parallel edges specified with different vertex orderings */
    igraph_small(&graph, 3, 0,
                 0, 1, 1, 2, 2, 0, 2, 2, 2, 2, 2, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Undirected graph 2");

    /* directed version of the same as above */
    igraph_small(&graph, 3, 1,
                 0, 1, 1, 2, 2, 0, 2, 2, 2, 2, 2, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Directed graph 1");

    /* isolated vertices */
    igraph_small(&graph, 4, 1,
                 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Directed graph 2");

    igraph_vector_int_destroy(&vcol);
    igraph_vector_int_destroy(&ecol);

    VERIFY_FINALLY_STACK();

    return 0;
}
