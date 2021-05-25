
#include <igraph.h>
#include "test_utilities.inc"


#define PRINT_DESTROY(name) \
    printf(name "\n"); \
    print_graph_canon(&graph); \
    igraph_destroy(&graph); \
    printf("\n");


int main() {
    igraph_t graph;

    igraph_tree(&graph, 0, 1, IGRAPH_TREE_UNDIRECTED);
    PRINT_DESTROY("Null graph");

    igraph_tree(&graph, 0, 1, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Directed null graph");

    igraph_tree(&graph, 1, 1, IGRAPH_TREE_UNDIRECTED);
    PRINT_DESTROY("Singleton graph");

    igraph_tree(&graph, 3, 1, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Path graph");

    igraph_tree(&graph, 3, 2, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Binary out-tree, n=3");

    igraph_tree(&graph, 3, 2, IGRAPH_TREE_IN);
    PRINT_DESTROY("Binary in-tree, n=3");

    igraph_tree(&graph, 14, 3, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Ternary out-tree, n=14");

    VERIFY_FINALLY_STACK();

    return 0;
}
