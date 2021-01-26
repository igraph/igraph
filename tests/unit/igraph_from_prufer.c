
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"


int main() {
    igraph_t graph;
    igraph_integer_t prufer1[] = {2, 3, 2, 3};
    igraph_integer_t prufer2[] = {0, 2, 4, 1, 1, 0};
    igraph_vector_int_t prufer;
    igraph_bool_t tree;

    igraph_vector_int_view(&prufer, prufer1, sizeof(prufer1) / sizeof(igraph_integer_t));
    igraph_from_prufer(&graph, &prufer);
    igraph_is_tree(&graph, &tree, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(tree);
    print_graph(&graph);
    igraph_destroy(&graph);

    igraph_vector_int_view(&prufer, prufer2, sizeof(prufer2) / sizeof(igraph_integer_t));
    igraph_from_prufer(&graph, &prufer);
    igraph_is_tree(&graph, &tree, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(tree);
    print_graph(&graph);
    igraph_destroy(&graph);

    /* For a zero-length array, we cannot use the same pattern as above because
       standard C does not allow raw zero-length arrays. */
    igraph_vector_int_init(&prufer, 0);
    igraph_from_prufer(&graph, &prufer);
    igraph_is_tree(&graph, &tree, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(tree);
    print_graph(&graph);
    igraph_destroy(&graph);
    igraph_vector_int_destroy(&prufer);

    VERIFY_FINALLY_STACK();

    return 0;
}
