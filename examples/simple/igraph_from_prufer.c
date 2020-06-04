
#include <igraph.h>
#include <stdio.h>
#include <assert.h>

void print_edges(const igraph_t *graph) {
    long ecount = igraph_ecount(graph);
    long i;

    for (i = 0; i < ecount; ++i) {
        printf("%d %d\n", IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
    }
    printf("\n");
}


int main() {
    igraph_t graph;
    igraph_integer_t prufer1[] = {2, 3, 2, 3};
    igraph_integer_t prufer2[] = {0, 2, 4, 1, 1, 0};
    igraph_integer_t prufer3[] = {};
    igraph_vector_int_t prufer;
    igraph_bool_t tree;

    igraph_vector_int_view(&prufer, prufer1, sizeof(prufer1) / sizeof(igraph_integer_t));
    igraph_from_prufer(&graph, &prufer);
    igraph_is_tree(&graph, &tree, NULL, IGRAPH_ALL);
    assert(tree);
    print_edges(&graph);
    igraph_destroy(&graph);

    igraph_vector_int_view(&prufer, prufer2, sizeof(prufer2) / sizeof(igraph_integer_t));
    igraph_from_prufer(&graph, &prufer);
    igraph_is_tree(&graph, &tree, NULL, IGRAPH_ALL);
    assert(tree);
    print_edges(&graph);
    igraph_destroy(&graph);

    igraph_vector_int_view(&prufer, prufer3, sizeof(prufer3) / sizeof(igraph_integer_t));
    igraph_from_prufer(&graph, &prufer);
    igraph_is_tree(&graph, &tree, NULL, IGRAPH_ALL);
    assert(tree);
    print_edges(&graph);
    igraph_destroy(&graph);

    return 0;
}
