#include <igraph.h>
#include <stdio.h>

#include "../unit/test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_matrix_t layout;

    if (igraph_empty(&graph, 2, 0)) {
        return 1;
    }

    if (igraph_add_edge(&graph, 0, 1)) {
        return 2;
    }

    if (igraph_matrix_init(&layout, 0, 0)) {
        return 3;
    }

    if (igraph_layout_kamada_kawai_3d(&graph, &layout, 0, 200, 0, 2, 0, 0, 0, 0, 0, 0, 0)) {
        return 4;
    }

    igraph_matrix_printf(&layout, "%.2f");

    igraph_matrix_destroy(&layout);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
