#include <igraph.h>
#include <stdio.h>

#include "../unit/test_utilities.h"

void snap_to_zero(igraph_real_t* value) {
    if (fabs(*value) < 1e-5) {
        *value = 0.0;
    }
}

int main(void) {
    igraph_t graph;
    igraph_matrix_t layout;
    igraph_integer_t i;

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

    /* Snap numbers close to zero in the layout; there are false failures on
     * MinGW if we don't do so */
    for (i = 0; i < 2; i++) {
        snap_to_zero(&MATRIX(layout, i, 0));
        snap_to_zero(&MATRIX(layout, i, 1));
        snap_to_zero(&MATRIX(layout, i, 2));
    }
    print_matrix_format(&layout, stdout, "%.2f");

    igraph_matrix_destroy(&layout);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
