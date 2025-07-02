#include <igraph.h>
#include "test_utilities.h"

int main(void) {

    igraph_t graph;

    igraph_real_t pointArray[8] = {
        0, 0,
        1, 1,
        0, 1,
        1, 0
    };

    igraph_matrix_t points;



    igraph_matrix_init_array(&points, &pointArray[0], 4, 2, IGRAPH_ROW_MAJOR);

    igraph_nearest_neighbor_graph(&graph, &points, IGRAPH_METRIC_L2, 10, 10);

    print_graph_canon(&graph);

    VERIFY_FINALLY_STACK();
    return 1;
}
