#include <igraph.h>
#include "igraph_interface.h"
#include "igraph_matrix.h"
#include "igraph_spatial.h"
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

    igraph_nearest_neighbor_graph(&graph, &points, IGRAPH_METRIC_L2, 1, 1);

    print_graph_canon(&graph);
    igraph_destroy(&graph);
    igraph_matrix_destroy(&points);
    VERIFY_FINALLY_STACK();
    return 1;
}
