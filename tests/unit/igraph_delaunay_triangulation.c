#include "igraph_interface.h"
#include "igraph_matrix.h"
#include "test_utilities.h"

#include <igraph.h>
#include <stdbool.h>

int main(void) {
    igraph_matrix_t points;

    igraph_t g;

    igraph_real_t points_raw[10] = {
        1, 1,
        4, 3,
        2, 3,
        6, 4,
        5, 3,
    };

    igraph_matrix_init_array(&points, &points_raw[0], 5, 2, false);
    igraph_delaunay_triangulation(&g, &points);

    print_graph_canon(&g);
    igraph_matrix_destroy(&points);
    igraph_destroy(&g);
}
