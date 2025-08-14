
#include <igraph.h>
#include "igraph_interface.h"
#include "igraph_matrix.h"
#include "igraph_spatial.h"
#include "test_utilities.h"

int main(void) {

    igraph_t graph;

    igraph_real_t points[] = {
0.474217, 0.0314797, 0.208089, 0.439308, 0.967367, 0.530466, \
0.177005, 0.426713, 0.568462, 0.57507, 0.441834, 0.284514, 0.479224, \
0.817988, 0.720209, 0.225744, 0.204941, 0.44297, 0.285318, 0.912984, \
0.831097, 0.0176603, 0.827154, 0.472702, 0.173059, 0.561858, \
0.156276, 0.88019, 0.65935, 0.538207, 0.570379, 0.518081, 0.900553, \
0.656416, 0.726631, 0.863709, 0.380264, 0.287159, 0.31098, 0.230773, \
0.243089, 0.164584, 0.967974, 0.524992, 0.726605, 0.0724703, \
0.739752, 0.447069, 0.0443581, 0.444839
    };

    igraph_matrix_t point_mat;

    igraph_matrix_init_array(&point_mat, points, 25, 2, false);

    igraph_lune_beta_skeleton(&graph, &point_mat, 2);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    igraph_lune_beta_skeleton(&graph, &point_mat, 1);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    igraph_lune_beta_skeleton(&graph, &point_mat, 0.5);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    igraph_circle_beta_skeleton(&graph, &point_mat, 1.1);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    igraph_matrix_destroy(&point_mat);
    VERIFY_FINALLY_STACK();
    return 0;
}
