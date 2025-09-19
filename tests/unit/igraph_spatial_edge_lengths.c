/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_matrix_t points;
    igraph_vector_t lengths;

    igraph_matrix_init(&points, 0, 0);
    igraph_vector_init(&lengths, 0);

    printf("K_4, 2D square\n");
    igraph_full(&graph, 4, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_layout_grid(&graph, &points, 2);
    printf("L2: ");
    igraph_spatial_edge_lengths(&graph, &lengths, &points, IGRAPH_METRIC_L2);
    print_vector(&lengths);
    printf("L1: ");
    igraph_spatial_edge_lengths(&graph, &lengths, &points, IGRAPH_METRIC_L1);
    print_vector(&lengths);

    igraph_matrix_resize(&points, 3, 2);
    CHECK_ERROR(igraph_spatial_edge_lengths(&graph, &lengths, &points, IGRAPH_METRIC_L2), IGRAPH_EINVAL);

    igraph_destroy(&graph);

    /* 0-by-0 matrix must be accepted for the null graph only. */
    printf("\nNull graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_matrix_resize(&points, 0, 0);
    printf("L2: ");
    igraph_spatial_edge_lengths(&graph, &lengths, &points, IGRAPH_METRIC_L2);
    print_vector(&lengths);
    printf("L1: ");
    igraph_spatial_edge_lengths(&graph, &lengths, &points, IGRAPH_METRIC_L1);
    print_vector(&lengths);
    igraph_destroy(&graph);

    igraph_vector_destroy(&lengths);
    igraph_matrix_destroy(&points);

    VERIFY_FINALLY_STACK();

    return 0;
}
