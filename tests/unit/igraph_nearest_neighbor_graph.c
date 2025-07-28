/*
   IGraph library.
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

#include <math.h>

double sqr(double x) { return x*x; }

igraph_error_t RKNN_neighbors(
        igraph_real_t *arr,
        igraph_integer_t num_points,
        igraph_integer_t dims,
        igraph_integer_t neighbors,
        igraph_real_t cutoff) {

    igraph_t graph;
    igraph_matrix_t points;
    igraph_matrix_t adj_mat;
    igraph_metric_t metric = IGRAPH_METRIC_L2;
    igraph_matrix_t dist_mat;

    IGRAPH_MATRIX_INIT_FINALLY(&dist_mat, num_points, num_points);

    IGRAPH_FINALLY(igraph_matrix_destroy, &points);
    igraph_matrix_init_array(&points, &arr[0], num_points, dims, IGRAPH_ROW_MAJOR);

    for (igraph_integer_t start = 0; start < num_points-1; start++ ) {
        for (igraph_integer_t end = start + 1; end < num_points; end++) {
            igraph_real_t distance = 0;
            switch (metric) {
                case IGRAPH_METRIC_L2:
                    for (igraph_integer_t i = 0; i < dims; i++) {
                        distance += sqr(MATRIX(points, start, i) - MATRIX(points, end, i));
                    }
                    break;
                case IGRAPH_METRIC_L1:
                    for (igraph_integer_t i = 0; i < dims; i++) {
                        distance += fabs(MATRIX(points, start, i) - MATRIX(points, end, i));
                    }
                    break;
            }
            MATRIX(dist_mat, start,end) = distance;
            MATRIX(dist_mat, end,start) = distance;
        }
    }


    IGRAPH_CHECK(igraph_nearest_neighbor_graph(&graph, &points, metric, neighbors, cutoff, IGRAPH_DIRECTED));
    print_matrix(&points);
    print_graph_canon(&graph);
    igraph_matrix_init(&adj_mat, 0, 0);
    IGRAPH_FINALLY(igraph_matrix_destroy, &adj_mat);
    IGRAPH_CHECK(igraph_get_adjacency(&graph, &adj_mat, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_destroy, &graph);

    for (igraph_integer_t start = 0; start < num_points; start ++) {

            igraph_real_t min_missing=INFINITY;
            igraph_real_t max_present=0;
        for (igraph_integer_t end = 0; end < num_points; end ++) {
            if (start == end) continue;
            if (MATRIX(adj_mat, start, end) == 1) {
                if (max_present < MATRIX(dist_mat, start, end)) max_present = MATRIX(dist_mat, start, end);
            } else {
                if (min_missing > MATRIX(dist_mat, start, end)) min_missing = MATRIX(dist_mat, start, end);
            }
        }
        IGRAPH_ASSERT(min_missing >= max_present);
        IGRAPH_ASSERT(max_present <= cutoff*cutoff);
    }

    igraph_destroy(&graph);
    igraph_matrix_destroy(&adj_mat);
    igraph_matrix_destroy(&points);
    igraph_matrix_destroy(&dist_mat);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_real_t points1d[5] = {
        12,
        8,
        5,
        10,
        12
    };
    igraph_real_t points2d[10] = {
        12  , 8,
        8   , 6,
        5   , 12,
        10  , 1,
        12  , 2
    };
    igraph_real_t points3d[15] = {
        1,6,4,
        6,2,3,
        3,6,6,
        3,2,2,
        2,3,3
    };

    igraph_real_t points4d[20] = {
        1,6,4,4,
        6,2,3,3,
        3,6,6,5,
        3,2,2,3,
        2,3,3,4
    };

    printf("1d 2 neighbors, cutoff 3\n");
    RKNN_neighbors(&points1d[0], 5, 1, 2, 3);
    printf("1d 1 neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points1d[0], 5, 1, 1, INFINITY);
    printf("1d unlimited neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points1d[0], 5, 1, -1, INFINITY);
    printf("1d unlimited neighbors, cutoff 7\n");
    RKNN_neighbors(&points1d[0], 5, 1, -1, 3);

    printf("2d 2 neighbors, cutoff 5\n");
    RKNN_neighbors(&points2d[0], 5, 2, 2, 5);
    printf("2d 1 neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points2d[0], 5, 2, 1, INFINITY);
    printf("2d unlimited neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points2d[0], 5, 2, -1, INFINITY);
    printf("2d unlimited neighbors, cutoff 7\n");
    RKNN_neighbors(&points2d[0], 5, 2, -1, 7);

    printf("3d, 2 neighbors, cutoff 4\n");
    RKNN_neighbors(&points3d[0], 5, 3, 2, 4);
    printf("3d, unlimited neighbors, cutoff 4\n");
    RKNN_neighbors(&points3d[0], 5, 3, -1, 4);
    printf("3d, 2 neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points3d[0], 5, 3, 2, INFINITY);
    printf("3d, unlimited neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points3d[0], 5, 3, -1, INFINITY);

    printf("4d 2 neighbors, cutoff 5\n");
    RKNN_neighbors(&points4d[0], 5, 4, 2, 5);
    printf("4d 1 neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points4d[0], 5, 4, 1, INFINITY);
    printf("4d unlimited neighbors, cutoff INFINITY\n");
    RKNN_neighbors(&points4d[0], 5, 4, -1, INFINITY);
    printf("4d unlimited neighbors, cutoff 4\n");
    RKNN_neighbors(&points4d[0], 5, 4, -1, 4);

    printf("3d, no points\n");
    RKNN_neighbors(&points3d[0], 0, 3, -1, 100);
    printf("3d, 1 point\n");
    RKNN_neighbors(&points3d[0], 1, 3, -1, 100);

    printf("0d, should error\n");
    CHECK_ERROR(RKNN_neighbors(&points3d[0], 0, 0, 10, INFINITY), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();
    return 0;
}
