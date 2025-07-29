/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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
#include "bench.h"
#include "igraph_interface.h"
#include "igraph_random.h"
#include "igraph_spatial.h"


void bench_pointcloud(const char *message, igraph_metric_t metric, igraph_integer_t point_count, igraph_integer_t dimensions, igraph_integer_t neighbors, igraph_real_t cutoff) {
    igraph_matrix_t points;
    igraph_t g;
    igraph_matrix_init(&points, point_count, dimensions);

    for (igraph_integer_t point = 0; point < point_count; point++) {
        for (igraph_integer_t dim = 0; dim < dimensions; dim++) {
            MATRIX(points, point, dim) = igraph_rng_get_unif(igraph_rng_default(), 0, 1);
        }
    }
    BENCH(message, igraph_nearest_neighbor_graph(&g, &points, metric, neighbors, cutoff, false));
    igraph_destroy(&g);
    igraph_matrix_destroy(&points);
}

int main(void) {
    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 2025-07-28);

    bench_pointcloud(" 1 1D RKNN construction with 100000 points in a cube k=10, r=0.1", IGRAPH_METRIC_L2, 100000, 1, 10, 0.1);
    bench_pointcloud(" 2 2D RKNN construction with 100000 points in a cube k=10, r=0.1", IGRAPH_METRIC_L2, 100000, 2, 10, 0.1);
    bench_pointcloud(" 3 3D RKNN construction with 100000 points in a cube k=10, r=0.1", IGRAPH_METRIC_L2, 100000, 3, 10, 0.1);
    bench_pointcloud(" 4 4D RKNN construction with 100000 points in a cube k=10, r=0.1", IGRAPH_METRIC_L2, 100000, 4, 10, 0.1);

    bench_pointcloud(" 5 1D RKNN construction with 100000 points in a cube k=-1, r=0.1", IGRAPH_METRIC_L2, 10000, 1, -1, 0.001);
    bench_pointcloud(" 6 2D RKNN construction with 100000 points in a cube k=-1, r=0.1", IGRAPH_METRIC_L2, 10000, 2, -1, 0.01);
    bench_pointcloud(" 7 3D RKNN construction with 100000 points in a cube k=-1, r=0.1", IGRAPH_METRIC_L2, 10000, 3, -1, 0.01);
    bench_pointcloud(" 8 4D RKNN construction with 100000 points in a cube k=-1, r=0.1", IGRAPH_METRIC_L2, 10000, 4, -1, 0.1);

    bench_pointcloud(" 9 1D RKNN construction with 100000 points in a cube k=-1, r=-1",  IGRAPH_METRIC_L2, 1000, 1, -1, -1);
    bench_pointcloud("10 2D RKNN construction with 100000 points in a cube k=-1, r=-1",  IGRAPH_METRIC_L2, 1000, 2, -1, -1);
    bench_pointcloud("11 3D RKNN construction with 100000 points in a cube k=-1, r=-1",  IGRAPH_METRIC_L2, 1000, 3, -1, -1);
    bench_pointcloud("12 4D RKNN construction with 100000 points in a cube k=-1, r=-1",  IGRAPH_METRIC_L2, 1000, 4, -1, -1);

    bench_pointcloud("13 1D RKNN construction with 100000 points in a cube k=10, r=-1",  IGRAPH_METRIC_L2, 100000, 1, 10, -1);
    bench_pointcloud("14 2D RKNN construction with 100000 points in a cube k=10, r=-1",  IGRAPH_METRIC_L2, 100000, 2, 10, -1);
    bench_pointcloud("15 3D RKNN construction with 100000 points in a cube k=10, r=-1",  IGRAPH_METRIC_L2, 100000, 3, 10, -1);
    bench_pointcloud("16 4D RKNN construction with 100000 points in a cube k=10, r=-1",  IGRAPH_METRIC_L2, 100000, 4, 10, -1);

    
}
