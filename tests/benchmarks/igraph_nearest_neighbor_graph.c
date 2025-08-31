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

#include "bench.h"

void bench_pointcloud(const char *message, igraph_metric_t metric, igraph_int_t point_count, igraph_int_t dimensions, igraph_int_t neighbors, igraph_real_t cutoff) {
    igraph_matrix_t points;
    igraph_t g;
    char msg[200];
    const char *metricname;

    igraph_matrix_init(&points, point_count, dimensions);

    for (igraph_int_t point = 0; point < point_count; point++) {
        for (igraph_int_t dim = 0; dim < dimensions; dim++) {
            MATRIX(points, point, dim) = RNG_UNIF01();
        }
    }

    switch (metric) {
    case IGRAPH_METRIC_L2: metricname = "L2"; break;
    case IGRAPH_METRIC_L1: metricname = "L1"; break;
    }

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%s RKNN with %s, %dD, n=%6" IGRAPH_PRId ", k=%2d, r=%g",
             message, metricname, (int) dimensions, point_count, (int) neighbors, cutoff);

    BENCH(msg, igraph_nearest_neighbor_graph(&g, &points, metric, neighbors, cutoff, false));

    igraph_destroy(&g);
    igraph_matrix_destroy(&points);
}

int main(void) {
    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 20250728);

    bench_pointcloud(" 1", IGRAPH_METRIC_L2, 100000, 1, 10, 0.1);
    bench_pointcloud(" 2", IGRAPH_METRIC_L2, 100000, 2, 10, 0.1);
    bench_pointcloud(" 3", IGRAPH_METRIC_L2, 100000, 3, 10, 0.1);
    bench_pointcloud(" 4", IGRAPH_METRIC_L2, 100000, 4, 10, 0.1);

    bench_pointcloud(" 5", IGRAPH_METRIC_L2, 10000, 1, -1, 0.001);
    bench_pointcloud(" 6", IGRAPH_METRIC_L2, 10000, 2, -1, 0.01);
    bench_pointcloud(" 7", IGRAPH_METRIC_L2, 10000, 3, -1, 0.01);
    bench_pointcloud(" 8", IGRAPH_METRIC_L2, 10000, 4, -1, 0.1);

    bench_pointcloud(" 9",  IGRAPH_METRIC_L2, 1000, 1, -1, -1);
    bench_pointcloud("10",  IGRAPH_METRIC_L2, 1000, 2, -1, -1);
    bench_pointcloud("11",  IGRAPH_METRIC_L2, 1000, 3, -1, -1);
    bench_pointcloud("12",  IGRAPH_METRIC_L2, 1000, 4, -1, -1);

    bench_pointcloud("13",  IGRAPH_METRIC_L2, 100000, 1, 10, -1);
    bench_pointcloud("14",  IGRAPH_METRIC_L2, 100000, 2, 10, -1);
    bench_pointcloud("15",  IGRAPH_METRIC_L2, 100000, 3, 10, -1);
    bench_pointcloud("16",  IGRAPH_METRIC_L2, 100000, 4, 10, -1);

    return 0;
}
