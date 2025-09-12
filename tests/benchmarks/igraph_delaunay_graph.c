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

void bench_delaunay(const char *message, igraph_int_t point_count, igraph_int_t dimensions) {
    igraph_matrix_t points;
    igraph_t g;
    char msg[200];

    igraph_matrix_init(&points, point_count, dimensions);

    for (igraph_int_t point = 0; point < point_count; point++) {
        for (igraph_int_t dim = 0; dim < dimensions; dim++) {
            MATRIX(points, point, dim) = RNG_UNIF01();
        }
    }

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%s Delaunay triangulation in %dD, n=%6" IGRAPH_PRId,
             message, (int) dimensions, point_count);

    BENCH(msg, igraph_delaunay_graph(&g, &points));

    igraph_destroy(&g);
    igraph_matrix_destroy(&points);
}

int main(void) {
    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 20250814);

    bench_delaunay(" 1", 1000000, 1);
    bench_delaunay(" 2", 100000, 2);
    bench_delaunay(" 3", 10000, 2);
    bench_delaunay(" 4", 1000, 2);
    bench_delaunay(" 5", 100, 2);
    bench_delaunay(" 6", 10000, 3);
    bench_delaunay(" 7", 1000, 3);
    bench_delaunay(" 8", 100, 3);
    bench_delaunay(" 9", 2000, 4);
    bench_delaunay("10", 500, 5);
    bench_delaunay("11", 200, 6);

    return 0;
}
