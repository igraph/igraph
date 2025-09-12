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

#include "bench.h"
#include "igraph_spatial.h"

#include <float.h>

void bench_lune(igraph_int_t test_nr, igraph_int_t point_count, igraph_int_t dimensions, igraph_real_t beta) {
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
             "%3"IGRAPH_PRId" Lune based beta skeleton in %dD, beta=%.3f, n=%6" IGRAPH_PRId,
             test_nr, (int) dimensions, beta, point_count);

    BENCH(msg, igraph_lune_beta_skeleton(&g, &points, beta));

    igraph_destroy(&g);
    igraph_matrix_destroy(&points);

}

void bench_circle(igraph_int_t test_nr, igraph_int_t point_count, igraph_int_t dimensions, igraph_real_t beta) {
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
             "%3"IGRAPH_PRId" Circle based beta skeleton in %dD, beta=%3f, n=%6" IGRAPH_PRId,
             test_nr, (int) dimensions, beta, point_count);

    BENCH(msg, igraph_circle_beta_skeleton(&g, &points, beta));

    igraph_destroy(&g);
    igraph_matrix_destroy(&points);

}

void bench_gabriel(igraph_int_t test_nr, igraph_int_t point_count, igraph_int_t dimensions, igraph_real_t beta_cutoff) {
    igraph_matrix_t points;
    igraph_vector_t edge_weights;
    igraph_t g;
    char msg[200];

    igraph_matrix_init(&points, point_count, dimensions);
    igraph_vector_init(&edge_weights, 0);

    for (igraph_int_t point = 0; point < point_count; point++) {
        for (igraph_int_t dim = 0; dim < dimensions; dim++) {
            MATRIX(points, point, dim) = RNG_UNIF01();
        }
    }


    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%3"IGRAPH_PRId" Beta weighted gabriel graph in %dD, cutoff=%3f, n=%6" IGRAPH_PRId,
             test_nr, (int) dimensions, beta_cutoff, point_count);

    BENCH(msg, igraph_beta_weighted_gabriel_graph(&g, &edge_weights,  &points, beta_cutoff));

    igraph_vector_destroy(&edge_weights);
    igraph_destroy(&g);
    igraph_matrix_destroy(&points);

}

int main(void) {
    igraph_int_t test_nr = 0;
    bench_lune(++test_nr, 100, 2, 2);
    bench_lune(++test_nr, 1000, 2, 2);
    bench_lune(++test_nr, 10000, 2, 2);
    bench_lune(++test_nr, 100000, 2, 2);
    bench_lune(++test_nr, 100, 3, 2);
    bench_lune(++test_nr, 1000, 3, 2);
    bench_lune(++test_nr, 10000, 3, 2);
    bench_lune(++test_nr, 100000, 3, 2);
    bench_lune(++test_nr, 1000, 4, 2);

    bench_lune(++test_nr, 10000, 2, 1);
    bench_lune(++test_nr, 10000, 2, 2);
    bench_lune(++test_nr, 10000, 2, 4);
    bench_lune(++test_nr, 10000, 2, 8);
    bench_lune(++test_nr, 10000, 2, 16);
    bench_lune(++test_nr, 10000, 2, 32);
    bench_lune(++test_nr, 10000, 2, 64);
    bench_lune(++test_nr, 10000, 2, 128);
    bench_lune(++test_nr, 10000, 2, 256);
    bench_lune(++test_nr, 10000, 2, 512);
    bench_lune(++test_nr, 10000, 2, 1024);
    bench_lune(++test_nr, 1000, 2, 0.9);
    bench_lune(++test_nr, 1000, 2, 0.5);
    bench_lune(++test_nr, 1000, 2, 0.1);


    bench_circle(++test_nr, 100, 2, 2);
    bench_circle(++test_nr, 1000, 2, 2);
    bench_circle(++test_nr, 10000, 2, 2);
    bench_circle(++test_nr, 100000, 2, 2);

    bench_circle(++test_nr, 10000, 2, 2);
    bench_circle(++test_nr, 10000, 2, 3);
    bench_circle(++test_nr, 10000, 2, 5);
    bench_circle(++test_nr, 10000, 2, 15);
    bench_circle(++test_nr, 10000, 2, 50);

    bench_gabriel(++test_nr, 100, 2, 2);
    bench_gabriel(++test_nr, 1000, 2, 2);
    bench_gabriel(++test_nr, 10000, 2, 2);
    bench_gabriel(++test_nr, 100000, 2, 2);

    bench_gabriel(++test_nr, 100, 3, 2);
    bench_gabriel(++test_nr, 1000, 3, 2);
    bench_gabriel(++test_nr, 10000, 3, 2);
    bench_gabriel(++test_nr, 100000, 3, 2);
    return 0;
}
