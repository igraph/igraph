/*
   igraph library.
   Copyright (C) 2013-2026  The igraph development team <igraph@igraph.org>

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
#include "core/barnes_hut.h"


/* * Dummy force calculation for benchmarking.
 * We use a simple inverse distance approximation to simulate math overhead.
 * The geometry (dx, dy, dz, dist_sq) is now pre-calculated by the BH Engine.
 */
static void bench_repulsion_kernel(
    const igraph_bh_point_t *p1,
    const igraph_bh_point_t *p2,
    igraph_real_t dx,
    igraph_real_t dy,
    igraph_real_t dz,
    igraph_real_t dist_sq,
    igraph_real_t force[3],
    void *user_data
) {
    IGRAPH_UNUSED(p1);
    IGRAPH_UNUSED(p2);
    IGRAPH_UNUSED(user_data);

    /* Add a tiny epsilon to simulate the old benchmark math overhead */
    dist_sq += 1e-6;
    igraph_real_t f = 1.0 / dist_sq;

    force[0] = f * dx;
    force[1] = f * dy;
    force[2] = f * dz;
}

/* * Used to prevent aggressive compiler optimizations from removing
 * the force calculations during the benchmark loop.
 */
volatile igraph_real_t dummy_sum = 0;

void do_benchmarks(const char *name, igraph_integer_t n_points, igraph_integer_t dim, int repeat) {
    igraph_bh_tree_t tree;
    igraph_matrix_t coords;
    igraph_matrix_t forces;

    /* 1. Generate Random Coordinates */
    igraph_matrix_init(&coords, n_points, dim);
    igraph_matrix_init(&forces, n_points, dim);

    for (igraph_integer_t i = 0; i < n_points; i++) {
        for (igraph_integer_t d = 0; d < dim; d++) {
            /* Random values between -100.0 and 100.0 */
            MATRIX(coords, i, d) = igraph_rng_get_unif(igraph_rng_default(), -100.0, 100.0);
        }
    }

    printf("\n--- %s (%d points, %dD) ---\n", name, (int)n_points, (int)dim);

    /* --- Benchmark 1: Memory Management Overhead --- */
    BENCH("1 init/destroy tree (allocation baseline)",
          REPEAT(
              do {
                  igraph_bh_tree_init(&tree, dim, 0.6, 15, 1);
                  igraph_bh_tree_destroy(&tree);
              } while (0),
          repeat);
    );

    /* --- Benchmark 2: Tree Construction (Default Capacity) --- */
    BENCH("2 tree build (leaf_cap = 1)",
          REPEAT(
              do {
                  igraph_bh_tree_init(&tree, dim, 0.6, 15, 1);
                  igraph_bh_tree_build(&tree, &coords, NULL);
                  igraph_bh_tree_destroy(&tree);
              } while (0),
          repeat);
    );

    /* --- Benchmark 3: Tree Construction (High Capacity) ---
     * Higher capacity means fewer nodes and shallower trees, which
     * drastically speeds up array partitioning and memory use.
     */
    BENCH("3 tree build (leaf_cap = 8)",
          REPEAT(
              do {
                  igraph_bh_tree_init(&tree, dim, 0.6, 15, 8);
                  igraph_bh_tree_build(&tree, &coords, NULL);
                  igraph_bh_tree_destroy(&tree);
              } while (0),
          repeat);
    );

    /* --- Benchmark 4: Repulsive Force Calculation (Theta = 0.6) --- */
    igraph_bh_tree_init(&tree, dim, 0.6, 15, 1);
    igraph_bh_tree_build(&tree, &coords, NULL);

    BENCH("4 repulsive forces (theta = 0.6, cap = 1)",
          REPEAT(
              do {
                  igraph_bh_apply_repulsion_from_tree(&tree, &forces, bench_repulsion_kernel, NULL);
                  dummy_sum += MATRIX(forces, 0, 0); /* Force evaluation */
              } while (0),
          repeat);
    );
    igraph_bh_tree_destroy(&tree);

    /* --- Benchmark 5: Repulsive Force Calculation (Theta = 0.8) ---
     * A higher theta approximates more aggressively, skipping deeper tree descents.
     * This should run significantly faster than theta = 0.6.
     */
    igraph_bh_tree_init(&tree, dim, 0.8, 15, 1);
    igraph_bh_tree_build(&tree, &coords, NULL);

    BENCH("5 repulsive forces (theta = 0.8, cap = 1)",
          REPEAT(
              do {
                  igraph_bh_apply_repulsion_from_tree(&tree, &forces, bench_repulsion_kernel, NULL);
                  dummy_sum += MATRIX(forces, 0, 0);
              } while (0),
          repeat);
    );
    igraph_bh_tree_destroy(&tree);

    /* --- Benchmark 6: Vectorized Leaf Calculation Potential --- */
    igraph_bh_tree_init(&tree, dim, 0.6, 15, 8);
    igraph_bh_tree_build(&tree, &coords, NULL);

    BENCH("6 repulsive forces (theta = 0.6, cap = 8)",
          REPEAT(
              do {
                  igraph_bh_apply_repulsion_from_tree(&tree, &forces, bench_repulsion_kernel, NULL);
                  dummy_sum += MATRIX(forces, 0, 0);
              } while (0),
          repeat);
    );
    igraph_bh_tree_destroy(&tree);

    /* Cleanup */
    igraph_matrix_destroy(&coords);
    igraph_matrix_destroy(&forces);
}


int main(void) {
    /* Ensure reproducible random states for coordinates */
    igraph_rng_seed(igraph_rng_default(), 42);

    BENCH_INIT();

    /* * Test Scale 1: Small Graph (e.g., standard hairball layout)
     * High repeat count to get accurate millisecond timings
     */
    do_benchmarks("Small Scale 2D", 1000, 2, 100);
    do_benchmarks("Small Scale 3D", 1000, 3, 100);

    /* * Test Scale 2: Medium Graph
     * Moderate repeat count
     */
    do_benchmarks("Medium Scale 2D", 10000, 2, 10);
    do_benchmarks("Medium Scale 3D", 10000, 3, 10);

    /* * Test Scale 3: Large Graph (e.g., massive social network layout)
     * Low repeat count due to O(N log N) scaling
     */
    do_benchmarks("Large Scale 2D", 50000, 2, 2);
    do_benchmarks("Large Scale 3D", 50000, 3, 2);

    return 0;
}
