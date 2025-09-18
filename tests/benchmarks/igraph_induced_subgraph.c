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

// Function to generate random vertex indices for subgraph extraction
void generate_random_vertices(igraph_vector_int_t *vertices,
                              igraph_int_t graph_size,
                              igraph_int_t subset_size) {
    igraph_vector_int_init(vertices, 0);
    igraph_vector_int_resize(vertices, subset_size);

    // Create a vector to track used vertices
    igraph_vector_int_t used_vertices;
    igraph_vector_int_init(&used_vertices, graph_size);
    for (igraph_int_t i = 0; i < graph_size; i++) {
        VECTOR(used_vertices)[i] = i;
    }

    // Randomly select subset_size unique vertices
    for (igraph_int_t i = 0; i < subset_size; i++) {
        igraph_int_t index = RNG_INTEGER(0, graph_size - i - 1);
        VECTOR(*vertices)[i] = VECTOR(used_vertices)[index];

        // Remove the selected vertex from used_vertices
        VECTOR(used_vertices)[index] = VECTOR(used_vertices)[graph_size - i - 1];
    }

    // Sort the vertices to simulate real-world scenarios
    igraph_vector_int_sort(vertices);

    igraph_vector_int_destroy(&used_vertices);
}

// Benchmark function
void bench_induced_subgraph(igraph_t *graph, igraph_vector_int_t *vertices,
                            igraph_subgraph_implementation_t impl) {
    igraph_t subgraph;
    igraph_vs_t vs;

    igraph_vs_vector(&vs, vertices);
    igraph_induced_subgraph(graph, &subgraph, vs, impl);
    igraph_vs_destroy(&vs);
    igraph_destroy(&subgraph);
}

void run_bench(int i, int n, int m, int subset_percentage) {
    igraph_t graph;
    igraph_vector_int_t vertices;

    // Calculate subset size based on percentage
    igraph_int_t subset_size = (n * subset_percentage) / 100;

    // Prepare random graph and vertices
    igraph_erdos_renyi_game_gnm(&graph, n, m, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    generate_random_vertices(&vertices, n, subset_size);

    char msg[255];
    int rep = 300000000 / (n * subset_percentage);

    // Benchmark method 1 (COPY_AND_DELETE)
    snprintf(msg, sizeof(msg),
             "Method 1 (COPY_AND_DELETE):     n=%5d, subset=%3d%%, %dx", n,
             subset_percentage, rep);
    BENCH(msg, REPEAT(bench_induced_subgraph(&graph, &vertices,
                                             IGRAPH_SUBGRAPH_COPY_AND_DELETE),
                      rep));

    // Benchmark method 2 (CREATE_FROM_SCRATCH)
    snprintf(msg, sizeof(msg),
             "Method 2 (CREATE_FROM_SCRATCH): n=%5d, subset=%3d%%, %dx", n,
             subset_percentage, rep);
    BENCH(msg, REPEAT(bench_induced_subgraph(&graph, &vertices,
                                             IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH),
                      rep));

    // Cleanup
    igraph_vector_int_destroy(&vertices);
    igraph_destroy(&graph);
}

int main(void) {
    int i = 0;

    // Initialize random number generator
    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

// Macro to run benchmarks for different graph sizes and subset percentages
#define BENCHSET(n, m)                                                           \
    run_bench(++i, n, m, 20); /* 20% subset */                                   \
    run_bench(++i, n, m, 25); /* 25% subset */                                   \
    run_bench(++i, n, m, 30); /* 30% subset */                                   \
    run_bench(++i, n, m, 35); /* 35% subset */                                   \
    run_bench(++i, n, m, 40); /* 40% subset */                                   \
    run_bench(++i, n, m, 45); /* 45% subset */                                   \
    run_bench(++i, n, m, 50); /* 50% subset */                                   \
    run_bench(++i, n, m, 55); /* 55% subset */                                   \
    printf("\n");

    // Benchmark different graph sizes
    BENCHSET(100, 500);
    BENCHSET(1000, 5000);
    BENCHSET(10000, 50000);
    BENCHSET(100000, 500000);

    return 0;
}
