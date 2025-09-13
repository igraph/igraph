/*
   igraph library.
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

/*
 * Benchmark creating graphs from dense adjacency matrices.
 *
 * When there are a small number of non-zero elements (low mean degree),
 * iterating through the matrix dominates the timing. When there are
 * many non-zero elements, creating the graphfrom its edge list dominates.
 */

void adjacency(const igraph_matrix_t *adjmatrix, igraph_adjacency_t mode, igraph_loops_t loops) {
    igraph_t g;
    igraph_adjacency(&g, adjmatrix, mode, loops);
    igraph_destroy(&g);
}

void weighted_adjacency(const igraph_matrix_t *adjmatrix, igraph_adjacency_t mode, igraph_vector_t *weights, igraph_loops_t loops) {
    igraph_t g;
    igraph_weighted_adjacency(&g, adjmatrix, mode, weights, loops);
    igraph_destroy(&g);
}

void run_bench(igraph_int_t vcount, igraph_int_t meandeg, igraph_int_t rep) {
    igraph_t g;
    igraph_matrix_t mat;
    igraph_vector_t weights;
    char msg[128];
    igraph_adjacency_t types[] = {
        IGRAPH_ADJ_DIRECTED,
        IGRAPH_ADJ_MAX,
        IGRAPH_ADJ_PLUS,
        IGRAPH_ADJ_UPPER /* similar to DIRECTED when unweighted, similar to MAX when weighted */
    };
    const char *names[] = { "DIRECTED", "MAX", "PLUS", "UPPER" };

    igraph_matrix_init(&mat, 0, 0);

    igraph_erdos_renyi_game_gnm(&g, vcount, meandeg * vcount / 2, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW,
                                false);
    igraph_get_adjacency(&g, &mat, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_ONCE);

    igraph_vector_init(&weights, igraph_ecount(&g));

    igraph_destroy(&g);

    for (size_t i=0; i < sizeof(types) / sizeof(types[0]); i++) {
        snprintf(msg, sizeof(msg) / sizeof(msg[0]),
                 "%2d vcount=%" IGRAPH_PRId ", meandeg=%3" IGRAPH_PRId ", %8s, unweighted, %" IGRAPH_PRId "x",
                 (int) i+1, vcount, meandeg, names[i], rep);

        BENCH(msg, REPEAT(adjacency(&mat, types[i], IGRAPH_LOOPS_ONCE), rep));

        snprintf(msg, sizeof(msg) / sizeof(msg[0]),
                 "%2d vcount=%" IGRAPH_PRId ", meandeg=%3" IGRAPH_PRId ", %8s,   weighted, %" IGRAPH_PRId "x",
                 (int) i+1, vcount, meandeg, names[i], rep);

        BENCH(msg, REPEAT(weighted_adjacency(&mat, types[i], &weights, IGRAPH_LOOPS_ONCE), rep));
    }
    printf("\n");

    igraph_vector_destroy(&weights);
    igraph_matrix_destroy(&mat);
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    run_bench(100, 5, 10000);
    run_bench(100, 50, 10000);
    run_bench(1000, 5, 100);
    run_bench(1000, 50, 100);
    run_bench(1000, 500, 100);
    run_bench(10000, 5, 1);
    run_bench(10000, 50, 1);
    run_bench(10000, 500, 1);

    return 0;
}
