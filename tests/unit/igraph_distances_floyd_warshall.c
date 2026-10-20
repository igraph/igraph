/*
   igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

void test_small_graph_selectors(void) {
    const igraph_floyd_warshall_algorithm_t methods[] = {
        IGRAPH_FLOYD_WARSHALL_ORIGINAL,
        IGRAPH_FLOYD_WARSHALL_TREE,
        IGRAPH_FLOYD_WARSHALL_AUTOMATIC
    };
    igraph_matrix_t result;
    igraph_vector_t weights;
    igraph_vector_int_t sources;

    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init(&weights, 0);
    igraph_vector_int_init(&sources, 0);

    for (igraph_int_t n = 0; n <= 1; n++) {
        /* Repeated sources are permitted and must produce repeated rows. */
        igraph_vector_int_resize(&sources, 2 * n);
        igraph_vector_int_null(&sources);
        const igraph_vs_t from[] = {igraph_vss_all(), igraph_vss_none(), igraph_vss_vector(&sources)};
        const igraph_vs_t to[] = {igraph_vss_all(), igraph_vss_none()};
        const igraph_int_t nrows[] = {n, 0, 2 * n};
        const igraph_int_t ncols[] = {n, 0};

        for (int directed = 0; directed < 2; directed++) {
            igraph_t graph;
            igraph_empty(&graph, n, directed);
            for (int weighted = 0; weighted < 2; weighted++) {
                for (int method = 0; method < 3; method++) {
                    for (int f = 0; f < 3; f++) {
                        for (int t = 0; t < 2; t++) {
                            IGRAPH_ASSERT(igraph_distances_floyd_warshall(
                                &graph, &result, from[f], to[t], weighted ? &weights : NULL,
                                IGRAPH_OUT, methods[method]) == IGRAPH_SUCCESS);
                            IGRAPH_ASSERT(igraph_matrix_nrow(&result) == nrows[f]);
                            IGRAPH_ASSERT(igraph_matrix_ncol(&result) == ncols[t]);
                            for (igraph_int_t i = 0; i < nrows[f]; i++) {
                                for (igraph_int_t j = 0; j < ncols[t]; j++) {
                                    IGRAPH_ASSERT(MATRIX(result, i, j) == 0);
                                }
                            }
                        }
                    }
                }
            }
            igraph_destroy(&graph);
        }
    }

    igraph_vector_int_destroy(&sources);
    igraph_vector_destroy(&weights);
    igraph_matrix_destroy(&result);
}

int main(void) {
    igraph_t g;
    igraph_matrix_t d, d2;
    igraph_vector_t weights;

    igraph_matrix_init(&d, 0, 0);
    igraph_matrix_init(&d2, 0, 0);

    printf("Null graph\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);
    igraph_destroy(&g);

    printf("\nSingleton graph\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);
    igraph_destroy(&g);

    igraph_small(&g, 9, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,2, 2,3, 3,0, 0,4, 4,5, 3,0, 6,7, 5,4, 9,9,
                 -1);

    printf("\nUnweighted directed\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);

    igraph_distances(&g, NULL, &d2, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    printf("\nUnweighted directed, 'in' mode\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_IN, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);

    igraph_distances(&g, NULL, &d2, igraph_vss_all(), igraph_vss_all(), IGRAPH_IN);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    printf("\nUnweighted undirected\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_ALL, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);

    igraph_distances(&g, NULL, &d2, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    igraph_vector_init_int(&weights, igraph_ecount(&g),
                           2, 1, 5, 1, 2, 6, 8, 3, 3, 2, 3);

    printf("\nWeighted directed\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);

    igraph_distances_bellman_ford(&g, &d2, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    printf("\nWeighted undirected\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_ALL, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);

    igraph_distances_bellman_ford(&g, &d2, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_ALL);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    VECTOR(weights)[1] = -2;
    printf("\nNegative weight, directed\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL);
    print_matrix(&d);

    igraph_distances_bellman_ford(&g, &d2, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    /* Check bad inputs */

    /* Negative weight edge in undirected graph */
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_ALL, IGRAPH_FLOYD_WARSHALL_ORIGINAL), IGRAPH_ENEGCYCLE);

    /* Negative cycle in directed graph */
    VECTOR(weights)[1] = -10;
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), IGRAPH_ENEGCYCLE);

    /* Negative self-loop in directed graph */
    VECTOR(weights)[1] = 1;
    VECTOR(weights)[10] = -1;
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), IGRAPH_ENEGCYCLE);

    /* NaN weight */
    VECTOR(weights)[1] = IGRAPH_NAN;
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), IGRAPH_EINVAL);

    igraph_vector_destroy(&weights);

    igraph_destroy(&g);

    igraph_matrix_destroy(&d2);
    igraph_matrix_destroy(&d);

    test_small_graph_selectors();

    VERIFY_FINALLY_STACK();
    return 0;
}
