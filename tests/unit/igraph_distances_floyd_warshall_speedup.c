/*
   IGraph library.
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

int main(void) {
    igraph_t g;
    igraph_matrix_t d, d2;
    igraph_vector_t weights;

    igraph_matrix_init(&d, 0, 0);
    igraph_matrix_init(&d2, 0, 0);

    printf("Null graph\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);
    igraph_destroy(&g);

    printf("\nSingleton graph\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);
    igraph_destroy(&g);

    igraph_small(&g, 9, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,2, 2,3, 3,0, 0,4, 4,5, 3,0, 6,7, 5,4, 9,9,
                 -1);

    printf("\nUnweighted directed\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);

    igraph_distances(&g, &d2, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    printf("\nUnweighted directed, 'in' mode\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_IN, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);

    igraph_distances(&g, &d2, igraph_vss_all(), igraph_vss_all(), IGRAPH_IN);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    printf("\nUnweighted undirected\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_ALL, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);

    igraph_distances(&g, &d2, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    igraph_vector_init_int(&weights, igraph_ecount(&g),
                           2, 1, 5, 1, 2, 6, 8, 3, 3, 2, 3);

    printf("\nWeighted directed\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);

    igraph_distances_bellman_ford(&g, &d2, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    printf("\nWeighted undirected\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_ALL, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);

    igraph_distances_bellman_ford(&g, &d2, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_ALL);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    VECTOR(weights)[1] = -2;
    printf("\nNegative weight, directed\n");
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE);
    print_matrix(&d);

    igraph_distances_bellman_ford(&g, &d2, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    /* Check bad inputs */

    /* Negative weight edge in undirected graph */
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_ALL, IGRAPH_FLOYD_WARSHALL_TREE), IGRAPH_ENEGLOOP);

    /* Negative cycle in directed graph */
    VECTOR(weights)[1] = -10;
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), IGRAPH_ENEGLOOP);

    /* Negative self-loop in directed graph */
    VECTOR(weights)[1] = 1;
    VECTOR(weights)[10] = -1;
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), IGRAPH_ENEGLOOP);

    /* NaN weight */
    VECTOR(weights)[1] = IGRAPH_NAN;
    CHECK_ERROR(igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), IGRAPH_EINVAL);
    igraph_destroy(&g);

    /* Unweighted directed - larger graph */
    igraph_erdos_renyi_game_gnp(&g, 100, 0.1, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_distances_floyd_warshall(&g, &d, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE);
    igraph_distances_dijkstra(&g, &d2, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT);
    IGRAPH_ASSERT(igraph_matrix_all_e(&d, &d2));

    igraph_vector_destroy(&weights);

    igraph_destroy(&g);

    igraph_matrix_destroy(&d2);
    igraph_matrix_destroy(&d);

    VERIFY_FINALLY_STACK();
    return 0;
}
