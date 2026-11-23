/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

/*
 * Cross-algorithm consistency checks for shortest-distance functions.
 *
 * For graphs with non-negative weights, igraph_distances_dijkstra(),
 * igraph_distances_bellman_ford(), and igraph_distances_floyd_warshall()
 * must all return identical distance matrices.
 * For the unweighted case (NULL weights), igraph_distances() (BFS) must
 * also agree.
 */

#include <igraph.h>
#include "test_utilities.h"

/* Run Dijkstra, Bellman-Ford, and Floyd-Warshall and assert all three
 * produce the same distance matrix.  When weights is NULL (unweighted),
 * also compare against BFS via igraph_distances(). */
static void check_all_agree(const igraph_t *g,
                             const igraph_vector_t *weights,
                             igraph_neimode_t mode) {
    igraph_matrix_t d_dijk, d_bf, d_fw, d_bfs;

    igraph_matrix_init(&d_dijk, 0, 0);
    igraph_matrix_init(&d_bf,   0, 0);
    igraph_matrix_init(&d_fw,   0, 0);
    igraph_matrix_init(&d_bfs,  0, 0);

    igraph_distances_dijkstra(g, &d_dijk,
        igraph_vss_all(), igraph_vss_all(), weights, mode);
    igraph_distances_bellman_ford(g, &d_bf,
        igraph_vss_all(), igraph_vss_all(), weights, mode);
    igraph_distances_floyd_warshall(g, &d_fw,
        igraph_vss_all(), igraph_vss_all(), weights, mode,
        IGRAPH_FLOYD_WARSHALL_AUTOMATIC);

    IGRAPH_ASSERT(igraph_matrix_all_e(&d_dijk, &d_bf));
    IGRAPH_ASSERT(igraph_matrix_all_e(&d_dijk, &d_fw));

    if (weights == NULL) {
        igraph_distances(g, NULL, &d_bfs,
            igraph_vss_all(), igraph_vss_all(), mode);
        IGRAPH_ASSERT(igraph_matrix_all_e(&d_dijk, &d_bfs));
    }

    igraph_matrix_destroy(&d_bfs);
    igraph_matrix_destroy(&d_fw);
    igraph_matrix_destroy(&d_bf);
    igraph_matrix_destroy(&d_dijk);
}

int main(void) {
    igraph_t g;
    igraph_vector_t weights;
    igraph_integer_t ecount;

    igraph_rng_seed(igraph_rng_default(), 137);

    /* Null graph (0 vertices) */
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    check_all_agree(&g, NULL, IGRAPH_ALL);
    igraph_destroy(&g);

    /* Singleton graph (1 vertex, no edges) */
    igraph_empty(&g, 1, IGRAPH_DIRECTED);
    check_all_agree(&g, NULL, IGRAPH_OUT);
    igraph_destroy(&g);

    /* Edgeless directed graph */
    igraph_empty(&g, 5, IGRAPH_DIRECTED);
    check_all_agree(&g, NULL, IGRAPH_OUT);
    igraph_destroy(&g);

    /* Undirected ring with non-negative weights */
    igraph_ring(&g, 6, IGRAPH_UNDIRECTED, /*mutual=*/ false, /*circular=*/ true);
    igraph_vector_init_int(&weights, 6, 2, 5, 1, 3, 4, 1);
    check_all_agree(&g, &weights, IGRAPH_ALL);
    check_all_agree(&g, NULL, IGRAPH_ALL);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Directed graph with self-loops and parallel edges, non-negative weights.
     * Uses the same topology as the Johnson and Floyd-Warshall tests. */
    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0,1, 0,2, 1,1, 1,2, 1,3, 2,0, 2,3, 3,4, 3,4,
                 -1);
    igraph_vector_init_int(&weights, 9, 1, 4, 0, 2, 3, 1, 0, 5, 5);
    check_all_agree(&g, &weights, IGRAPH_OUT);
    check_all_agree(&g, &weights, IGRAPH_IN);
    check_all_agree(&g, NULL, IGRAPH_OUT);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Random sparse directed graph with positive weights */
    igraph_erdos_renyi_game_gnp(&g, 50, 0.12,
        IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    ecount = igraph_ecount(&g);
    igraph_vector_init(&weights, ecount);
    for (igraph_int_t i = 0; i < ecount; i++) {
        /* Integer weights are represented exactly as doubles, ensuring that all
         * algorithms produce bit-identical distance matrices. */
        VECTOR(weights)[i] = (igraph_real_t) RNG_INTEGER(1, 10);
    }
    check_all_agree(&g, &weights, IGRAPH_OUT);
    check_all_agree(&g, NULL, IGRAPH_OUT);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Random sparse undirected graph with positive weights */
    igraph_erdos_renyi_game_gnp(&g, 50, 0.12,
        IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    ecount = igraph_ecount(&g);
    igraph_vector_init(&weights, ecount);
    for (igraph_int_t i = 0; i < ecount; i++) {
        VECTOR(weights)[i] = (igraph_real_t) RNG_INTEGER(1, 10);
    }
    check_all_agree(&g, &weights, IGRAPH_ALL);
    check_all_agree(&g, NULL, IGRAPH_ALL);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
