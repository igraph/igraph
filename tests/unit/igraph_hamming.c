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

#include <igraph.h>

#include "test_utilities.h"

/*
 * This is an alternative Hamming graph generator based on the fact that
 * the Hamming graph H(n,q) is, equivalently, the cartesian product of n
 * copies of complete graphs K(q).
 */
igraph_error_t cartesian_hamming(igraph_t *graph, igraph_int_t n, igraph_int_t q) {
    igraph_t full;

    if (n < 0 || q < 0) {
        IGRAPH_ERROR("n and q must not be negative.", IGRAPH_EINVAL);
    }

    igraph_full(&full, q, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_empty(graph, 1, IGRAPH_UNDIRECTED);

    for (igraph_int_t i = 0; i < n; i++) {
        igraph_t temp;
        igraph_product(&temp, graph, &full, IGRAPH_PRODUCT_CARTESIAN);
        igraph_destroy(graph);
        igraph_copy(graph, &temp);
        igraph_destroy(&temp);
    }

    igraph_destroy(&full);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_t g1, g2, g3;
    igraph_bool_t iso;

    /* Compare H_n_1 to the singleton graph K_1. */
    igraph_hamming(&g1, 3, 1, IGRAPH_UNDIRECTED); // d = 3
    IGRAPH_ASSERT(igraph_vcount(&g1) == 1);
    IGRAPH_ASSERT(igraph_ecount(&g1) == 0);
    IGRAPH_ASSERT(!igraph_is_directed(&g1));
    igraph_destroy(&g1);

    /* Compare H_1_q to the complete(full) graph K_q. */
    igraph_full(&g2, 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_hamming(&g1, 1, 10, IGRAPH_UNDIRECTED);
    igraph_isomorphic(&g1, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    /* Compare H_n_2 to the hypercube graph Q_n. */
    /* undirected */
    igraph_hypercube(&g2, 10, IGRAPH_UNDIRECTED);
    igraph_hamming(&g1, 10, 2, IGRAPH_UNDIRECTED);
    igraph_isomorphic(&g2, &g1, &iso);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g1);
    igraph_destroy(&g2);
    /* directed */
    igraph_hypercube(&g2, 10, IGRAPH_DIRECTED);
    igraph_hamming(&g1, 10, 2, IGRAPH_DIRECTED);
    igraph_isomorphic(&g2, &g1, &iso);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    /* Compare H_2_q to the lattice(rook) graph L_q_q. */
    igraph_full(&g3, 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_product(&g2, &g3, &g3, IGRAPH_PRODUCT_CARTESIAN);
    igraph_hamming(&g1, 2, 10, IGRAPH_UNDIRECTED);
    igraph_isomorphic(&g2, &g1, &iso);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g1);
    igraph_destroy(&g2);
    igraph_destroy(&g3);

    /* Test edge cases with n==0 or q==0 */

    igraph_empty(&g2, 1, IGRAPH_UNDIRECTED); /* singleton */
    igraph_empty(&g3, 0, IGRAPH_UNDIRECTED); /* null graph */

    igraph_hamming(&g1, 0, 0, IGRAPH_UNDIRECTED);
    igraph_isomorphic(&g1, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g1);

    igraph_hamming(&g1, 0, 2, IGRAPH_UNDIRECTED);
    igraph_isomorphic(&g1, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g1);

    igraph_hamming(&g1, 1, 0, IGRAPH_UNDIRECTED);
    igraph_isomorphic(&g1, &g3, &iso);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g1);

    igraph_destroy(&g3);
    igraph_destroy(&g2);

    /* Compare to the implementation in terms of the graph product of cliques. */
    for (igraph_int_t n = 0; n < 4; n++) {
        for (igraph_int_t q = 0; q < 4; q++) {
            igraph_hamming(&g1, n, q, IGRAPH_UNDIRECTED);
            cartesian_hamming(&g2, n, q);
            igraph_isomorphic(&g1, &g2, &iso);
            IGRAPH_ASSERT(iso);
            igraph_destroy(&g2);
            igraph_destroy(&g1);
        }
    }

    CHECK_ERROR(igraph_hamming(&g1, -1, 0, IGRAPH_UNDIRECTED), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_hamming(&g1, -1, 10, IGRAPH_UNDIRECTED), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_hamming(&g1, 10, -1, IGRAPH_UNDIRECTED), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();

    return 0;
}
