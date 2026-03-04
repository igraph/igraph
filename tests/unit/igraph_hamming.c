/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.h"

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

    CHECK_ERROR(igraph_hamming(&g1, 0, 0, IGRAPH_UNDIRECTED), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_hamming(&g1, 0, 10, IGRAPH_UNDIRECTED), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_hamming(&g1, 10, 0, IGRAPH_UNDIRECTED), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();

    return 0;
}
