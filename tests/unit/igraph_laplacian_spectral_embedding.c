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

int main(void) {
    igraph_t graph;
    igraph_matrix_t X;

    /* Every vertex has positive in- and out-degree, as required by OAP. */
    igraph_small(&graph, 8, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4,
                 4, 5, 5, 6, 6, 7, 7, 0,
                 0, 2, 1, 4, 3, 6, 5, 7, 7, 3,
                 -1);
    igraph_matrix_init(&X, 0, 0);
    igraph_rng_seed(igraph_rng_default(), 42);

    /* The second half of the directed embedding is optional. */
    igraph_laplacian_spectral_embedding(&graph, 2, NULL,
                                        IGRAPH_EIGEN_LA, IGRAPH_EMBEDDING_OAP,
                                        /*scaled=*/ false, &X, /*Y=*/ NULL,
                                        /*D=*/ NULL, /*options=*/ NULL);
    IGRAPH_ASSERT(igraph_matrix_nrow(&X) == 8);
    IGRAPH_ASSERT(igraph_matrix_ncol(&X) == 2);

    igraph_matrix_destroy(&X);
    igraph_destroy(&graph);
    VERIFY_FINALLY_STACK();

    return 0;
}
