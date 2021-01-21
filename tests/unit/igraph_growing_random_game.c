/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.inc"


int main() {
    igraph_t g;
    igraph_vector_t degree;
    igraph_bool_t tree;
    int i;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* no vertices */

    igraph_growing_random_game(&g, /* n: vertices */ 0, /* m: edges_per_vertex */ 3, /* directed */ 0, /* citation */ 0);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    IGRAPH_ASSERT(!igraph_is_directed(&g));
    igraph_destroy(&g);

    /* 1 edge per vertex with citation makes a tree */

    igraph_growing_random_game(&g, /* n: vertices */ 20, /* m: edges_per_vertex */ 1, /* directed */ 0, /* citation */ 1);
    igraph_is_tree(&g, &tree, /* root*/ NULL, /*unused mode*/ 0);
    IGRAPH_ASSERT(tree);
    igraph_destroy(&g);

    /* out degree of citation equals edges per vertex */

    igraph_growing_random_game(&g, /* n: vertices */ 10, /* m: edges_per_vertex */ 7, /* directed */ 1, /* citation */ 1);
    igraph_vector_init(&degree, 0);
    igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    for(i = 1; i < 10; i++) {
        IGRAPH_ASSERT(VECTOR(degree)[i] == 7);
    }
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_vector_destroy(&degree);
    igraph_destroy(&g);

    /* total number of edges is (vertices - 1) * edges */

    igraph_growing_random_game(&g, /* n: vertices */ 10, /* m: edges_per_vertex */ 7, /* directed */ 1, /* citation */ 0);
    IGRAPH_ASSERT(igraph_ecount(&g) == 63);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
