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
    igraph_t g, g_copy;
    igraph_bool_t same;
    igraph_vector_t degrees;
    igraph_vs_t vertices;
    igraph_rng_seed(igraph_rng_default(), 42);

    /*No edges, should just return the same graph*/
    igraph_small(&g, 5, IGRAPH_DIRECTED, -1);
    IGRAPH_ASSERT(igraph_rewire_directed_edges(&g, /*probability*/ 0.1, /*loops*/ 0, /*mode*/ IGRAPH_ALL) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    IGRAPH_ASSERT(igraph_vcount(&g) == 5);
    igraph_destroy(&g);

    /*No rewire*/
    igraph_small(&g, 10, IGRAPH_DIRECTED, 0, 1, 0, 3, 5, 4, 4, 8, 9, 2, 9, 3, 9, 7, 7, 7, 7, 8, -1);
    igraph_copy(&g_copy, &g);
    IGRAPH_ASSERT(igraph_rewire_directed_edges(&g, /*probability*/ 0.0, /*loops*/ 0, /*mode*/ IGRAPH_ALL) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_copy, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);

    /*Rewire*/
    IGRAPH_ASSERT(igraph_rewire_directed_edges(&g, /*probability*/ 0.5, /*loops*/ 1, /*mode*/ IGRAPH_ALL) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_copy, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(!same); /*guaranteed for this seed*/
    IGRAPH_ASSERT(igraph_ecount(&g) == 9);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    igraph_destroy(&g);
    igraph_destroy(&g_copy);

    /*Out-star remains out-star if outs are moved*/
    igraph_small(&g, 10, IGRAPH_DIRECTED, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, -1);
    IGRAPH_ASSERT(igraph_rewire_directed_edges(&g, /*probability*/ 1.0, /*loops*/ 0, /*mode*/ IGRAPH_OUT) == IGRAPH_SUCCESS);
    igraph_vector_init(&degrees, 0);
    igraph_vs_1(&vertices, 0);
    igraph_degree(&g, &degrees, vertices, IGRAPH_ALL, 0);
    IGRAPH_ASSERT(VECTOR(degrees)[0] == 9);
    igraph_vector_destroy(&degrees);
    igraph_vs_destroy(&vertices);

    /*A few erroneous calls*/
    igraph_set_error_handler(igraph_error_handler_ignore);
    IGRAPH_ASSERT(igraph_rewire_directed_edges(&g, /*probability*/ -0.1, /*loops*/ 0, /*mode*/ IGRAPH_ALL) == IGRAPH_EINVAL);
    IGRAPH_ASSERT(igraph_rewire_directed_edges(&g, /*probability*/ 1.1, /*loops*/ 0, /*mode*/ IGRAPH_ALL) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
