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

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_bool_t result;

    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    igraph_empty(&g, 2, IGRAPH_UNDIRECTED);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    igraph_small(&g, 2, IGRAPH_UNDIRECTED,
                 0,1,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(result);
    igraph_destroy(&g);

    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,3, 3,0,
                 2,4, 4,5, 5,2,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(result);
    igraph_destroy(&g);

    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 1,3,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 1,3, 3,4, 4,2,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(result);
    igraph_destroy(&g);

    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 1,3, 3,4,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    /* Two disjoint cycles */
    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 3,4, 4,5, 5,3,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    /* Cycle + isolated vertex */
    igraph_small(&g, 4, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    /* Special case: the root is an articulation point */
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 0,3, 3,4, 4,0,
                 -1);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    /* Cache concistency checks */

    /* K_2 is a special graph: it is biconnected yet it has no cycles.
     * Make sure we don't accidentally set or interpret the cache
     * incorrectly. */

    igraph_bool_t acyclic;

    igraph_full(&g, 2, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(result);
    igraph_is_acyclic(&g, &acyclic);
    IGRAPH_ASSERT(acyclic);

    igraph_invalidate_cache(&g);

    igraph_is_acyclic(&g, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_is_biconnected(&g, &result);
    IGRAPH_ASSERT(result);

    igraph_destroy(&g);


    return 0;
}
