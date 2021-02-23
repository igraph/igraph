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
    /*igraph_subisomorphic now calls igraph_subisomorphic_vf2, most of
      the testing is done through calling that directly*/
    igraph_t g1, g2;
    igraph_bool_t result;

    printf("No vertices.\n");
    igraph_small(&g1, 0, 0, -1);
    igraph_small(&g2, 0, 0, -1);
    IGRAPH_ASSERT(igraph_subisomorphic(&g1, &g2, &result) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(result);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("Basic positive example, undirected.\n");
    igraph_small(&g1, 4, 0, 0,1, 1,2, 3,2, 3,1, -1);
    igraph_small(&g2, 3, 0, 0,1, 1,2, 2,0, -1);
    IGRAPH_ASSERT(igraph_subisomorphic(&g1, &g2, &result) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(result);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("Basic negative example, directed.\n");
    igraph_small(&g1, 4, 1, 0,1, 1,2, 3,2, 3,1, -1);
    igraph_small(&g2, 3, 1, 0,1, 1,2, 2,0, -1);
    IGRAPH_ASSERT(igraph_subisomorphic(&g1, &g2, &result) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(!result);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Mismatching directedness.\n");
    igraph_small(&g1, 4, 1, 0,1, 1,2, 3,2, 3,1, -1);
    igraph_small(&g2, 3, 0, 0,1, 1,2, 2,0, -1);
    IGRAPH_ASSERT(igraph_subisomorphic(&g1, &g2, &result) == IGRAPH_EINVAL);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    VERIFY_FINALLY_STACK();
    return 0;
}
