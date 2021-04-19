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
    igraph_t g, g_test;
    igraph_bool_t same;
    long int n_vertices = 4;

    /*    Undirected, should be a full graph    */
    IGRAPH_ASSERT(igraph_full_citation(&g, n_vertices, 0 /*undirected*/) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 4, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    /*    Directed, only edges from i->j if i > j    */
    IGRAPH_ASSERT(igraph_full_citation(&g, n_vertices, 1 /*directed*/) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 4, IGRAPH_DIRECTED, 1, 0, 2, 0, 3, 0, 2, 1, 3, 1, 3, 2, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    /*    Directed, 1 vertex, should be edgeless    */
    IGRAPH_ASSERT(igraph_full_citation(&g, 1 /*n_vertices*/, 1 /*directed*/) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    /*    Directed, 0 vertices, empty graph    */
    IGRAPH_ASSERT(igraph_full_citation(&g, 0 /*n_vertices*/, 1 /*directed*/) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;

}
