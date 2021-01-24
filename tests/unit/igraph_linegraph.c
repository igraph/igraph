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
    igraph_t g_start, g_line, g_test;
    igraph_bool_t same;

    /*    Undirected    */
    igraph_small(&g_start, 7, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 1, 3, 1, 3, 2, 2, 2, 4, 3, 4, 4, 5, -1);
    IGRAPH_ASSERT(igraph_linegraph(&g_start, &g_line) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 8, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 1, 4, 1, 4, 1, 5, 2, 3, 2, 3,
                 2, 6, 3, 6, 4, 5, 4, 5, 5, 6, 5, 7, 6, 7, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g_line, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g_start);
    igraph_destroy(&g_line);
    igraph_destroy(&g_test);

     /*    Directed    */
    igraph_small(&g_start, 7, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 1, 3, 3, 1, 2, 2, 2, 4, 3, 4, 4, 5, -1);
    IGRAPH_ASSERT(igraph_linegraph(&g_start, &g_line) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 8, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 1, 4, 1, 5, 2, 3, 2, 6, 3, 1, 3, 2, 4, 4, 4, 5,
                 5, 7, 6, 7, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g_line, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g_start);
    igraph_destroy(&g_line);
    igraph_destroy(&g_test);

    /*    No edges    */
    igraph_small(&g_start, 7, IGRAPH_DIRECTED, -1);
    IGRAPH_ASSERT(igraph_linegraph(&g_start, &g_line) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 0, IGRAPH_DIRECTED, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g_line, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g_start);
    igraph_destroy(&g_line);
    igraph_destroy(&g_test);

    VERIFY_FINALLY_STACK();
    return 0;
}
