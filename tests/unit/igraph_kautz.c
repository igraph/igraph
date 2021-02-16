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
    igraph_bool_t iso, same;

    /*   BA, AB, CB, etc.     */
    IGRAPH_ASSERT(igraph_kautz(&g, /* m */ 2, /* n */ 1) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 6, IGRAPH_DIRECTED, 0, 1, 0, 5, 1, 0, 1, 4, 2, 0, 2, 4, 3, 1, 3, 5, 4, 2, 4, 3, 5, 3, 5, 2, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&g, &g_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    /*  1 symbol, string length 11, should be empty graph   */
    IGRAPH_ASSERT(igraph_kautz(&g, /* m */ 0, /* n */ 10) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 0, IGRAPH_DIRECTED, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    /*  1 symbol, string length 1 should be single vertex   */
    IGRAPH_ASSERT(igraph_kautz(&g, /* m */ 0, /* n */ 0) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 1, IGRAPH_DIRECTED, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    /*  String length 1 should be full graph   */
    IGRAPH_ASSERT(igraph_kautz(&g, /* m */ 5, /* n */ 0) == IGRAPH_SUCCESS);
    igraph_full(&g_test, 6, IGRAPH_DIRECTED, /*loops*/ 0);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    VERIFY_FINALLY_STACK();
    return 0;
}
