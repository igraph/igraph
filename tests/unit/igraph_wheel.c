/* IGraph library.  Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

void call_and_print(igraph_t *graph, igraph_integer_t n, igraph_wheel_mode_t mode,
                igraph_integer_t center) {

    IGRAPH_ASSERT(igraph_wheel(graph, n, mode, center) == IGRAPH_SUCCESS);
    print_graph_canon(graph);
    igraph_destroy(graph);
    printf("\n");
}

int main() {
    igraph_t g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8;

    printf("-- Test graph with 1 vertex --\n");
    call_and_print(&g_1, 1, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("-- Test graph with 2 vertices --\n");
    call_and_print(&g_2, 2, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("-- Test graph with OUT mode --\n");
    call_and_print(&g_3, 4, IGRAPH_WHEEL_OUT, 0);
    printf("-- Test graph with IN mode --\n");
    call_and_print(&g_4, 4, IGRAPH_WHEEL_IN, 0);
    printf("-- Test graph with MUTUAL mode --\n");
    call_and_print(&g_5, 4, IGRAPH_WHEEL_MUTUAL, 0);
    printf("-- Test graph with UNDIRECTED mode \n");
    call_and_print(&g_6, 4, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("-- Test graph with center equal to n/2 --\n");
    call_and_print(&g_7, 4, IGRAPH_WHEEL_OUT, 2);
    printf("-- Test graph with center equal n - 1: --\n");
    call_and_print(&g_8, 4, IGRAPH_WHEEL_OUT, 3);

    igraph_destroy(&g_1);
    igraph_destroy(&g_2);
    igraph_destroy(&g_3);
    igraph_destroy(&g_4);
    igraph_destroy(&g_5);
    igraph_destroy(&g_6);
    igraph_destroy(&g_7);
    igraph_destroy(&g_8);

    VERIFY_FINALLY_STACK();

    return 0;
}
