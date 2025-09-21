/* igraph library.  Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

void call_and_print(igraph_t *graph, igraph_int_t n, igraph_wheel_mode_t mode,
                igraph_int_t center) {

    IGRAPH_ASSERT(igraph_wheel(graph, n, mode, center) == IGRAPH_SUCCESS);
    print_graph_canon(graph);
    igraph_destroy(graph);
    printf("\n");
}

int main(void) {
    igraph_t graph;

    printf("-- Test graph with 1 vertex --\n");
    call_and_print(&graph, 1, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("-- Test graph with 2 vertices --\n");
    call_and_print(&graph, 2, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("-- Test graph with OUT mode --\n");
    call_and_print(&graph, 4, IGRAPH_WHEEL_OUT, 0);
    printf("-- Test graph with IN mode --\n");
    call_and_print(&graph, 4, IGRAPH_WHEEL_IN, 0);
    printf("-- Test graph with MUTUAL mode --\n");
    call_and_print(&graph, 4, IGRAPH_WHEEL_MUTUAL, 0);
    printf("-- Test graph with UNDIRECTED mode \n");
    call_and_print(&graph, 4, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("-- Test graph with center equal to n/2 --\n");
    call_and_print(&graph, 4, IGRAPH_WHEEL_OUT, 2);
    printf("-- Test graph with center equal to n - 1: --\n");
    call_and_print(&graph, 4, IGRAPH_WHEEL_OUT, 3);
    printf("-- Test graph with center equal to 1: --\n");
    call_and_print(&graph, 4, IGRAPH_WHEEL_OUT, 1);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
