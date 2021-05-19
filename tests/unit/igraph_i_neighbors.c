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
#include "../../src/graph/neighbors.h"

void call_and_print(igraph_t *graph, igraph_integer_t pnode, igraph_neimode_t mode, igraph_loops_t loops, igraph_multiple_t multiple) {
    igraph_vector_t neis;
    igraph_vector_init(&neis, 0);
    IGRAPH_ASSERT(igraph_i_neighbors(graph, &neis, pnode, mode, loops, multiple) == IGRAPH_SUCCESS);
    print_vector(&neis);
    igraph_vector_destroy(&neis);
}


int main() {
    igraph_t g_1, g_lm, g_lmu, g_s1, g_s2;

    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_s1, 2, 1, 0,1, 0,1, 1,0, 1,0, -1);
    igraph_small(&g_s2, 2, 1, 0,1, 1,0, 1,0, -1);

    igraph_vector_t neis;
    igraph_vector_init(&neis, 0);

    printf("One vertex:\n");
    call_and_print(&g_1, 0, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Vertex with multiple edges, IGRAPH_IN, IGRAPH_MULTIPLE:\n");
    call_and_print(&g_lm, 0, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Vertex with multiple edges, IGRAPH_OUT, IGRAPH_MULTIPLE:\n");
    call_and_print(&g_lm, 0, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Vertex with multiple edges, IGRAPH_ALL, IGRAPH_MULTIPLE:\n");
    call_and_print(&g_lm, 0, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("Vertex with multiple edges, undirected, IGRAPH_MULTIPLE:\n");
    call_and_print(&g_lmu, 0, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Vertex with multiple edges, IGRAPH_IN, IGRAPH_NO_MULTIPLE:\n");
    call_and_print(&g_lm, 0, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);

    printf("Vertex with multiple edges, IGRAPH_OUT, IGRAPH_NO_MULTIPLE:\n");
    call_and_print(&g_lm, 0, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);

    printf("Vertex with multiple edges, IGRAPH_ALL, IGRAPH_NO_MULTIPLE:\n");
    call_and_print(&g_lm, 0, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    printf("Vertex with multiple edges, undirected, IGRAPH_NO_MULTIPLE:\n");
    call_and_print(&g_lmu, 0, IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);

    printf("Vertex 1 with loop, IGRAPH_OUT, IGRAPH_NO_LOOPS:\n");
    call_and_print(&g_lm, 1, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);

    printf("Vertex 1 with loop, IGRAPH_ALL, IGRAPH_NO_LOOPS:\n");
    call_and_print(&g_lm, 1, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);

    printf("Vertex 1 with loop, undirected, IGRAPH_NO_LOOPS:\n");
    call_and_print(&g_lmu, 1, IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);

    printf("Vertex 1 with loop, IGRAPH_OUT, IGRAPH_LOOPS_ONCE:\n");
    call_and_print(&g_lm, 1, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Vertex 1 with loop, IGRAPH_ALL, IGRAPH_LOOPS_ONCE:\n");
    call_and_print(&g_lm, 1, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Vertex 1 with loop, undirected, IGRAPH_LOOPS_ONCE:\n");
    call_and_print(&g_lmu, 1, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Vertex 1 with loop, IGRAPH_ALL, IGRAPH_LOOPS_TWICE:\n");
    call_and_print(&g_lm, 1, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("Vertex 1 with loop, undirected, IGRAPH_LOOPS_TWICE:\n");
    call_and_print(&g_lmu, 1, IGRAPH_IN, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("Graph with 2 edges from 0 to 1, and 2 from 1 to 0, IGRAPH_ALL, IGRAPH_MULTIPLE:\n");
    call_and_print(&g_s1, 0, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("Graph with 1 edge from 0 to 1, and 2 from 1 to 0, IGRAPH_ALL, IGRAPH_MULTIPLE:\n");
    call_and_print(&g_s2, 0, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Trying IGRAPH_LOOPS_TWICE with IGRAPH_OUT:\n");
    IGRAPH_ASSERT(igraph_i_neighbors(&g_lm, &neis, 0, IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE) == IGRAPH_EINVAL);

    printf("Trying invalid vertex ID:\n");
    IGRAPH_ASSERT(igraph_i_neighbors(&g_lm, &neis, 42, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE) == IGRAPH_EINVVID);

    printf("Trying invalid mode:\n");
    IGRAPH_ASSERT(igraph_i_neighbors(&g_lm, &neis, 0, (igraph_neimode_t) 42, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE) == IGRAPH_EINVMODE);

    igraph_destroy(&g_1);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lmu);
    igraph_destroy(&g_s1);
    igraph_destroy(&g_s2);
    igraph_vector_destroy(&neis);

    VERIFY_FINALLY_STACK();
    return 0;
}
