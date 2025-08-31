/*
   igraph library.
   Copyright (C) 2007-2024  The igraph development team <igraph@igraph.org>

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

void analyze_loops(const igraph_t *graph) {
    igraph_vector_bool_t is_loop;
    igraph_bool_t has_loop;
    igraph_int_t loop_count;

    igraph_has_loop(graph, &has_loop);
    printf("Has loops? %s\n", has_loop ? "Yes" : "No");

    igraph_count_loops(graph, &loop_count);
    printf("How many? %" IGRAPH_PRId "\n", loop_count);

    igraph_vector_bool_init(&is_loop, 0);
    igraph_is_loop(graph, &is_loop, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    printf("Loop positions: "); igraph_vector_bool_print(&is_loop);
    igraph_vector_bool_destroy(&is_loop);

    printf("\n");
}

int main(void) {
    igraph_t graph;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,1, 0,1, 1,0, 3,4, 11,10, -1);
    analyze_loops(&graph);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,0, 1,1, 2,2, 2,3, 2,4, 2,5, 2,6, 2,2, 0,0, -1);
    analyze_loops(&graph);
    igraph_destroy(&graph);

    return 0;
}
