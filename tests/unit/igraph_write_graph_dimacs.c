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

    igraph_t g;
    igraph_vector_t capacity;
    int source = 0;
    int target = 5;

    /*
    Expected output:

    ```
    c created by igraph
    p problem n_vertices n_edges
    n source s
    n target t
    a arc_node_1 arc_node2 capacity
    ```

    We always outout max as the problem.

    */

    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 3, 4, 3, 5, 4, 5, -1);
    igraph_vector_init_int_end(&capacity, -1, 5, 2, 2, 3, 4, 1, 2, 5, -1);

    printf("DIMACS graph output:\n");
    igraph_write_graph_dimacs(&g, stdout, source, target, &capacity);

    igraph_destroy(&g);
    igraph_vector_destroy(&capacity);

    igraph_small(&g, 0, IGRAPH_DIRECTED, -1);
    igraph_vector_init(&capacity, 0);

    /* Check that the function does not crash/misbehave on a null graph.
       Note that currently igraph outputs DIMACS files for the max-flow
       problem, which only makes sense if there are at least two vertices,
       a source and the target. Here we use dummy values for them. */
    printf("\nDIMACS graph output for null graph:\n");
    igraph_write_graph_dimacs(&g, stdout, source, target, &capacity);

    igraph_vector_destroy(&capacity);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
