/* IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

    printf("Graph with no vertices:\n");
    igraph_small(&g, IGRAPH_UNDIRECTED, 0, -1);
    igraph_connect_neighborhood(&g, 10, IGRAPH_ALL);
    print_graph(&g);
    igraph_destroy(&g);

    printf("Directed graph with loops and multiple edges, order 0, IGRAPH_OUT:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    igraph_connect_neighborhood(&g, 0, IGRAPH_OUT);
    print_graph(&g);
    printf("Same graph, order 1:\n");
    igraph_connect_neighborhood(&g, 1, IGRAPH_OUT);
    print_graph(&g);
    printf("Same graph, order 2:\n");
    igraph_connect_neighborhood(&g, 2, IGRAPH_OUT);
    print_graph(&g);
    igraph_destroy(&g);

    printf("Same starting graph, order 2, IGRAPH_IN:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    igraph_connect_neighborhood(&g, 2, IGRAPH_IN);
    print_graph(&g);
    igraph_destroy(&g);

    printf("Same starting graph, order 2, IGRAPH_ALL:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    igraph_connect_neighborhood(&g, 2, IGRAPH_ALL);
    print_graph(&g);
    igraph_destroy(&g);

    printf("Same starting graph, order 3, IGRAPH_OUT:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    igraph_connect_neighborhood(&g, 3, IGRAPH_OUT);
    print_graph(&g);
    igraph_destroy(&g);

    printf("Same starting graph, order 12, IGRAPH_IN:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    igraph_connect_neighborhood(&g, 12, IGRAPH_IN);
    print_graph(&g);
    igraph_destroy(&g);

    printf("Same starting graph, but undirected, order 3, IGRAPH_OUT:\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    igraph_connect_neighborhood(&g, 2, IGRAPH_OUT);
    print_graph(&g);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    printf("Check negative order error.\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    CHECK_ERROR(igraph_connect_neighborhood(&g, -1, IGRAPH_OUT), IGRAPH_EINVAL);
    igraph_destroy(&g);

    return 0;
}
