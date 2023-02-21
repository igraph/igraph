/*
   IGraph library.
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
    /*most functionality is currently tested in
      igraph_get_shortest_paths_bellman_ford(), which is called by
      igraph_get_shortest_path_bellman_ford()*/
    igraph_t g;
    igraph_vector_int_t vertices, edges;

    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);

    printf("Basic example, don't ask for vertices and edges.\n");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    igraph_get_shortest_path_bellman_ford(&g, NULL, NULL,
            0, 1, NULL, IGRAPH_OUT);

    printf("Basic example, ask for vertices and edges:\n");
    igraph_get_shortest_path_bellman_ford(&g, &vertices, &edges,
            0, 1, NULL, IGRAPH_OUT);
    printf("vertices: ");
    print_vector_int(&vertices);
    printf("edges: ");
    print_vector_int(&edges);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    printf("Check error when passing null graph.\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    CHECK_ERROR(
        igraph_get_shortest_path_bellman_ford(&g, NULL, NULL, 0, 0, NULL, IGRAPH_OUT),
        IGRAPH_EINVVID);

    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&edges);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
