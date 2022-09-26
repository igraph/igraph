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

    igraph_t g;
    igraph_vector_int_t center;

    igraph_vector_int_init(&center, 0);

    printf("Null graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nSingleton graph:\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nPath with isolated vertex:\n");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED,
                 0,2,
                 -1);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nDisconnected graph:\n");
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, -1);
    igraph_graph_center(&g, &center, IGRAPH_OUT);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nUndirected path graph:\n");
    igraph_ring(&g, 5, IGRAPH_UNDIRECTED, /* mutual */ 0, /* circular */ 0);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nUndirected Graph\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 0, 2, 3, 0, 4, 1, 2, 5, -1);
    igraph_graph_center(&g, &center, IGRAPH_OUT);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nDirected path graph:\n");
    igraph_ring(&g, 5, IGRAPH_DIRECTED, /* mutual */ 0, /* circular */ 0);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nUndirected star:\n");
    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nOut-star:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_ALL) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nIn-star:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    IGRAPH_ASSERT(igraph_graph_center(&g, &center, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_vector_int(&center);
    igraph_destroy(&g);

    igraph_vector_int_destroy(&center);

    VERIFY_FINALLY_STACK();

    return 0;
}
