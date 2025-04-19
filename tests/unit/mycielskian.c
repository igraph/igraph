/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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
    igraph_t g, res;
    igraph_t groetzsch;
    igraph_bool_t isomorphic;

    igraph_famous(&groetzsch, "groetzsch");

    // k == 0 testing
    printf("small graph, k=0\n");
    igraph_small(
        &g, 5, /* directed = */ false,
        /* edge 0 */ 0, 1,
        /* edge 1 */ 0, 3,
        /* edge 2 */ 1, 2,
        /* edge 3 */ 2, 4,
        /* edge 4 */ 2, 3,
        /* edge 5 */ 3, 4,
        -1
    );
    igraph_mycielskian(&g, &res, 0);
    print_graph(&res);
    igraph_destroy(&res);
    igraph_destroy(&g);
    printf("\n\n");

    // vcount == 0, k==0 testing
    printf("null graph, k=0\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&g, &res, 0);
    print_graph(&res);
    igraph_destroy(&res);
    igraph_destroy(&g);
    printf("\n\n");

    // vcount == 0, k==1 testing
    printf("null graph, k=1\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&g, &res, 1);
    print_graph(&res);
    igraph_destroy(&res);
    igraph_destroy(&g);
    printf("\n\n");

    // vcount == 0, k==3 testing
    printf("null graph, k=3\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&g, &res, 3);
    print_graph(&res);
    igraph_destroy(&res);
    igraph_destroy(&g);
    printf("\n\n");

    // vcount == 1, k==0 testing
    printf("singleton graph, k=0\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&g, &res, 0);
    print_graph(&res);
    igraph_destroy(&res);
    igraph_destroy(&g);
    printf("\n\n");

    // checking for directed graphs for directed attribute
    printf("small directed graph, k=0\n");
    igraph_small(
        &g, 5, /* directed = */ true,
        /* edge 0 */ 0, 1,
        /* edge 1 */ 0, 3,
        /* edge 2 */ 1, 2,
        /* edge 3 */ 2, 4,
        /* edge 4 */ 2, 3,
        /* edge 5 */ 3, 4,
        -1
    );
    igraph_mycielskian(&g, &res, 0);
    print_graph(&res);
    igraph_destroy(&res);
    igraph_destroy(&g);
    printf("\n\n");

    // vcount == 1, k==3 testing --> should output a Grötzsch graph
    printf("singleton graph, k=3\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&g, &res, 3);
    print_graph(&res);
    igraph_isomorphic(&groetzsch, &res, &isomorphic);
    IGRAPH_ASSERT(isomorphic);
    igraph_destroy(&res);
    igraph_destroy(&g);
    printf("\n\n");

    // P_2 graph with k=2 gives the Grötzsch graph
    printf("P_2, k=2\n");
    igraph_ring(&g, 2, IGRAPH_UNDIRECTED, false, false);
    igraph_mycielskian(&g, &res, 2);
    print_graph(&res);
    igraph_isomorphic(&groetzsch, &res, &isomorphic);
    IGRAPH_ASSERT(isomorphic);
    igraph_destroy(&res);
    igraph_destroy(&g);

    igraph_destroy(&groetzsch);

    VERIFY_FINALLY_STACK();

    return 0;
}
