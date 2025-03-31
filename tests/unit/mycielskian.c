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

    // k == 0 testing
    igraph_small(
        &g, 5, /* directed = */ 0,
        /* edge 0 */ 0, 1,
        /* edge 1 */ 0, 3,
        /* edge 2 */ 1, 2,
        /* edge 3 */ 2, 4,
        /* edge 4 */ 2, 3,
        /* edge 5 */ 3, 4,
        -1
    );
    igraph_mycielskian(&res, &g, 0);
    print_graph(&res);
    printf("\n\n");
    igraph_destroy(&g);
    igraph_destroy(&res);

    // vcount == 0, k==0 testing
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&res, &g, 0);
    print_graph(&res);
    printf("\n\n");
    igraph_destroy(&g);
    igraph_destroy(&res);

    // vcount == 0, k==1 testing
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&res, &g, 1);
    print_graph(&res);
    printf("\n\n");
    igraph_destroy(&g);
    igraph_destroy(&res);

    // vcount == 0, k==3 testing
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&res, &g, 3);
    print_graph(&res);
    printf("\n\n");
    igraph_destroy(&g);
    igraph_destroy(&res);

    // vcount == 1, k==0 testing
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&res, &g, 0);
    print_graph(&res);
    printf("\n\n");
    igraph_destroy(&g);
    igraph_destroy(&res);

    // checking for directed graphs for directed attribute
    igraph_small(
        &g, 5, /* directed = */ 1,
        /* edge 0 */ 0, 1,
        /* edge 1 */ 0, 3,
        /* edge 2 */ 1, 2,
        /* edge 3 */ 2, 4,
        /* edge 4 */ 2, 3,
        /* edge 5 */ 3, 4,
        -1
    );
    igraph_mycielskian(&res, &g, 0);
    print_graph(&res);
    printf("\n\n");
    igraph_destroy(&g);
    igraph_destroy(&res);

    // vcount == 1, k==3 testing --> should output a Grötzsch graph
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_mycielskian(&res, &g, 3);
    igraph_t grostzsch;
    igraph_bool_t res_bool;
    igraph_famous(&grostzsch, "grotzsch");
    igraph_isomorphic(&grostzsch, &res, &res_bool);
    IGRAPH_ASSERT(res_bool);
    igraph_destroy(&grostzsch);
    igraph_destroy(&g);
    igraph_destroy(&res);

    // path graph gives k=2 a grotzsch graph
    igraph_small(
        &g, 2, /* directed = */ 0,
        /* edge 0 */ 0, 1,
        -1
    );
    igraph_mycielskian(&res, &g, 2);
    igraph_famous(&grostzsch, "grotzsch");
    igraph_isomorphic(&grostzsch, &res, &res_bool);
    IGRAPH_ASSERT(res_bool);
    igraph_destroy(&grostzsch);
    igraph_destroy(&g);
    igraph_destroy(&res);

    VERIFY_FINALLY_STACK();

    return 0;
}
