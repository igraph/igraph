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
    igraph_t res;
    igraph_t expected_res;
    igraph_bool_t is_isomorph;

    //should be a vertex
    igraph_mycielski_graph(&res, 1);
    igraph_small(
        &expected_res, 1, /* directed = */ 0,
        -1
    );
    igraph_isomorphic(&res, &expected_res, &is_isomorph);
    IGRAPH_ASSERT(is_isomorph);
    igraph_destroy(&res);
    igraph_destroy(&expected_res);

    // should be a path
    igraph_mycielski_graph(&res, 2);
    igraph_small(
        &expected_res, 2, /* directed = */ 0,
        /* edge 0 */ 0, 1,
        -1
    );
    igraph_isomorphic(&res, &expected_res, &is_isomorph);
    IGRAPH_ASSERT(is_isomorph);
    igraph_destroy(&res);
    igraph_destroy(&expected_res);

    // should a 5-cycle
    igraph_mycielski_graph(&res, 3);
    igraph_small(
        &expected_res, 5, /* directed = */ 0,
        /* edge 0 */ 0, 1,
        /* edge 1 */ 1, 2,
        /* edge 2 */ 2, 3,
        /* edge 3 */ 3, 4,
        /* edge 4 */ 4, 0,
        -1
    );
    igraph_isomorphic(&res, &expected_res, &is_isomorph);
    IGRAPH_ASSERT(is_isomorph);
    igraph_destroy(&res);
    igraph_destroy(&expected_res);

    // should be a grotzsch graph
    igraph_mycielski_graph(&res, 4);
    igraph_famous(&expected_res, "grotzsch");
    igraph_isomorphic(&res, &expected_res, &is_isomorph);
    IGRAPH_ASSERT(is_isomorph);
    igraph_destroy(&res);
    igraph_destroy(&expected_res);

    VERIFY_FINALLY_STACK();

    return 0;
}
