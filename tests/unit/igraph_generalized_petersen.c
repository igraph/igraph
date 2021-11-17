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
    igraph_t graph, graph_test;
    igraph_bool_t iso;

    /* Compares G(5,2) with Petersen Graph in igraph_famous */
    IGRAPH_ASSERT(igraph_generalized_petersen(&graph, /* n */ 5, /* k */ 2) == IGRAPH_SUCCESS);
    igraph_famous(&graph_test, "petersen");
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);

    /* Compares G(4,2) with expected small graph */
    IGRAPH_ASSERT(igraph_generalized_petersen(&graph, /* n */ 4, /* k */ 2) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 8, IGRAPH_UNDIRECTED, 0, 2, 0, 4, 1, 3, 1, 5, 2, 6, 3, 7, 4, 5, 4, 7, 5, 6, 6, 7, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);

    /* Compares G(10,2) with Dodecahedral Graph in igraph_famous */
    IGRAPH_ASSERT(igraph_generalized_petersen(&graph, /* n */ 10, /* k */ 2) == IGRAPH_SUCCESS);
    igraph_famous(&graph_test, "dodecahedron");
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);

    VERIFY_FINALLY_STACK();
    return 0;
}
