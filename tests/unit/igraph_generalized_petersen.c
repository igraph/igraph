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

    /* Compares G(5,2) with Petersen Graph */
    IGRAPH_ASSERT(igraph_generalized_petersen(&graph, /* n */ 5, /* k */ 2) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 10, IGRAPH_UNDIRECTED, 0, 2, 0, 3, 0, 5, 1, 3, 1, 4, 1, 6, 2, 4, 2, 7, 3, 8, 4, 9, 5, 6, 5, 9, 6, 7, 7, 8, 8, 9, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);

    /* Compares G(5,2) with Petersen Graph in igraph_famous */
    IGRAPH_ASSERT(igraph_generalized_petersen(&graph, /* n */ 5, /* k */ 2) == IGRAPH_SUCCESS);
    igraph_famous(&graph_test, "petersen");
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);

    VERIFY_FINALLY_STACK();
    return 0;
}
