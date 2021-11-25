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
    igraph_vector_int_t l;
    igraph_bool_t iso;
    igraph_integer_t i;

    /* Testing typical use case */
    /* Test n = 5, undirected, l = [1, 2] */
    igraph_vector_int_init(&l, 0);
    igraph_vector_int_push_back(&l, 1);
    igraph_vector_int_push_back(&l, 2);

    IGRAPH_ASSERT(igraph_circulant(&graph, 5, &l, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 2,
            2, 4, 4, 1, 1, 3, 3, 0, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&l);

    /* Testing simplification */
    /* n = 6, undirected, l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] */
    igraph_vector_int_init(&l, 0);
    for (i = 1; i <= 15; i++) {
        igraph_vector_int_push_back(&l, i);
    }
    IGRAPH_ASSERT(igraph_circulant(&graph, 6, &l, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_full(&graph_test, 6, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&l);

    // todo: add more tests

    VERIFY_FINALLY_STACK();
    return 0;
}
