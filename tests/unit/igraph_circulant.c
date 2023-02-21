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

#include "test_utilities.h"

int main(void) {

    igraph_t graph, graph_test;
    igraph_vector_int_t shifts;
    igraph_bool_t iso, same;
    igraph_integer_t i;

    /* Testing invalid values */

    /* n = -3, shifts = [1], undirected */
    igraph_vector_int_init_int(&shifts, 1, 1);
    CHECK_ERROR(igraph_circulant(&graph, -3, &shifts, IGRAPH_UNDIRECTED), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&shifts);

    /* Testing n = 0 case */
    /* n = 0, shifts = [1], undirected */
    igraph_vector_int_init_int(&shifts, 1, 1);
    IGRAPH_ASSERT(igraph_circulant(&graph, 0, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 0);
    IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
    igraph_vector_int_destroy(&shifts);
    igraph_destroy(&graph);

    /* Testing n = 1 case */
    /* n = 1, shifts = [1], undirected */
    igraph_vector_int_init_int(&shifts, 1, 1);
    IGRAPH_ASSERT(igraph_circulant(&graph, 1, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 1, IGRAPH_UNDIRECTED, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&graph, &graph_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing n = 2 case */
    /* n = 2, shifts = [1], undirected */
    igraph_vector_int_init_int(&shifts, 1, 1);
    IGRAPH_ASSERT(igraph_circulant(&graph, 2, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 2, IGRAPH_UNDIRECTED, 0, 1, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&graph, &graph_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing empty list case */
    /* n = 5, shifts = [], undirected */
    igraph_vector_int_init(&shifts, 0);
    IGRAPH_ASSERT(igraph_circulant(&graph, 5, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 5, IGRAPH_UNDIRECTED, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&graph, &graph_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing typical use case */
    /* Test n = 5, shifts = [1, 2], undirected */
    igraph_vector_int_init_int(&shifts, 2, 1, 2);
    IGRAPH_ASSERT(igraph_circulant(&graph, 5, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 2,
            2, 4, 4, 1, 1, 3, 3, 0, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing simplification */
    /* n = 6, shifts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], undirected */
    igraph_vector_int_init(&shifts, 0);
    for (i = 1; i <= 15; i++) {
        igraph_vector_int_push_back(&shifts, i);
    }
    IGRAPH_ASSERT(igraph_circulant(&graph, 6, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_full(&graph_test, 6, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing simplification with negative values */
    /* n = 7, shifts = [-15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5,
            -4, -3, -2, -1, 0, 1, 2, 3, 4], undirected */
    igraph_vector_int_init(&shifts, 0);
    for (i = -15; i <= 4; i++) {
        igraph_vector_int_push_back(&shifts, i);
    }
    IGRAPH_ASSERT(igraph_circulant(&graph, 7, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_full(&graph_test, 7, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing simplification when n is even, the offset = n/2, and the graph is undirected */
    /* n = 6, shifts = [3], undirected */
    igraph_vector_int_init_int(&shifts, 1, 3);
    IGRAPH_ASSERT(igraph_circulant(&graph, 6, &shifts, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 6, IGRAPH_UNDIRECTED, 0, 3, 1, 4, 2, 5, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing directed graph */
    /* n = 7, shifts = [1, -2], directed */
    igraph_vector_int_init_int(&shifts, 2, 1, -2);
    IGRAPH_ASSERT(igraph_circulant(&graph, 7, &shifts, IGRAPH_DIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 7, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 0,
            0, 5, 1, 6, 2, 0, 3, 1, 4, 2, 5, 3, 6, 4, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing simplification works for directed too */
    /* n = 5, shifts = [-4, 1, 6], directed */
    igraph_vector_int_init_int(&shifts, 3, -4, 1, 6);
    IGRAPH_ASSERT(igraph_circulant(&graph, 5, &shifts, IGRAPH_DIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, -1);
    IGRAPH_ASSERT(igraph_isomorphic(&graph, &graph_test, &iso) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* Testing that things that are normally simplified for undirected are not simplified for directed */
    /* n = 5, shifts = [1, -1], directed */
    igraph_vector_int_init_int(&shifts, 2, 1, -1);
    IGRAPH_ASSERT(igraph_circulant(&graph, 5, &shifts, IGRAPH_DIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0,
            1, 0, 2, 1, 3, 2, 4, 3, 0, 4, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&graph, &graph_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    /* n = 6, shifts = [3], directed */
    igraph_vector_int_init_int(&shifts, 1, 3);
    IGRAPH_ASSERT(igraph_circulant(&graph, 6, &shifts, IGRAPH_DIRECTED) == IGRAPH_SUCCESS);
    igraph_small(&graph_test, 6, IGRAPH_DIRECTED, 0, 3, 1, 4, 2, 5, 3, 0, 4, 1, 5, 2, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&graph, &graph_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph);
    igraph_destroy(&graph_test);
    igraph_vector_int_destroy(&shifts);

    VERIFY_FINALLY_STACK();
    return 0;
}
