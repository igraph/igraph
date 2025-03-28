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

#include <stdio.h>
#include <igraph.h>
#include "test_utilities.h"

void null_graph(void) {
    igraph_t graph;
    igraph_vector_t result;
    igraph_vector_int_t vertexOrder;

    igraph_vector_init(&result, 0);
    igraph_vector_int_init(&vertexOrder, 0);
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED); // null graph

    /* output */
    printf("Null graph\n");
    igraph_rich_club_density_sequence(&graph, &vertexOrder, 0, 0, 0, &result);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void singleton_graph(void) {
    igraph_t graph;
    igraph_vector_t result;
    igraph_vector_int_t vertexOrder;

    igraph_vector_init(&result, 1);
    igraph_vector_int_init(&vertexOrder, 1); // vertexOrder: [0]
    VECTOR(vertexOrder)[0] = 0;
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED); // singleton

    /* output */
    printf("Test 2: singleton graph\n");
    igraph_rich_club_density_sequence(&graph, &vertexOrder, 0, 0, 0, &result);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void more_complex_graph(void) {
    igraph_t graph;
    igraph_vector_t result;
    igraph_vector_int_t vertexOrder;
    const igraph_integer_t numVertices = 7;

    igraph_vector_init(&result, numVertices);
    igraph_vector_int_init(&vertexOrder, numVertices);
    for (int i = 0; i < numVertices; i++){
        VECTOR(vertexOrder)[i] = i;
    }
    igraph_small(&graph, numVertices, IGRAPH_UNDIRECTED,
                 0,3, 1,3, 2,3, 4,3, 5,3, 5,6, 1,2, 2,5, -1);

    /* output */
    printf("Test 3a: more complex graph (in-order vertex removal)\n"); // vertexOrder: in order 0-6
    igraph_rich_club_density_sequence(&graph, &vertexOrder, 0, 0, 0, &result);
    print_vector(&result);
    printf("\n");

    igraph_vector_int_reverse(&vertexOrder);
    printf("Test 3b: more complex graph (reverse vertex removal)\n"); // vertexOrder: reverse 6-0
    igraph_rich_club_density_sequence(&graph, &vertexOrder, 0, 0, 0, &result);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

int main(void) {
    null_graph();      // (NaN)
    printf("\n");

    singleton_graph(); // (NaN)
    printf("\n");

    /* in-order: (0.380952 0.466667 0.5 0.5 0.333333 1 NaN)
     * reverse: (0.380952 0.466667 0.5 0.666667 0.333333 1 NaN)
     */
    more_complex_graph();

    VERIFY_FINALLY_STACK();
    return 0;
}
