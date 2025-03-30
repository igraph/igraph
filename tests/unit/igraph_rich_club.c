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
    printf("Test 1: null graph\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, false, false);
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
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, false, false);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void undirected_no_loop_graph(void) {
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
    printf("Test 3a: undirected, no-loop graph (in-order vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, false, false);
    print_vector(&result);
    printf("\n");

    igraph_vector_int_reverse(&vertexOrder);
    printf("Test 3b: undirected, no-loop graph (reverse vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, false, false);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void directed_no_loop_graph(void) {
    igraph_t graph;
    igraph_vector_t result;
    igraph_vector_int_t vertexOrder;
    const igraph_integer_t numVertices = 7;

    igraph_vector_init(&result, numVertices);
    igraph_vector_int_init(&vertexOrder, numVertices);
    for (int i = 0; i < numVertices; i++){
        VECTOR(vertexOrder)[i] = i;
    }
    igraph_small(&graph, numVertices, IGRAPH_DIRECTED,
                 0,2, 1,2, 2,3, 1,3, 3,5, 3,4, 5,6, 6,5, -1);

    /* output */
    printf("Test 4a: directed, no-loop graph (in-order vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, false, true);
    print_vector(&result);
    printf("\n");

    igraph_vector_int_reverse(&vertexOrder);
    printf("Test 4b: directed, no-loop graph (reverse vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, false, true);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void undirected_loop_graph(void) {
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
                 0,3, 1,3, 2,3, 4,4, 5,3, 5,6, 1,2, 2,5, 6,4, -1);

    /* output */
    printf("Test 5a: undirected, loop graph (in-order vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, true, false);
    print_vector(&result);
    printf("\n");

    igraph_vector_int_reverse(&vertexOrder);
    printf("Test 5b: undirected, loop graph (reverse vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, true, false);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void directed_loop_graph(void) {
    igraph_t graph;
    igraph_vector_t result;
    igraph_vector_int_t vertexOrder;
    const igraph_integer_t numVertices = 7;

    igraph_vector_init(&result, numVertices);
    igraph_vector_int_init(&vertexOrder, numVertices);
    for (int i = 0; i < numVertices; i++){
        VECTOR(vertexOrder)[i] = i;
    }
    igraph_small(&graph, numVertices, IGRAPH_DIRECTED,
                 0,2, 1,2, 2,3, 1,3, 3,5, 3,4, 5,6, 6,5, 4,4, -1);

    /* output */
    printf("Test 6a: directed, loop graph (in-order vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, true, true);
    print_vector(&result);
    printf("\n");

    igraph_vector_int_reverse(&vertexOrder);
    printf("Test 6b: directed, loop graph (reverse vertex removal)\n");
    igraph_rich_club_density_sequence(&graph, NULL, &result, &vertexOrder, true, true);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void weighted_graph(void) {
    // same undirected, no-loop graph as Test 3
    igraph_t graph;
    igraph_vector_t result, weights;
    igraph_vector_int_t vertexOrder;
    const igraph_integer_t numVertices = 7;

    igraph_small(&graph, numVertices, IGRAPH_UNDIRECTED,
                 0,3, 1,3, 2,3, 4,3, 5,3, 5,6, 1,2, 2,5, -1);

    igraph_vector_init(&result, numVertices);
    igraph_vector_int_init(&vertexOrder, numVertices);
    for (int i = 0; i < numVertices; i++){
        VECTOR(vertexOrder)[i] = i;
    }

    igraph_vector_init(&weights, igraph_ecount(&graph));
    for (int i = 0; i < igraph_ecount(&graph); i++) {
        // all "edges" have weight 2, output should be double of Test 3a
        VECTOR(weights)[i] = 2;
    }

    /* output */
    printf("Test 7a: weighted graph\n");
    igraph_rich_club_density_sequence(&graph, &weights, &result,
                                      &vertexOrder, false, false);
    print_vector(&result);
    printf("\n");

    // change weights to include some non-integer weights
    // for this specific graph, weights should now be the same as that of Test 3a
    VECTOR(weights)[0] = 1;   // (0,3)
    VECTOR(weights)[1] = 0.5; // (1,3)
    VECTOR(weights)[2] = 0.5; // (2,3)
    VECTOR(weights)[3] = 1;   // (4,3)
    VECTOR(weights)[4] = 1;   // (5,3)
    VECTOR(weights)[5] = 1;   // (5,6)
    VECTOR(weights)[6] = 1.5; // (1,2)
    VECTOR(weights)[7] = 1.5; // (2,5)

    printf("Test 7b: weighted graph (non-integer weights)\n");
    igraph_rich_club_density_sequence(&graph, &weights, &result,
                                      &vertexOrder, false, false);
    print_vector(&result);

    igraph_vector_int_destroy(&vertexOrder);
    igraph_vector_destroy(&result);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
}

int main(void) {
    null_graph();
    printf("\n");

    singleton_graph();
    printf("\n");

    undirected_no_loop_graph();
    printf("\n");

    directed_no_loop_graph();
    printf("\n");

    undirected_loop_graph();
    printf("\n");

    directed_loop_graph();
    printf("\n");

    weighted_graph();

    VERIFY_FINALLY_STACK();
    return 0;
}
