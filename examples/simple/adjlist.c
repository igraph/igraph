/*
   igraph library.

   This example shows how to build and use an adjacency list in igraph.
   An adjacency list stores, for each vertex, the vertices that are
   directly connected to it by edges.

   Adjacency lists are useful when you want to iterate over the neighbors
   of vertices efficiently.

   Copyright (C) 2008–2012  Gábor Csárdi
   Copyright (C) 2025      igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <igraph.h>
#include <stdio.h>

int main(void) {
    igraph_t graph;
    igraph_adjlist_t adjacency;
    igraph_integer_t vertex;
    igraph_integer_t neighbor_index;

    /* Initialize the igraph library. */
    igraph_setup();

    /* Create a small directed graph with five vertices.
       The edges are:
         0 -> 1
         0 -> 2
         1 -> 2
         2 -> 3
         3 -> 4
    */
    igraph_small(
        &graph, 5, IGRAPH_DIRECTED,
        0, 1,
        0, 2,
        1, 2,
        2, 3,
        3, 4,
        -1
    );

    /* Build an adjacency list containing the outgoing neighbors
       of each vertex. */
    igraph_adjlist_init(
        &graph,
        &adjacency,
        IGRAPH_OUT,        /* outgoing edges */
        IGRAPH_NO_LOOPS,   /* ignore self-loops */
        IGRAPH_NO_MULTIPLE /* ignore parallel edges */
    );

    /* Print the adjacency list. */
    printf("Adjacency list (outgoing neighbors):\n");
    for (vertex = 0; vertex < igraph_vcount(&graph); vertex++) {
        const igraph_vector_int_t *neighbors =
            igraph_adjlist_get(&adjacency, vertex);

        printf("Vertex %" IGRAPH_PRId ": ", vertex);
        for (neighbor_index = 0;
             neighbor_index < igraph_vector_int_size(neighbors);
             neighbor_index++) {
            printf("%" IGRAPH_PRId " ",
                   VECTOR(*neighbors)[neighbor_index]);
        }
        printf("\n");
    }

    /* Free allocated data structures. */
    igraph_adjlist_destroy(&adjacency);
    igraph_destroy(&graph);

    return 0;
}
