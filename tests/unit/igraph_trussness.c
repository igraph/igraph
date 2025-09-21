/*

  Copyright 2017 The Johns Hopkins University Applied Physics Laboratory LLC. All Rights Reserved.
  Copyright 2021 The igraph team

  Truss algorithm for cohesive subgroups.

  Author: Alex Perrone
  Date: 2017-08-03
  Minor edits: The igraph team, 2021

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA

*/

#include <stdio.h>
#include <igraph.h>

#include "test_utilities.h"

void print_and_destroy(igraph_t *graph, igraph_vector_int_t *trussness) {
    igraph_int_t i, n = igraph_vector_int_size(trussness);

    printf("fromNode, toNode, trussness\n");
    for (i=0; i < n; i++) {
        igraph_int_t from, to;
        igraph_edge(graph, i, &from, &to);
        printf("%" IGRAPH_PRId ", %" IGRAPH_PRId ", %" IGRAPH_PRId "\n", from, to, VECTOR(*trussness)[i]);
    }

    igraph_vector_int_destroy(trussness);
    igraph_destroy(graph);
}

int main(void) {

    igraph_t graph;
    igraph_vector_int_t trussness;

    /* Create actual graph */
    igraph_int_t edges[] = { 0,1, 0,2, 0,3, 0,4,
      1,2, 1,3, 1,4, 2,3, 2,4, 3,4, 3,6, 3,11,
      4,5, 4,6, 5,6, 5,7, 5,8, 5,9, 6,7, 6,10, 6,11,
      7,8, 7,9, 8,9, 8,10 };
    const igraph_vector_int_t v = igraph_vector_int_view(edges, sizeof(edges) / sizeof(edges[0]));

    igraph_create(&graph, &v, 0, IGRAPH_UNDIRECTED);

    /* Compute the trussness of the edges. */
    printf("Simple graph:\n");
    igraph_vector_int_init(&trussness, 0);
    igraph_trussness(&graph, &trussness);
    print_and_destroy(&graph, &trussness);

    /* Add some loop edges -- they should have trussness = 2 */
    printf("\nGraph with loops:\n");
    igraph_create(&graph, &v, 0, IGRAPH_UNDIRECTED);
    igraph_add_edge(&graph, 0, 0);
    igraph_add_edge(&graph, 7, 7);
    igraph_add_edge(&graph, 5, 5);
    igraph_vector_int_init(&trussness, 0);
    igraph_trussness(&graph, &trussness);
    print_and_destroy(&graph, &trussness);

    /* Null graph trivial case */
    printf("\nNull graph:\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&trussness, 0);
    igraph_trussness(&graph, &trussness);
    print_and_destroy(&graph, &trussness);

    /* Singleton graph trivial case */
    printf("\nSingleton graph:\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&trussness, 0);
    igraph_trussness(&graph, &trussness);
    print_and_destroy(&graph, &trussness);

    /* Graph with no edges trivial case */
    printf("\nGraph with no edges:\n");
    igraph_empty(&graph, 10, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&trussness, 0);
    igraph_trussness(&graph, &trussness);
    print_and_destroy(&graph, &trussness);

    VERIFY_FINALLY_STACK();

    /* Multigraph */
    printf("\nTrying multigraph:\n");
    igraph_create(&graph, &v, 0, IGRAPH_UNDIRECTED);
    igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_MUTUAL);
    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_EACH, 0);
    igraph_vector_int_init(&trussness, 0);
    CHECK_ERROR(igraph_trussness(&graph, &trussness), IGRAPH_UNIMPLEMENTED);
    igraph_vector_int_destroy(&trussness);
    igraph_destroy(&graph);

    return 0;
}
