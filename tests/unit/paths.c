/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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
    igraph_t graph;
    igraph_vector_int_t edge_walk, vertex_walk;

    igraph_vector_int_init(&vertex_walk, 0);
    igraph_vector_int_init(&edge_walk, 0);

    igraph_ring(&graph, 4, IGRAPH_DIRECTED, false, true);
    igraph_vector_int_range(&edge_walk, 0, igraph_ecount(&graph));

    printf("Cycle, detect start:\n");
    igraph_vertex_path_from_edge_path(&graph, -1, &edge_walk, &vertex_walk, IGRAPH_ALL);
    print_vector_int(&edge_walk);
    print_vector_int(&vertex_walk);

    printf("\nCycle plus one vertex, start given:\n");
    igraph_vector_int_push_back(&edge_walk, 0);
    igraph_vertex_path_from_edge_path(&graph, 0, &edge_walk, &vertex_walk, IGRAPH_ALL);
    print_vector_int(&edge_walk);
    print_vector_int(&vertex_walk);

    /* Inconsistent start vertex */
    CHECK_ERROR(
        igraph_vertex_path_from_edge_path(&graph, 1, &edge_walk, &vertex_walk, IGRAPH_ALL),
        IGRAPH_EINVAL
    );

    /* Invalid start vertex */
    CHECK_ERROR(
        igraph_vertex_path_from_edge_path(&graph, 10, &edge_walk, &vertex_walk, IGRAPH_ALL),
        IGRAPH_EINVVID
    );

    printf("\nSingle edge, detect start:\n");
    igraph_vector_int_clear(&edge_walk);
    igraph_vector_int_push_back(&edge_walk, 2);
    igraph_vertex_path_from_edge_path(&graph, -1, &edge_walk, &vertex_walk, IGRAPH_ALL);
    print_vector_int(&edge_walk);
    print_vector_int(&vertex_walk);

    printf("\nSingle edge, directed OUT, detect start:\n");
    igraph_vector_int_clear(&edge_walk);
    igraph_vector_int_push_back(&edge_walk, 2);
    igraph_vertex_path_from_edge_path(&graph, -1, &edge_walk, &vertex_walk, IGRAPH_OUT);
    print_vector_int(&edge_walk);
    print_vector_int(&vertex_walk);

    printf("\nSingle edge, directed IN, detect start:\n");
    igraph_vector_int_clear(&edge_walk);
    igraph_vector_int_push_back(&edge_walk, 2);
    igraph_vertex_path_from_edge_path(&graph, -1, &edge_walk, &vertex_walk, IGRAPH_IN);
    print_vector_int(&edge_walk);
    print_vector_int(&vertex_walk);

    printf("\nZero edges, start given:\n");
    igraph_vector_int_clear(&edge_walk);
    igraph_vertex_path_from_edge_path(&graph, 2, &edge_walk, &vertex_walk, IGRAPH_ALL);
    print_vector_int(&edge_walk);
    print_vector_int(&vertex_walk);

    /* Zero edges, start not given: */
    CHECK_ERROR(
        igraph_vertex_path_from_edge_path(&graph, -1, &edge_walk, &vertex_walk, IGRAPH_ALL),
        IGRAPH_EINVAL
    );

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&edge_walk);
    igraph_vector_int_destroy(&vertex_walk);

    VERIFY_FINALLY_STACK();

    return 0;
}
