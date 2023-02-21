/* -*- mode: C -*-  */
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

void call_and_print(igraph_t *g, igraph_vector_t *weights, igraph_bool_t unconn, igraph_bool_t directed) {
    igraph_vector_int_t path_vertex, path_edge;
    igraph_real_t result;
    igraph_integer_t from, to;

    igraph_vector_int_init(&path_edge, 0);
    igraph_vector_int_init(&path_vertex, 0);

    igraph_diameter_dijkstra(g, weights, &result, &from, &to, &path_vertex, &path_edge, directed, unconn);

    printf("Diameter: ");
    print_real(stdout, result, "%g");
    printf(", from %" IGRAPH_PRId " to %" IGRAPH_PRId "\n", from, to);
    printf("Edges:\n");
    print_vector_int(&path_edge);
    printf("Vertices:\n");
    print_vector_int(&path_vertex);

    igraph_vector_int_destroy(&path_edge);
    igraph_vector_int_destroy(&path_vertex);
    printf("\n");
}

int main(void) {
    igraph_t g_ring, g_0, g_1, g_2, g_lm;
    igraph_vector_t weights, weights_neg, weights_0;

    igraph_vector_init_int(&weights, 9, 1, 2, 3, 4, 5, 1, 1, 1, 1);
    igraph_vector_init_int(&weights_neg, 9, -1, 2, 3, 4, 5, 1, 1, 1, 1);
    igraph_vector_init(&weights_0, 0);

    igraph_empty(&g_0, 0, IGRAPH_DIRECTED);
    igraph_empty(&g_1, 1, IGRAPH_DIRECTED);
    igraph_empty(&g_2, 2, IGRAPH_DIRECTED);
    igraph_ring(&g_ring, 10, IGRAPH_DIRECTED, 0, 0);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 4,3, 4,3, -1);

    printf("Graph with zero nodes:\n");
    call_and_print(&g_0, NULL, 1, 1);

    printf("Graph with zero nodes, weighted:\n");
    call_and_print(&g_0, &weights_0, 1, 1);

    printf("Graph with one node:\n");
    call_and_print(&g_1, NULL, 1, 1);

    printf("Graph with one node, weighted:\n");
    call_and_print(&g_1, &weights_0, 1, 1);

    printf("Graph with one node, returns inf for unconnected:\n");
    call_and_print(&g_1, NULL, 0, 1);

    printf("Graph with two nodes:\n");
    call_and_print(&g_2, NULL, 1, 1);

    printf("Graph with two nodes, returns inf for unconnected:\n");
    call_and_print(&g_2, NULL, 0, 1);

    printf("Ring without weights:\n");
    call_and_print(&g_ring, NULL, 1, 1);

    printf("Ring with weights:\n");
    call_and_print(&g_ring, &weights, 1, 1);

    printf("Graph with loops and multiple edges:\n");
    call_and_print(&g_lm, NULL, 1, 1);

    printf("Graph with loops and multiple edges, direction ignored:\n");
    call_and_print(&g_lm, NULL, 1, 0);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Ring with some negative weights:\n");
    igraph_diameter_dijkstra(&g_ring, &weights_neg, NULL, NULL, NULL, NULL, NULL, IGRAPH_DIRECTED, 1);

    printf("Ring with wrong weight vector size:\n");
    igraph_diameter_dijkstra(&g_ring, &weights_0, NULL, NULL, NULL, NULL, NULL, IGRAPH_DIRECTED, 1);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_2);
    igraph_destroy(&g_ring);
    igraph_destroy(&g_lm);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&weights_neg);
    igraph_vector_destroy(&weights_0);

    VERIFY_FINALLY_STACK();
    return 0;
}
