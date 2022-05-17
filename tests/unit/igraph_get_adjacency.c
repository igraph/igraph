/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>

#include "test_utilities.h"

void test_undirected() {
    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    igraph_vector_t weights;
    igraph_matrix_t m;

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&graph));

    igraph_matrix_init(&m, 2, 2);

    printf("Undirected, unweighted, upper:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, NULL);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, unweighted, lower:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, NULL);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, unweighted, both:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, upper:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, &weights);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, lower:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, &weights);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, both:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights);
    igraph_matrix_print(&m);
    printf("========\n");

    igraph_matrix_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

void test_directed() {
    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    igraph_vector_t weights;
    igraph_matrix_t m;

    igraph_small(&graph, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&graph));

    igraph_matrix_init(&m, 2, 2);

    printf("Directed, unweighted:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Directed, weighted:\n");
    igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights);
    igraph_matrix_print(&m);
    printf("========\n");

    igraph_matrix_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

void test_errors() {
    igraph_t graph;
    igraph_matrix_t m;

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_matrix_init(&m, 2, 2);

    igraph_set_error_handler(&igraph_error_handler_ignore);
    IGRAPH_ASSERT(igraph_get_adjacency(&graph, &m, 42, NULL) == IGRAPH_EINVAL);
    igraph_set_error_handler(&igraph_error_handler_abort);

    igraph_matrix_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

int main() {

    test_undirected();
    test_directed();
    test_errors();

    return 0;
}
