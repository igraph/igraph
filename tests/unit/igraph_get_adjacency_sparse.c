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

void print_sparse_matrix(const igraph_sparsemat_t* sm) {
    igraph_matrix_t m;

    igraph_matrix_init(&m, 0, 0);
    igraph_sparsemat_as_matrix(&m, sm);
    igraph_matrix_print(&m);
    igraph_matrix_destroy(&m);
}


void test_undirected(void) {
    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    igraph_vector_t weights;
    igraph_sparsemat_t m;

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&graph));

    igraph_sparsemat_init(&m, 2, 2, 0);

    printf("Undirected, unweighted, upper, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, NULL, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, upper, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, NULL, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, upper, loops twice:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, NULL, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, lower, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, NULL, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, lower, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, NULL, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, lower, loops twice:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, NULL, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, both, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, both, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, unweighted, both, loops twice:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, upper, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, &weights, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, upper, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, &weights, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, upper, loops twice:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_UPPER, &weights, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, lower, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, &weights, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, lower, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, &weights, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, lower, loops twice:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_LOWER, &weights, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, both, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, both, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Undirected, weighted, both, loops twice:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    igraph_sparsemat_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

void test_directed(void) {
    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    igraph_vector_t weights;
    igraph_sparsemat_t m;

    igraph_small(&graph, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&graph));

    igraph_sparsemat_init(&m, 2, 2, 0);

    printf("Directed, unweighted, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Directed, unweighted, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Directed, unweighted, loops twice (same as once):\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Directed, weighted, no loops:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_NO_LOOPS);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Directed, weighted, loops once:\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_LOOPS_ONCE);
    print_sparse_matrix(&m);
    printf("========\n");

    printf("Directed, weighted, loops twice (same as once):\n");
    igraph_get_adjacency_sparse(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_LOOPS_TWICE);
    print_sparse_matrix(&m);
    printf("========\n");

    igraph_sparsemat_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

void test_errors(void) {
    igraph_t graph;
    igraph_sparsemat_t m;

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_sparsemat_init(&m, 2, 2, 0);

    CHECK_ERROR(igraph_get_adjacency_sparse(&graph, &m, (igraph_get_adjacency_t) 42, NULL, IGRAPH_LOOPS_ONCE), IGRAPH_EINVAL);

    igraph_sparsemat_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

int main(void) {

    test_undirected();
    test_directed();
    test_errors();

    return 0;
}
