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

void test_undirected_with_subsets(void) {


    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    igraph_vector_t weights;
    igraph_matrix_t m;
    igraph_vs_t from;
    igraph_vs_t to;

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&graph));
    igraph_vs_vector_small(&from, 0, 1, 2, -1);
    igraph_vs_vector_small(&to, 0, 2, -1);

    igraph_matrix_init(&m, 2, 2);

    printf("Undirected, unweighted, no loops, from [0,1,2] to [0,2]:\n");
    igraph_get_adjacency_submatrix(&graph, &m, from, to, NULL, IGRAPH_NO_LOOPS);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, unweighted, loops once, from [0,1,2] to [0,2]:\n");
    igraph_get_adjacency_submatrix(&graph, &m, from, to, NULL, IGRAPH_LOOPS_ONCE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, unweighted, loops twice, from [0,1,2] to [0,2]:\n");
    igraph_get_adjacency_submatrix(&graph, &m, from, to, NULL, IGRAPH_LOOPS_TWICE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, no loops, from [0,1,2] to [0,2]:\n");
    igraph_get_adjacency_submatrix(&graph, &m, from, to, &weights, IGRAPH_NO_LOOPS);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, loops once, from [0,1,2] to [0,2]:\n");
    igraph_get_adjacency_submatrix(&graph, &m, from, to, &weights, IGRAPH_LOOPS_ONCE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, loops twice, from [0,1,2] to [0,2]:\n");
    igraph_get_adjacency_submatrix(&graph, &m, from, to, &weights, IGRAPH_LOOPS_TWICE);
    igraph_matrix_print(&m);
    printf("========\n");

    igraph_matrix_destroy(&m);
    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

}

void test_undirected(void) {
    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    igraph_vector_t weights;
    igraph_matrix_t m;

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&graph));

    igraph_matrix_init(&m, 2, 2);

    printf("Undirected, unweighted, no loops:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_NO_LOOPS);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, unweighted, loops once:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_LOOPS_ONCE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, unweighted, loops twice:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_LOOPS_TWICE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, no loops:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_NO_LOOPS);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, loops once:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_LOOPS_ONCE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Undirected, weighted, loops twice:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_LOOPS_TWICE);
    igraph_matrix_print(&m);
    printf("========\n");

    igraph_matrix_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

void test_directed(void) {
    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    igraph_vector_t weights;
    igraph_matrix_t m;

    igraph_small(&graph, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&graph));

    igraph_matrix_init(&m, 2, 2);

    printf("Directed, unweighted, no loops:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_NO_LOOPS);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Directed, unweighted, loops once:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_LOOPS_ONCE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Directed, unweighted, loops twice (same as once):\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_LOOPS_TWICE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Directed, weighted, no loops:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_NO_LOOPS);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Directed, weighted, loops once:\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_LOOPS_ONCE);
    igraph_matrix_print(&m);
    printf("========\n");

    printf("Directed, weighted, loops twice (same as once):\n");
    igraph_get_adjacency_submatrix(&graph, &m, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_LOOPS_TWICE);
    igraph_matrix_print(&m);
    printf("========\n");

    igraph_matrix_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

int main(void) {

    test_undirected();
    test_undirected_with_subsets();
    test_directed();

    return 0;
}
