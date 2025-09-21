/*
   igraph library.
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

void validate_sums(igraph_matrix_t* m, igraph_bool_t column_wise) {
    igraph_vector_t v;
    igraph_vector_t expected;

    igraph_vector_init(&v, 0);
    igraph_vector_init(&expected, 0);

    if (column_wise) {
        igraph_matrix_colsum(m, &v);
    } else {
        igraph_matrix_rowsum(m, &v);
    }

    igraph_vector_resize(&expected, igraph_matrix_nrow(m));
    igraph_vector_fill(&expected, 1);
    VECTOR(expected)[igraph_vector_size(&expected) - 1] = 0;

    IGRAPH_ASSERT(igraph_vector_all_almost_e(&v, &expected, 1e-7));

    igraph_vector_destroy(&expected);
    igraph_vector_destroy(&v);
}

void test_graph(igraph_bool_t directed) {
    igraph_t graph;
    igraph_real_t weights_array[] = { 5, 4, 3, 2, 1, 6, 3, 2 };
    const igraph_vector_t weights =
        igraph_vector_view(weights_array, sizeof(weights_array) / sizeof(weights_array[0]));
    igraph_matrix_t m;
    const char* prefix = directed ? "Directed" : "Undirected";

    igraph_small(
        &graph, 6, directed ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED,
        0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 3, 2, 2, 0, 1, -1
    );

    igraph_matrix_init(&m, 2, 2);

    printf("%s, unweighted, row-wise:\n", prefix);
    igraph_get_stochastic(&graph, &m, /* column_wise = */ 0, NULL);
    igraph_matrix_print(&m);
    validate_sums(&m, /* column_wise = */ 0);
    printf("========\n");

    printf("%s, unweighted, column-wise:\n", prefix);
    igraph_get_stochastic(&graph, &m, /* column_wise = */ 1, NULL);
    igraph_matrix_print(&m);
    validate_sums(&m, /* column_wise = */ 1);
    printf("========\n");

    printf("%s, weighted, row-wise:\n", prefix);
    igraph_get_stochastic(&graph, &m, /* column_wise = */ 0, &weights);
    igraph_matrix_print(&m);
    validate_sums(&m, /* column_wise = */ 0);
    printf("========\n");

    printf("%s, weighted, column-wise:\n", prefix);
    igraph_get_stochastic(&graph, &m, /* column_wise = */ 1, &weights);
    igraph_matrix_print(&m);
    validate_sums(&m, /* column_wise = */ 1);
    printf("========\n");

    igraph_matrix_destroy(&m);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

int main(void) {

    test_graph(/* directed = */ 0);
    test_graph(/* directed = */ 1);

    return 0;
}
