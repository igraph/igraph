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
#include "test_utilities.h"

void check_graph(const igraph_t *graph, const igraph_vector_int_t *terminals, const igraph_vector_t *weights) {
    igraph_real_t value;
    igraph_vector_int_t tree_edges;

    igraph_vector_int_init(&tree_edges, 0);

    igraph_error_t x = igraph_steiner_dreyfus_wagner(graph, terminals, weights, &value, &tree_edges);
    IGRAPH_ASSERT(x == IGRAPH_SUCCESS);

    print_vector_int(&tree_edges);
    printf("Total Steiner tree weight: %g\n", value);

    /* Check total tree weight. */
    igraph_real_t value2 = 0.0;
    igraph_integer_t tree_size = igraph_vector_int_size(&tree_edges);
    if (weights) {
        for (igraph_integer_t i = 0; i < tree_size; i++) {
            value2 += VECTOR(*weights)[VECTOR(tree_edges)[i]];
        }
    } else {
        value2 = tree_size;
    }
    IGRAPH_ASSERT(value == value2);

    /* Check that the result is indeed a tree. */
    if (igraph_vector_int_size(terminals) > 1) {
        igraph_t tree;
        igraph_subgraph_from_edges(graph, &tree, igraph_ess_vector(&tree_edges), /* delete_vertices= */ true);

        igraph_bool_t is_tree;
        igraph_is_tree(&tree, &is_tree, NULL, IGRAPH_ALL);
        IGRAPH_ASSERT(is_tree);

        igraph_destroy(&tree);
    }

    igraph_vector_int_destroy(&tree_edges);
}

int main(void) {
    igraph_t g_null, g_k7, g_k6_k1, g_k7_n, g_k7_n1, g_k7_real, g_k7_non_simple, g_k1, g_k7_self_loop, g_10_9, g_new, g_new_1, g_new_2, g_new_3, g_new_4;
    igraph_vector_int_t terminals_null, terminals_k7, terminals_k6_k1, terminals_k7_real, terminals_10_9, terminals_1, terminals_2, terminals_3, terminals_8, terminals_7, terminals_7_special, terminals_square, terminals_square_2, terminals_fpw;
    igraph_vector_t weights_null, weights_k7, weights_k6_k1, weights_k7_n, weights_k7_n1, weights_k7_real, weights_k7_non_simple, weights_k7_loop, weights_square_2, weights_fp;

    /* Null graph */
    igraph_empty(&g_null, 0, 0);
    igraph_vector_init(&weights_null, 0);
    igraph_vector_int_init(&terminals_null, 0);

    /* K_7 complete graph with a specific edge ordering. */
    igraph_small(&g_k7, 7, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,
                 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
                 2, 3, 2, 4, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6,
                 5, 6,
                 -1);
    // 1 , 1
    igraph_vector_init_int(&weights_k7, 21,
                           2, 2, 2, 1, 1, 2,
                           2, 2, 2, 1, 2,
                           2, 2, 2, 1,
                           1, 2, 1,
                           2, 1,
                           1);

    igraph_vector_int_init_int(&terminals_k7, 4,
                               0, 1, 2, 3);

    /* K_7 non-complete graph with edges 0-4,0-5 removed */
    igraph_small(&g_k7_n, 7, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 6,
                 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
                 2, 3, 2, 4, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6,
                 5, 6,
                 -1);

    igraph_vector_init_int(&weights_k7_n, 19,
                           2, 2, 2, 2,
                           2, 2, 2, 1, 2,
                           2, 2, 2, 1,
                           1, 2, 1,
                           2, 1,
                           1);

    /* K_7 non-complete graph with edges 0-4 removed */
    igraph_small(&g_k7_n1, 7, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 6,
                 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
                 2, 3, 2, 4, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6,
                 5, 6,
                 -1);

    igraph_vector_init_int(&weights_k7_n1, 20,
                           2, 2, 2, 1, 2,
                           2, 2, 2, 1, 2,
                           2, 2, 2, 1,
                           1, 2, 1,
                           2, 1,
                           1);

    /* Disconnected graph: K_6 with a specific edge order and an additional isolated vertex. */
    igraph_small(&g_k6_k1, 7, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5,
                 1, 2, 1, 3, 1, 4, 1, 5,
                 2, 3, 2, 4, 2, 5,
                 3, 4, 3, 5,
                 4, 5,
                 -1);

    igraph_vector_init_int(&weights_k6_k1, 15,
                           2, 2, 2, 1, 1,
                           2, 2, 2, 1,
                           2, 2, 2,
                           1, 2,
                           2);

    igraph_vector_int_init_int(&terminals_k6_k1, 4,
                               0, 1, 2, 6);

    igraph_small(&g_k7_real, 7, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,
                 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
                 2, 3, 2, 4, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6,
                 5, 6,
                 -1);
    // 1 , 1
    igraph_vector_init_real(&weights_k7_real, 21,
                            2.5, 2.5, 2.5, 0.5, 0.5, 2.5,
                            2.5, 2.5, 2.5, 0.5, 2.5,
                            2.5, 2.5, 2.5, 0.5,
                            0.5, 2.5, 0.5,
                            2.5, 0.5,
                            0.5);

    igraph_vector_int_init_int(&terminals_k7_real, 4,
                               0, 1, 2, 3);

    /* K_7 complete graph with a specific edge ordering. */
    igraph_small(&g_k7_non_simple, 7, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 5, 0, 6,
                 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
                 2, 3, 2, 4, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6,
                 5, 6,
                 -1);
    // 1 , 1
    igraph_vector_init_int(&weights_k7_non_simple, 22,
                           2, 2, 2, 1, 1, 2, 2,
                           2, 2, 2, 1, 2,
                           2, 2, 2, 1,
                           1, 2, 1,
                           2, 1,
                           1);

    // A graph with 1 node.
    igraph_small(&g_k1, 1, IGRAPH_UNDIRECTED, -1);

    /*
        K_7 graph with a specific edge ordering and a self loop
        edge 0-1 is removed so it's not a complate graph
    */
    igraph_small(&g_k7_self_loop, 7, IGRAPH_UNDIRECTED,
                 0, 0, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,
                 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
                 2, 3, 2, 4, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6,
                 5, 6,
                 -1);
    // 1 , 1
    igraph_vector_init_int(&weights_k7_loop, 21,
                           2, 2, 2, 1, 1, 2,
                           2, 2, 2, 1, 2,
                           2, 2, 2, 1,
                           1, 2, 1,
                           2, 1,
                           1);

    igraph_full(&g_10_9, 10, IGRAPH_UNDIRECTED, 0);
    igraph_vector_int_init_int(&terminals_10_9, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8);

    igraph_vector_int_init_int(&terminals_1, 1, 0);
    igraph_vector_int_init_int(&terminals_2, 2, 1, 5);
    igraph_vector_int_init_int(&terminals_3, 3, 1, 3, 5);

    igraph_small(&g_new, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,
                 2, 7, 2, 8, 2, 9,
                 -1);

    igraph_vector_int_init_int(&terminals_8, 8, 0, 2, 3, 4, 5, 6, 7, 8);
    igraph_vector_int_init_int(&terminals_7_special, 7, 0, 3, 4, 5, 6, 7, 8);
    igraph_vector_int_init_int(&terminals_7, 7, 0, 1, 2, 3, 4, 5, 7);

    igraph_small(&g_new_1, 9, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 6,
                 2, 7, 2, 9,
                 3, 4,
                 4, 5,
                 7, 8,
                 -1
                );

    igraph_small(&g_new_2, 6, IGRAPH_UNDIRECTED,
                 0, 1,
                 1, 2, 1, 3,
                 3, 4,
                 2, 4,
                 4, 5,
                 -1
                );

    igraph_vector_int_init_int(&terminals_square, 4, 0, 3, 4, 5);

    igraph_small(&g_new_3, 9, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2,
                 0, 3, 1, 4, 2, 5,
                 3, 4, 4, 5,
                 3, 6, 4, 7, 5, 8,
                 6, 7, 7, 8,
                 -1);

    igraph_vector_init(&weights_square_2, 12);
    igraph_vector_fill(&weights_square_2, 100);

    igraph_vector_int_init_int_end(&terminals_square_2, -1, 0, 4, 2, 8, -1);

    igraph_vector_init_real(&weights_fp, 8, 0.218959, 0.383502, 0.5297, 0.417486, 0.526929, 0.910321, 0.328234, 0.247039);
    igraph_small(&g_new_4, 6, IGRAPH_UNDIRECTED,
                 5, 2,
                 5, 2,
                 2, 0,
                 4, 1,
                 5, 3,
                 2, 0,
                 2, 1,
                 3, 1,
                 -1);
    igraph_vector_int_init_int(&terminals_fpw, 3, 2, 0, 1);

    printf("Null graph:\n");
    check_graph(&g_null, &terminals_null, &weights_null);

    printf("\nK_7 graph:\n");
    check_graph(&g_k7, &terminals_k7, &weights_k7);

    printf("\nK_7 Non Complete graph-1:\n");
    check_graph(&g_k7_n, &terminals_k7, &weights_k7_n);

    printf("\nK_7 Non Complete graph-2:\n");
    check_graph(&g_k7_n1, &terminals_k7, &weights_k7_n1);

    igraph_real_t value;
    igraph_vector_int_t tree_edges;
    igraph_vector_int_init(&tree_edges, 0);

    CHECK_ERROR(igraph_steiner_dreyfus_wagner(&g_k6_k1, &terminals_k6_k1, &weights_k6_k1, &value, &tree_edges), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&tree_edges);

    igraph_real_t value_unit;
    igraph_vector_int_t tree_edges_unit;
    igraph_vector_int_init(&tree_edges_unit, 0);

    printf("\nK_7_Unit graph:\n");
    igraph_steiner_dreyfus_wagner(&g_k7, &terminals_k7, NULL, &value_unit, &tree_edges_unit);
    print_vector_int(&tree_edges_unit);
    printf("Total Steiner tree weight: %g\n", value_unit);

    igraph_vector_int_destroy(&tree_edges_unit);

    printf("\nK_7_Real graph:\n");
    check_graph(&g_k7_real, &terminals_k7_real, &weights_k7_real);

    printf("\nK_7 graph with parallel edges:\n");
    check_graph(&g_k7_non_simple, &terminals_k7, &weights_k7_non_simple);

    printf("\nA graph with single node:\n");
    check_graph(&g_k1, &terminals_null, &weights_null);

    printf("\nK_7 graph with self loop:\n");
    check_graph(&g_k7_self_loop, &terminals_k7, &weights_k7_loop);

    printf("\nA graph with n-1 terminals:\n");
    check_graph(&g_10_9, &terminals_10_9, NULL);

    printf("\nA graph with 1 terminal:\n");
    check_graph(&g_10_9, &terminals_1, NULL);

    printf("\nA graph with 2 terminals:\n");
    check_graph(&g_10_9, &terminals_2, NULL);

    printf("\na graph with 3 terminals:\n");
    check_graph(&g_10_9, &terminals_3, NULL);

    printf("\na graph:\n");
    check_graph(&g_new, &terminals_10_9, NULL);

    printf("\na graph with 8 terminals:\n");
    check_graph(&g_new, &terminals_8, NULL);

    printf("\na graph with 7 nonconsecutive terminals:\n");
    check_graph(&g_new, &terminals_7_special, NULL);

    printf("\nanother graph with 7 terminals:\n");
    check_graph(&g_new_1, &terminals_7, NULL);

    printf("\nA simple square graph with few more edges outside the square\n");
    check_graph(&g_new_2, &terminals_square, NULL);

    printf("\nA 3 vertices square graph with weights\n");
    check_graph(&g_new_3, &terminals_square_2, &weights_square_2);

    printf("\nA different graph than before with floating point weights\n");
    check_graph(&g_new_4, &terminals_fpw, &weights_fp);

    igraph_destroy(&g_null);
    igraph_destroy(&g_k7);
    igraph_destroy(&g_k6_k1);
    igraph_destroy(&g_k7_real);
    igraph_destroy(&g_k7_n1);
    igraph_destroy(&g_k7_n);
    igraph_destroy(&g_k7_non_simple);
    igraph_destroy(&g_k1);
    igraph_destroy(&g_k7_self_loop);
    igraph_destroy(&g_10_9);
    igraph_destroy(&g_new);
    igraph_destroy(&g_new_1);
    igraph_destroy(&g_new_2);
    igraph_destroy(&g_new_3);
    igraph_destroy(&g_new_4);

    igraph_vector_destroy(&weights_null);
    igraph_vector_destroy(&weights_k7);
    igraph_vector_destroy(&weights_k6_k1);
    igraph_vector_destroy(&weights_k7_n);
    igraph_vector_destroy(&weights_k7_n1);
    igraph_vector_destroy(&weights_k7_real);
    igraph_vector_destroy(&weights_k7_non_simple);
    igraph_vector_destroy(&weights_k7_loop);
    igraph_vector_destroy(&weights_square_2);
    igraph_vector_destroy(&weights_fp);

    igraph_vector_int_destroy(&terminals_k7);
    igraph_vector_int_destroy(&terminals_null);
    igraph_vector_int_destroy(&terminals_k6_k1);
    igraph_vector_int_destroy(&terminals_k7_real);
    igraph_vector_int_destroy(&terminals_10_9);
    igraph_vector_int_destroy(&terminals_1);
    igraph_vector_int_destroy(&terminals_2);
    igraph_vector_int_destroy(&terminals_3);
    igraph_vector_int_destroy(&terminals_8);
    igraph_vector_int_destroy(&terminals_7);
    igraph_vector_int_destroy(&terminals_7_special);
    igraph_vector_int_destroy(&terminals_square);
    igraph_vector_int_destroy(&terminals_square_2);
    igraph_vector_int_destroy(&terminals_fpw);


    VERIFY_FINALLY_STACK();

    return 0;
}
