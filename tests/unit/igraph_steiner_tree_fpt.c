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

int main(void) {
    igraph_t g_null, g_k7, g_k6_k1;
    igraph_vector_int_t terminals_null, terminals_k7, terminals_k6_k1;
    igraph_vector_t weights_null, weights_k7, weights_k6_k1;

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

    igraph_vector_init_int(&weights_k7, 21,
                           2, 2, 2, 1, 1, 2,
                           2, 2, 2, 1, 2,
                           2, 2, 2, 1,
                           1, 2, 1,
                           2, 1,
                           1);

    igraph_vector_int_init_int(&terminals_k7, 4,
                               0, 1, 2, 3);

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

    printf("No vertices, not directed:\n");
    igraph_real_t val1, val2,val3;
    igraph_vector_int_t res_tree, res_tree_1,res_tree_2;

    igraph_vector_int_init(&res_tree, 0);
    igraph_vector_int_init(&res_tree_1, 0);
    igraph_vector_int_init(&res_tree_2, 0);

    IGRAPH_ASSERT(igraph_steiner_dreyfus_wagner(&g_null, &terminals_null, &weights_null, &val1, &res_tree) == IGRAPH_SUCCESS);
    printf("%.2f\n", val1);
    IGRAPH_ASSERT(val1 == 0);
    printf("Un-Directed graph with loops and multi-edges, select none:\n");
    IGRAPH_ASSERT(igraph_steiner_dreyfus_wagner(&g_k7, &terminals_k7, &weights_k7, &val2, &res_tree_1) == IGRAPH_SUCCESS);
    printf("%.2f\n", val2);
    IGRAPH_ASSERT(val2 == 5);

    print_vector_int(&res_tree_1);

    CHECK_ERROR(igraph_steiner_dreyfus_wagner(&g_k6_k1, &terminals_k6_k1, &weights_k6_k1, &val3, &res_tree_2), IGRAPH_EINVAL);

    igraph_destroy(&g_null);
    igraph_destroy(&g_k7);
    igraph_destroy(&g_k6_k1);

    igraph_vector_destroy(&weights_null);
    igraph_vector_destroy(&weights_k7);
    igraph_vector_destroy(&weights_k6_k1);

    igraph_vector_int_destroy(&terminals_k7);
    igraph_vector_int_destroy(&terminals_null);
    igraph_vector_int_destroy(&terminals_k6_k1);

    igraph_vector_int_destroy(&res_tree);
    igraph_vector_int_destroy(&res_tree_1);
    igraph_vector_int_destroy(&res_tree_2);

    VERIFY_FINALLY_STACK();

    return 0;
}
