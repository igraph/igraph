/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2013  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

#include <igraph.h>

#include "test_utilities.h"

igraph_error_t validate_tree(const igraph_t *graph, const igraph_t *tree,
                  const igraph_vector_t *flow, const igraph_vector_t *capacity) {
    igraph_integer_t n = igraph_vcount(graph);
    igraph_integer_t no_of_clusters, min_weight_edge_index;
    igraph_vector_int_t edges;
    igraph_vector_int_t membership;
    igraph_real_t min_weight, flow_value;
    igraph_t copy;
    igraph_integer_t i, j, k, m;

    if (igraph_vcount(tree) != n) {
        printf("Gomory-Hu tree should have %" IGRAPH_PRId " vertices\n", n);
        return IGRAPH_EINVAL;
    }

    if (igraph_ecount(tree) != n - 1) {
        printf("Gomory-Hu tree should have %" IGRAPH_PRId " edges\n", n - 1);
        return IGRAPH_EINVAL;
    }

    if (igraph_is_directed(tree)) {
        printf("Gomory-Hu tree should be undirected\n");
        return IGRAPH_EINVAL;
    }

    if (n < 2) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            IGRAPH_CHECK(igraph_get_shortest_path(tree, 0, &edges, i, j, IGRAPH_ALL));
            m = igraph_vector_int_size(&edges);
            if (m == 0) {
                continue;
            }

            /* first, check whether the minimum weight along the shortest path
             * from i to j is the same as the maximum flow between i and j in
             * the original graph */
            min_weight = VECTOR(*flow)[VECTOR(edges)[0]];
            min_weight_edge_index = VECTOR(edges)[0];
            for (k = 1; k < m; k++) {
                if (VECTOR(*flow)[VECTOR(edges)[k]] < min_weight) {
                    min_weight = VECTOR(*flow)[VECTOR(edges)[k]];
                    min_weight_edge_index = VECTOR(edges)[k];
                }
            }

            IGRAPH_CHECK(igraph_maxflow(graph, &flow_value, 0, 0, 0, 0, i, j, capacity, 0));
            if (flow_value != min_weight) {
                printf("Min weight of path %" IGRAPH_PRId " -- %" IGRAPH_PRId " in Gomory-Hu tree is %.4f, "
                       "expected %.4f from flow calculation\n", i, j, min_weight, flow_value);
                return IGRAPH_EINVAL;
            }

            /* next, check whether removing an edge s-t from the Gomory-Hu tree would
             * partition it exactly the same way as a minimum cut between s and t in
             * the original graph */
            IGRAPH_CHECK(igraph_copy(&copy, tree));
            IGRAPH_FINALLY(igraph_destroy, &copy);

            IGRAPH_CHECK(igraph_delete_edges(&copy, igraph_ess_1(min_weight_edge_index)));
            IGRAPH_CHECK(igraph_connected_components(&copy, &membership, 0, &no_of_clusters, IGRAPH_WEAK));

            if (no_of_clusters != 2) {
                printf(
                    "Removing edge %" IGRAPH_PRId " -- %" IGRAPH_PRId
                    " (index %" IGRAPH_PRId ") from the Gomory-Hu tree cuts it "
                    "in %" IGRAPH_PRId " clusters, expected 2\n",
                    IGRAPH_FROM(tree, min_weight_edge_index),
                    IGRAPH_TO(tree, min_weight_edge_index),
                    min_weight_edge_index,
                    no_of_clusters
                );
                return IGRAPH_EINVAL;
            }

            /* finally, check the total capacity of the edges that go between the
             * partitions in the original graph; it should be the same as the
             * weight of the edge in the Gomory-Hu tree that corresponds to the
             * minimum weight along the path we found above */
            m = igraph_ecount(graph);
            flow_value = 0.0;
            for (j = 0; j < m; j++) {
                if (VECTOR(membership)[IGRAPH_FROM(graph, j)] != VECTOR(membership)[IGRAPH_TO(graph, j)]) {
                    flow_value += capacity ? VECTOR(*capacity)[j] : 1;
                }
            }

            if (flow_value != VECTOR(*flow)[min_weight_edge_index]) {
                printf(
                    "Edge %" IGRAPH_PRId " -- %" IGRAPH_PRId
                    " (index %" IGRAPH_PRId ") in the Gomory-Hu tree has weight = %.2f, but "
                    "the corresponding flow in the original graph has value = %.2f\n",
                    IGRAPH_FROM(tree, min_weight_edge_index),
                    IGRAPH_TO(tree, min_weight_edge_index),
                    min_weight_edge_index,
                    VECTOR(*flow)[min_weight_edge_index], flow_value
                );
                printf("Edge list of original graph:\n");
                print_graph(graph);
                if (capacity) {
                    printf("Capacities of original graph: ");
                    print_vector(capacity);
                } else {
                    printf("All edges have capacity = 1\n");
                }
                printf("Edge list of Gomory-Hu tree:\n");
                print_graph(tree);
                printf("Weights of the Gomory-Hu tree: ");
                print_vector(flow);
                printf("Partition of original graph corresponding to this cut:\n");
                print_vector_int(&membership);
                return IGRAPH_EINVAL;
            }

            igraph_destroy(&copy);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

int main(void) {

    igraph_t g;
    igraph_t tree;
    igraph_vector_t flow;
    igraph_vector_t capacity;

    /* initialize flow and capacity vectors */
    igraph_vector_init(&capacity, 0);
    igraph_vector_init(&flow, 0);

    /* empty undirected graph */
    igraph_empty(&g, 0, 0);
    IGRAPH_ASSERT(igraph_gomory_hu_tree(&g, &tree, &flow, &capacity) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&tree) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&flow) == 0);
    igraph_destroy(&tree);
    igraph_destroy(&g);

    /* simple undirected graph */
    igraph_small(&g, 6, 0, 0, 1, 0, 2, 1, 2, 1, 3, 1, 4, 2, 4, 3, 4, 3, 5, 4, 5, -1);
    igraph_vector_resize(&capacity, 9);
    VECTOR(capacity)[0] = 1;
    VECTOR(capacity)[1] = 7;
    VECTOR(capacity)[2] = 1;
    VECTOR(capacity)[3] = 3;
    VECTOR(capacity)[4] = 2;
    VECTOR(capacity)[5] = 4;
    VECTOR(capacity)[6] = 1;
    VECTOR(capacity)[7] = 6;
    VECTOR(capacity)[8] = 2;
    IGRAPH_ASSERT(igraph_gomory_hu_tree(&g, &tree, &flow, &capacity) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(validate_tree(&g, &tree, &flow, &capacity) == IGRAPH_SUCCESS);
    igraph_destroy(&tree);

    /* Make sure we don't blow up without an outgoing flow vector */
    IGRAPH_ASSERT(igraph_gomory_hu_tree(&g, &tree, 0, &capacity) == IGRAPH_SUCCESS);
    igraph_destroy(&tree);
    igraph_destroy(&g);

    /* example from Github issue #1810 */
    igraph_full(&g, 4, /* directed = */ 0, /* loops = */ 0);
    IGRAPH_ASSERT(igraph_gomory_hu_tree(&g, &tree, &flow, 0) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(validate_tree(&g, &tree, &flow, 0) == IGRAPH_SUCCESS);
    igraph_destroy(&tree);
    igraph_destroy(&g);

    /* simple directed graph - should throw an error */
    igraph_small(&g, 6, 1, 0, 1, 0, 2, 1, 2, 1, 3, 1, 4, 2, 4, 3, 4, 3, 5, 4, 5, -1);
    igraph_set_error_handler(igraph_error_handler_ignore);
    VERIFY_FINALLY_STACK();
    IGRAPH_ASSERT(igraph_gomory_hu_tree(&g, &tree, &flow, &capacity) == IGRAPH_EINVAL);
    igraph_set_error_handler(igraph_error_handler_abort);
    igraph_destroy(&g);

    /* destroy flow and capacity vectors */
    igraph_vector_destroy(&flow);
    igraph_vector_destroy(&capacity);

    VERIFY_FINALLY_STACK();

    return 0;
}
