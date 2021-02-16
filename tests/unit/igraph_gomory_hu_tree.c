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

#include "test_utilities.inc"

int validate_tree(const igraph_t *graph, const igraph_t *tree,
                  const igraph_vector_t *flow, const igraph_vector_t *capacity) {
    igraph_integer_t n = igraph_vcount(graph);
    igraph_vector_t edges;
    igraph_real_t min_weight, flow_value;
    long int i, j, k, m;

    if (igraph_vcount(tree) != n) {
        printf("Gomory-Hu tree should have %ld vertices\n", (long int)n);
        return IGRAPH_EINVAL;
    }

    if (igraph_ecount(tree) != n - 1) {
        printf("Gomory-Hu tree should have %ld edges\n", (long int)n - 1);
        return IGRAPH_EINVAL;
    }

    if (igraph_is_directed(tree)) {
        printf("Gomory-Hu tree should be undirected\n");
        return IGRAPH_EINVAL;
    }

    if (n < 2) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            IGRAPH_CHECK(igraph_get_shortest_path(tree, 0, &edges, i, j, IGRAPH_ALL));
            m = igraph_vector_size(&edges);
            if (m == 0) {
                continue;
            }

            min_weight = VECTOR(*flow)[(long int)VECTOR(edges)[0]];
            for (k = 1; k < m; k++) {
                if (VECTOR(*flow)[(long int)VECTOR(edges)[k]] < min_weight) {
                    min_weight = VECTOR(*flow)[(long int)VECTOR(edges)[k]];
                }
            }

            IGRAPH_CHECK(igraph_maxflow(graph, &flow_value, 0, 0, 0, 0, i, j, capacity, 0));
            if (flow_value != min_weight) {
                printf("Min weight of path %ld --> %ld in Gomory-Hu tree is %.4f, "
                       "expected %.4f from flow calculation\n", i, j, min_weight, flow_value);
                return IGRAPH_EINVAL;
            }
        }
    }

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

int main() {

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

    /* simple directed graph - should throw an error */
    igraph_small(&g, 6, 1, 0, 1, 0, 2, 1, 2, 1, 3, 1, 4, 2, 4, 3, 4, 3, 5, 4, 5, -1);
    igraph_set_error_handler(igraph_error_handler_ignore);
    IGRAPH_ASSERT(igraph_gomory_hu_tree(&g, &tree, &flow, &capacity) == IGRAPH_EINVAL);
    igraph_set_error_handler(igraph_error_handler_abort);
    igraph_destroy(&g);

    /* destroy flow and capacity vectors */
    igraph_vector_destroy(&flow);
    igraph_vector_destroy(&capacity);

    VERIFY_FINALLY_STACK();

    return 0;
}
