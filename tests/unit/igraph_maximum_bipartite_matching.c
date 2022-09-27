/* IGraph library.
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

int test_weighted_graph_from_mit_notes(void) {
    /* Test graph from the following lecture notes:
     * http://math.mit.edu/~goemans/18433S07/matching-notes.pdf
     */
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_vector_int_t matching;
    igraph_vector_t weights;
    igraph_integer_t matching_size;
    igraph_real_t matching_weight;
    igraph_bool_t is_matching;
    igraph_real_t weight_array[] = { 2, 7, 2, 3,
                                     1, 3, 9, 3, 3,
                                     1, 3, 3, 1, 2,
                                     4, 1, 2,
                                     3
                                   };
    int i;

    igraph_small(&graph, 0, 0,
                 0, 6, 0, 7, 0, 8, 0, 9,
                 1, 5, 1, 6, 1, 7, 1, 8, 1, 9,
                 2, 5, 2, 6, 2, 7, 2, 8, 2, 9,
                 3, 5, 3, 7, 3, 9,
                 4, 7, -1);
    igraph_vector_bool_init(&types, 10);
    for (i = 0; i < 10; i++) {
        VECTOR(types)[i] = (i >= 5);
    }
    igraph_vector_int_init(&matching, 0);
    igraph_vector_init_array(&weights, weight_array,
                             sizeof(weight_array) / sizeof(weight_array[0]));

    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, &weights, 0);
    if (matching_size != 4) {
        printf("matching_size is %" IGRAPH_PRId ", expected: 4\n", matching_size);
        return 1;
    }
    if (matching_weight != 19) {
        printf("matching_weight is %" IGRAPH_PRId ", expected: 19\n", (igraph_integer_t) matching_weight);
        return 2;
    }
    igraph_is_maximal_matching(&graph, &types, &matching, &is_matching);
    if (!is_matching) {
        printf("not a matching: ");
        igraph_vector_int_print(&matching);
        return 3;
    }

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&matching);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    return 0;
}

int test_weighted_graph_generated(void) {
    /* Several randomly generated small test graphs */
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_vector_int_t matching;
    igraph_vector_t weights;
    igraph_integer_t matching_size;
    igraph_real_t matching_weight;
    igraph_real_t weight_array_1[] = { 8, 5, 9, 18, 20, 13 };
    igraph_real_t weight_array_2[] = { 20, 4, 20, 3, 13, 1 };
    int i;

    igraph_vector_bool_init(&types, 10);
    for (i = 0; i < 10; i++) {
        VECTOR(types)[i] = (i >= 5);
    }
    igraph_vector_int_init(&matching, 0);

    /* Case 1 */

    igraph_small(&graph, 0, 0, 0, 8, 2, 7, 3, 7, 3, 8, 4, 5, 4, 9, -1);
    igraph_vector_init_array(&weights, weight_array_1,
                             sizeof(weight_array_1) / sizeof(weight_array_1[0]));
    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, &weights, 0);
    if (matching_weight != 43) {
        printf("matching_weight is %" IGRAPH_PRId ", expected: 43\n", (igraph_integer_t)matching_weight);
        return 2;
    }
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    /* Case 2 */

    igraph_small(&graph, 0, 0, 0, 5, 0, 6, 1, 7, 2, 5, 3, 5, 3, 9, -1);
    igraph_vector_init_array(&weights, weight_array_2,
                             sizeof(weight_array_2) / sizeof(weight_array_2[0]));
    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, &weights, 0);
    if (matching_weight != 41) {
        printf("matching_weight is %" IGRAPH_PRId ", expected: 41\n", (igraph_integer_t)matching_weight);
        return 2;
    }
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    igraph_vector_int_destroy(&matching);
    igraph_vector_bool_destroy(&types);

    return 0;
}

int main(void) {
    igraph_t g;
    igraph_vector_bool_t types;
    igraph_vector_t weights;

    igraph_integer_t matching_size;
    igraph_real_t weighted_size;

    igraph_vector_int_t matching;

    igraph_small(&g, 4, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3,
                 -1);

    igraph_vector_bool_init(&types, 4);
    VECTOR(types)[0] = 0;
    VECTOR(types)[1] = 1;
    VECTOR(types)[2] = 0;
    VECTOR(types)[3] = 1;

    igraph_vector_int_init(&matching, 0);

    igraph_vector_init(&weights, igraph_vcount(&g));
    igraph_vector_fill(&weights, 1.0);

    // Test incorrect types
    CHECK_ERROR(igraph_maximum_bipartite_matching(&g, &types, &matching_size, NULL, &matching, NULL, 0), IGRAPH_EINVAL);

    // Test incorrect types for weighted graph
    VECTOR(types)[2] = 0;
    CHECK_ERROR(igraph_maximum_bipartite_matching(&g, &types, &matching_size, &weighted_size, &matching, &weights, 0), IGRAPH_EINVAL);

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&matching);

    igraph_vector_bool_destroy(&types);
    igraph_destroy(&g);

    if (test_weighted_graph_generated()) {
        return 1;
    }

    if (test_weighted_graph_from_mit_notes()) {
        return 2;
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
