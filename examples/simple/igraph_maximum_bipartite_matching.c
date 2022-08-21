/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2012  Tamas Nepusz <ntamas@gmail.com>

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

int test_graph_from_leda_tutorial() {
    /* Test graph from the LEDA tutorial:
     * http://www.leda-tutorial.org/en/unofficial/ch05s03s05.html
     */
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_vector_int_t matching;
    igraph_integer_t matching_size;
    igraph_real_t matching_weight;
    igraph_bool_t is_matching;
    int i;

    igraph_small(&graph, 0, 0,
                 0, 8, 0, 12, 0, 14,
                 1, 9, 1, 10, 1, 13,
                 2, 8, 2, 9,
                 3, 10, 3, 11, 3, 13,
                 4, 9, 4, 14,
                 5, 14,
                 6, 9, 6, 14,
                 7, 8, 7, 12, 7, 14
                 , -1);
    igraph_vector_bool_init(&types, 15);
    for (i = 0; i < 15; i++) {
        VECTOR(types)[i] = (i >= 8);
    }
    igraph_vector_int_init(&matching, 0);

    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, 0, 0);
    if (matching_size != 6) {
        printf("matching_size is %" IGRAPH_PRId ", expected: 6\n", matching_size);
        return 1;
    }
    if (matching_weight != 6) {
        printf("matching_weight is %" IGRAPH_PRId ", expected: 6\n", (igraph_integer_t) matching_weight);
        return 2;
    }
    igraph_is_maximal_matching(&graph, &types, &matching, &is_matching);
    if (!is_matching) {
        printf("not a matching: ");
        igraph_vector_int_print(&matching);
        return 3;
    }

    igraph_vector_int_destroy(&matching);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    return 0;
}

int test_weighted_graph_from_mit_notes() {
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

int test_weighted_graph_generated() {
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

int main() {
    if (test_graph_from_leda_tutorial()) {
        return 1;
    }

    if (test_weighted_graph_from_mit_notes()) {
        return 2;
    }

    if (test_weighted_graph_generated()) {
        return 3;
    }

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        printf("Finally stack still has %d elements.\n", IGRAPH_FINALLY_STACK_SIZE());
        return 5;
    }

    return 0;
}
