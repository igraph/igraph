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
    igraph_vector_long_t matching;
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
    igraph_vector_long_init(&matching, 0);

    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, 0, 0);
    if (matching_size != 6) {
        printf("matching_size is %ld, expected: 6\n", (long)matching_size);
        return 1;
    }
    if (matching_weight != 6) {
        printf("matching_weight is %ld, expected: 6\n", (long)matching_weight);
        return 2;
    }
    igraph_is_maximal_matching(&graph, &types, &matching, &is_matching);
    if (!is_matching) {
        printf("not a matching: ");
        igraph_vector_long_print(&matching);
        return 3;
    }

    igraph_vector_long_destroy(&matching);
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
    igraph_vector_long_t matching;
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
    igraph_vector_long_init(&matching, 0);
    igraph_vector_init_copy(&weights, weight_array,
                            sizeof(weight_array) / sizeof(weight_array[0]));

    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, &weights, 0);
    if (matching_size != 4) {
        printf("matching_size is %ld, expected: 4\n", (long)matching_size);
        return 1;
    }
    if (matching_weight != 19) {
        printf("matching_weight is %ld, expected: 19\n", (long)matching_weight);
        return 2;
    }
    igraph_is_maximal_matching(&graph, &types, &matching, &is_matching);
    if (!is_matching) {
        printf("not a matching: ");
        igraph_vector_long_print(&matching);
        return 3;
    }

    igraph_vector_destroy(&weights);
    igraph_vector_long_destroy(&matching);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    return 0;
}

int test_weighted_graph_generated() {
    /* Several randomly generated small test graphs */
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_vector_long_t matching;
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
    igraph_vector_long_init(&matching, 0);

    /* Case 1 */

    igraph_small(&graph, 0, 0, 0, 8, 2, 7, 3, 7, 3, 8, 4, 5, 4, 9, -1);
    igraph_vector_init_copy(&weights, weight_array_1,
                            sizeof(weight_array_1) / sizeof(weight_array_1[0]));
    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, &weights, 0);
    if (matching_weight != 43) {
        printf("matching_weight is %ld, expected: 43\n", (long)matching_weight);
        return 2;
    }
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    /* Case 2 */

    igraph_small(&graph, 0, 0, 0, 5, 0, 6, 1, 7, 2, 5, 3, 5, 3, 9, -1);
    igraph_vector_init_copy(&weights, weight_array_2,
                            sizeof(weight_array_2) / sizeof(weight_array_2[0]));
    igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
                                      &matching_weight, &matching, &weights, 0);
    if (matching_weight != 41) {
        printf("matching_weight is %ld, expected: 41\n", (long)matching_weight);
        return 2;
    }
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    igraph_vector_long_destroy(&matching);
    igraph_vector_bool_destroy(&types);

    return 0;
}

int test_weighted_graph_from_file(const char* fname, int type1_count, long exp_weight) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_vector_long_t matching;
    igraph_vector_t weights;
    igraph_real_t matching_weight;
    FILE* f;
    int i, n;

    f = fopen(fname, "r");
    if (!f) {
        fprintf(stderr, "No such file: %s\n", fname);
        return 1;
    }
    igraph_read_graph_ncol(&graph, f, 0, 1, IGRAPH_ADD_WEIGHTS_YES, 0);
    fclose(f);

    n = igraph_vcount(&graph);
    igraph_vector_bool_init(&types, n);
    for (i = 0; i < n; i++) {
        VECTOR(types)[i] = (i >= type1_count);
    }

    igraph_vector_long_init(&matching, 0);

    igraph_vector_init(&weights, 0);
    EANV(&graph, "weight", &weights);
    igraph_maximum_bipartite_matching(&graph, &types, 0, &matching_weight,
                                      &matching, &weights, 0);
    igraph_vector_destroy(&weights);

    igraph_vector_long_print(&matching);
    if (matching_weight != exp_weight) {
        printf("matching_weight is %ld, expected: %ld\n", (long)matching_weight,
               (long)exp_weight);
        return 2;
    }

    igraph_vector_destroy(&weights);
    igraph_vector_long_destroy(&matching);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    return 0;
}

// This test addresses issue #1110, where an incorrect
// types vector (i.e. that doesn't correspond to a bipartite
// labelling of the graph) would cause a possible infinite loop.
int test_incorrect_types() {
    igraph_t g;
    igraph_vector_bool_t types;
    igraph_vector_t weights;

    igraph_integer_t matching_size;
    igraph_real_t weighted_size;

    igraph_vector_long_t matching;

    igraph_error_type_t err;


    igraph_small(&g, 4, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3,
                 -1);

    igraph_vector_bool_init(&types, 4);
    VECTOR(types)[0] = 0;
    VECTOR(types)[1] = 1;
    VECTOR(types)[2] = 0;
    VECTOR(types)[3] = 1;

    igraph_vector_long_init(&matching, 0);

    igraph_vector_init(&weights, igraph_vcount(&g));
    igraph_vector_fill(&weights, 1.0);

    igraph_set_error_handler(&igraph_error_handler_ignore);

    // Test incorrect types
    err = igraph_maximum_bipartite_matching(&g, &types, &matching_size, NULL, &matching, NULL, 0);
    if (err != IGRAPH_EINVAL) {
        return 3;
    }

    // Test correct types
    VECTOR(types)[2] = 1;
    err = igraph_maximum_bipartite_matching(&g, &types, &matching_size, NULL, &matching, NULL, 0);
    if (err == IGRAPH_EINVAL) {
        return 4;
    }

    // Test incorrect types for weighted graph
    VECTOR(types)[2] = 0;
    err = igraph_maximum_bipartite_matching(&g, &types, &matching_size, &weighted_size, &matching, &weights, 0);
    if (err != IGRAPH_EINVAL) {
        return 5;
    }

    // Test correct types for weighted graph
    VECTOR(types)[2] = 1;
    err = igraph_maximum_bipartite_matching(&g, &types, &matching_size, &weighted_size, &matching, &weights, 0);
    if (err == IGRAPH_EINVAL) {
        return 6;
    }

    igraph_vector_destroy(&weights);
    igraph_vector_long_destroy(&matching);

    igraph_vector_bool_destroy(&types);
    igraph_destroy(&g);

    return 0;
}

int main() {
    igraph_set_attribute_table(&igraph_cattribute_table);

    if (test_graph_from_leda_tutorial()) {
        return 1;
    }

    if (test_weighted_graph_from_mit_notes()) {
        return 2;
    }

    if (test_weighted_graph_generated()) {
        return 3;
    }

    if (test_incorrect_types()) {
        return 4;
    }

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        printf("Finally stack still has %d elements.\n", IGRAPH_FINALLY_STACK_SIZE());
        return 5;
    }

    return 0;
}
