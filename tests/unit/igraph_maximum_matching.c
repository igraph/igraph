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

void test_trivial_graphs(void) {
    igraph_t graph;
    igraph_vector_int_t edges;
    igraph_vector_int_t matching;

    igraph_vector_int_init(&edges,0);
    igraph_vector_int_init(&matching,0);

    //test null graph
    igraph_integer_t null_matching_size;
    igraph_create(&graph, &edges, 0, 0);
    igraph_maximum_matching(&graph, &null_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(null_matching_size == 0);

    //test singleton graph
    igraph_integer_t single_matching_size;
    igraph_create(&graph, &edges, 1, false);
    igraph_maximum_matching(&graph, &single_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(single_matching_size == 0);

    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_path_graphs(void) {
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);
    igraph_bool_t is_matching;

    //test path graphs lengths 1-3
    igraph_integer_t path1_matching_size;
    igraph_ring(&graph, 2, IGRAPH_UNDIRECTED, 0, false);
    igraph_maximum_matching(&graph, &path1_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(path1_matching_size == 1);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t path2_matching_size;
    igraph_ring(&graph, 3, IGRAPH_UNDIRECTED, 0, false);
    igraph_maximum_matching(&graph, &path2_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(path2_matching_size == 1);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t path3_matching_size;
    igraph_ring(&graph, 4, IGRAPH_UNDIRECTED, 0, false);
    igraph_maximum_matching(&graph, &path3_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(path3_matching_size == 2);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_bipartite_graphs(void) {
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);
    igraph_bool_t is_matching;

    //test two bipartite graphs
    igraph_integer_t bi1_matching_size;
    igraph_small(&graph, 9, IGRAPH_UNDIRECTED, 1,2, 2,3, 2,4, 2,5, 3,6, 4,7, 5,8, 6,9, -1);
    igraph_maximum_matching(&graph, &bi1_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(bi1_matching_size == 4);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t bi2_matching_size;
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,4, 3,4, 3,5, 4,6, 5,6, -1);
    igraph_maximum_matching(&graph, &bi2_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(bi2_matching_size == 3);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_general_graphs(void) {
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);
    igraph_bool_t is_matching;

    //test two non-bipartite graphs
    igraph_integer_t gen1_matching_size;
    igraph_small(&graph, 11, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,9, 3,4, 3,5, 4,6, 5,7, 5,8, 7,8, 8,9, 8,10, -1);
    igraph_maximum_matching(&graph, &gen1_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(gen1_matching_size == 5);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t gen2_matching_size;
    igraph_small(&graph, 7, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,4, 3,5, 4,5, 5,6, -1);
    igraph_maximum_matching(&graph, &gen2_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(gen2_matching_size == 3);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    // next three are designed to test blossom contraction
    igraph_integer_t gen3_matching_size;
    igraph_small(&graph, 13, IGRAPH_UNDIRECTED,
                 3,1, 1,2, 2,4, 2,5, 4,6, 5,7, 6,7,
                 5,8, 8,9, 9,10, 10,11, 11,12, -1);
    igraph_maximum_matching(&graph, &gen3_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(gen3_matching_size == 6);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t gen4_matching_size;
    igraph_small(&graph, 13, IGRAPH_UNDIRECTED,
                 3,1, 1,2, 2,4, 2,5, 4,6, 5,7, 6,7,
                 7,8, 8,9, 9,10, 10,11, 11,12, -1);
    igraph_maximum_matching(&graph, &gen4_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(gen4_matching_size == 6);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t gen5_matching_size;
    igraph_small(&graph, 21, IGRAPH_UNDIRECTED,
                 3,1, 1,2, 2,4, 2,5, 4,6, 5,7, 6,7,
                 5,8, 8,9, 9,10, 9,11, 10,12, 11,13, 12,13,
                 11,14, 14,15, 15,16, 16,17, 17,18, 18,19, 19,20, -1);
    igraph_maximum_matching(&graph, &gen5_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(gen5_matching_size == 10);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t gen6_matching_size;
    igraph_small(&graph, 19, IGRAPH_UNDIRECTED,
                 1,2, 2,3, 2,9, 3,4, 4,5, 4,6, 5,7, 6,8, 7,8,
                 8,12, 9,10, 10,11, 11,12, 1,13, 7,14, 14,15,
                 15,16, 16,17, 17,18, -1);
    igraph_maximum_matching(&graph, &gen6_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(gen6_matching_size == 9);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_integer_t gen7_matching_size;
    igraph_small(&graph, 34, IGRAPH_UNDIRECTED,
                 1,2, 2,3, 2,4, 3,5, 4,6, 5,6, 4,7, 7,8, 8,9,
                 8,10, 9,11, 10,12, 11,12, 11,13, 13,14, 14,15,
                 14,16, 15,17, 16,18, 17,18, 15,19, 19,20, 20,21,
                 20,22, 21,23, 22,24, 23,24, 24,25, 25,26, 27,28,
                 28,29, 29,30, 29,31, 30, 32, 31,33, 32,33, 33,1, -1);
    igraph_maximum_matching(&graph, &gen7_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(gen7_matching_size == 16);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_petersen_graph(void) {
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);
    igraph_bool_t is_matching;

    //test petersen graph
    igraph_integer_t petersen_matching_size;
    igraph_famous(&graph, "Petersen");
    igraph_maximum_matching(&graph, &petersen_matching_size, NULL, &matching, NULL, 0);
    IGRAPH_ASSERT(petersen_matching_size == 5);
    igraph_is_maximal_matching(&graph, NULL, &matching, &is_matching);
    IGRAPH_ASSERT(is_matching);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

int main() {
    test_trivial_graphs();
    test_path_graphs();
    test_bipartite_graphs();
    test_general_graphs();
    test_petersen_graph();
    VERIFY_FINALLY_STACK();
    return 0;
}
