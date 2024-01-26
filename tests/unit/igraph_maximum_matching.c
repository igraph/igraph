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

void test_trivial_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t edges;
    igraph_vector_int_t matching;

    igraph_vector_int_init(&edges,0);
    igraph_vector_int_init(&matching,0);

    //test null graph
    igraph_integer_t null_matching_size;
    igraph_create(&graph, &edges, 0, 0);
    igraph_maximum_matching(&graph, &null_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(null_matching_size == 0);

    //test singleton graph
    igraph_integer_t single_matching_size;
    igraph_create(&graph, &edges, 1, false);
    igraph_maximum_matching(&graph, &single_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(single_matching_size == 0);

    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_path_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);

    //test path graphs lengths 1-3
    igraph_integer_t path1_matching_size;
    igraph_ring(&graph, 2, IGRAPH_UNDIRECTED, 0, 1);
    igraph_maximum_matching(&graph, &path1_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(path1_matching_size == 1);

    igraph_integer_t path2_matching_size;
    igraph_ring(&graph, 3, IGRAPH_UNDIRECTED, 0, 1);
    igraph_maximum_matching(&graph, &path2_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(path2_matching_size == 2);

    igraph_integer_t path3_matching_size;
    igraph_ring(&graph, 4, IGRAPH_UNDIRECTED, 0, 1);
    igraph_maximum_matching(&graph, &path3_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(path3_matching_size == 3);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_bipartite_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);

    //test two bipartite graphs
    igraph_integer_t bi1_matching_size;
    igraph_small(&graph, 9, IGRAPH_UNDIRECTED, 1,2, 2,3, 2,4, 2,5, 3,6, 4,7, 5,8, 6,9, -1);
    igraph_maximum_matching(&graph, &bi1_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(bi1_matching_size == 4);

    igraph_integer_t bi2_matching_size;
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,4, 3,4, 3,5, 4,6, 5,6, -1);
    igraph_maximum_matching(&graph, &bi2_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(bi2_matching_size == 3);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_general_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);

    //test two non-bipartite graphs
    igraph_integer_t gen1_matching_size;
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,9, 3,4, 3,5, 4,6, 5,7, 5,8, 7,8, 8,9, 8,10, -1);
    igraph_maximum_matching(&graph, &gen1_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(gen1_matching_size == 5);

    igraph_integer_t gen2_matching_size;
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,4, 3,5, 4,5, 5,6, -1);
    igraph_maximum_matching(&graph, &gen2_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(gen2_matching_size == 3);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

void test_petersen_graph() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_vector_int_init(&matching,0);

    //test petersen graph
    igraph_integer_t petersen_matching_size;
    igraph_famous(&graph, "Petersen");
    igraph_maximum_matching(&graph, &petersen_matching_size, NULL, &matching, NULL);
    IGRAPH_ASSERT(petersen_matching_size == 5);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
}

int main() {
    test_trivial_graphs();
    test_path_graphs();
    test_bipartite_graphs();
    test_general_graphs();
    test_petersen_graph();
    return 0;
}
