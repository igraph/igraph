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

int test_trivial_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t edges;
    igraph_vector_int_t matching;
    igraph_integer_t matching_size;

    igraph_vector_int_init(&edges,0);
    igraph_vector_int_init(&matching,0);
    //test null graph
    igraph_create(&graph, &edges, 0, 0);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 0) {
        printf("Null graph had %" IGRAPH_PRId " matching size.\n", matching_size);
        result = 1;
    }
    //test singleton graph
    igraph_create(&graph, &edges, 1, false);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 0) {
        printf("Singleton graph had %" IGRAPH_PRId " matching size.\n", matching_size);
        result = 2;
    }
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
    return result;
}

int test_path_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_integer_t matching_size;
    igraph_vector_int_init(&matching,0);
    //test path graphs lengths 1-3
    igraph_ring(&graph, 2, IGRAPH_UNDIRECTED, 0, 1);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 1) {
        printf("Length 1 path graph had %" IGRAPH_PRId " matching size.\n", matching_size);
        result = 1;
    }
    igraph_ring(&graph, 3, IGRAPH_UNDIRECTED, 0, 1);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 2) {
        printf("Length 2 path graph had %" IGRAPH_PRId " matching size.\n", matching_size);
        result = 2;
    }
    igraph_ring(&graph, 4, IGRAPH_UNDIRECTED, 0, 1);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 3) {
        printf("Length 3 path graph had %" IGRAPH_PRId " matching size.\n", matching_size);
        result = 3;
    }
    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
    return result;
}

int test_bipartite_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_integer_t matching_size;
    igraph_vector_int_init(&matching,0);
    //test two bipartite graphs
    igraph_small(&graph, 9, IGRAPH_UNDIRECTED, 1,2, 2,3, 2,4, 2,5, 3,6, 4,7, 5,8, 6,9, -1);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 4) {
        printf("Bipartite graph 1 had %" IGRAPH_PRId " matching size, expected 4.\n", matching_size);
        result = 1;
    }
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,4, 3,4, 3,5, 4,6, 5,6, -1);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 3) {
        printf("Bipartite graph 2 had %" IGRAPH_PRId " matching size, expected 3.\n", matching_size);
        result = 2;
    }
    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
    return result;
}

int test_general_graphs() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_integer_t matching_size;
    igraph_vector_int_init(&matching,0);
    //test two non-bipartite graphs
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,9, 3,4, 3,5, 4,6, 5,7, 5,8, 7,8, 8,9, 8,10, -1);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 5) {
        printf("General graph 1 had %" IGRAPH_PRId " matching size, expected 5.\n", matching_size);
        result = 1;
    }
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED, 1,2, 1,3, 2,4, 3,5, 4,5, 5,6, -1);
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 3) {
        printf("General graph 2 had %" IGRAPH_PRId " matching size, expected 3.\n", matching_size);
        result = 2;
    }
    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
    return result;
}

int test_petersen_graph() {
    int result = 0;
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_integer_t matching_size;
    igraph_vector_int_init(&matching,0);
    //test petersen graph
    igraph_famous(&graph, "Petersen");
    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL);
    if (matching_size != 5) {
        printf("Petersen graph had %" IGRAPH_PRId " matching size, expected 5.\n", matching_size);
        result = 1;
    }
    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);
    return result;
}

int main() {
    int result = 0;
    if (test_trivial_graphs()) {
        result = 1;
    }
    if (test_path_graphs()) {
        result = 1;
    }
    if (test_bipartite_graphs()) {
        result = 1;
    }
    if (test_general_graphs()) {
        result = 1;
    }
    if (test_petersen_graph()) {
        result = 1;
    }
    return result;
}
