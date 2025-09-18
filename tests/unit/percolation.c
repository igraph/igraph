/*
    igraph library.
    Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

igraph_error_t percolate_bond(igraph_t *graph, igraph_vector_int_t *edge_indices, igraph_bool_t printing) {
    igraph_vector_int_t giant_size, vertex_count;
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_int_t number_percolated;

    if (edge_indices != NULL) {
        number_percolated = igraph_vector_int_size(edge_indices);
    } else {
        number_percolated = ecount;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&giant_size, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_count, 0);

    IGRAPH_CHECK(igraph_bond_percolation(graph, &giant_size, &vertex_count, edge_indices));

    if (printing) {
        print_vector_int(&giant_size);
        print_vector_int(&vertex_count);
    }
    if (number_percolated == ecount) {
        // It's only guaranteed to have the same component ecount if all edges are included.
        igraph_vector_int_t component_sizes;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&component_sizes, 0);

        IGRAPH_CHECK(igraph_connected_components(graph, NULL, &component_sizes, NULL, IGRAPH_WEAK));

        IGRAPH_ASSERT(igraph_vector_int_size(&giant_size) == ecount);
        IGRAPH_ASSERT(igraph_vector_int_size(&vertex_count) == ecount);
        if (ecount > 1) {
            IGRAPH_ASSERT(igraph_vector_int_max(&giant_size) == igraph_vector_int_max(&component_sizes));
        }
        igraph_vector_int_destroy(&component_sizes);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_int_t prev = 0;
    for (igraph_int_t i = 0; i < number_percolated; i++) {
        IGRAPH_ASSERT(VECTOR(giant_size)[i] > 0);      // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(giant_size)[i] >= prev);   // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(giant_size)[i] <= i + 2); // Largest component cannot be bigger than a tree with the same number of edges.
        prev = VECTOR(giant_size)[i];
    }
    prev = 0;
    for (igraph_int_t i = 0; i < number_percolated; i++) {
        IGRAPH_ASSERT(VECTOR(vertex_count)[i] > 0);      // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(vertex_count)[i] >= prev);   // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(vertex_count)[i] <= 2*i + 2); // Largest component cannot be bigger than a tree with the same number of edges.
        prev = VECTOR(vertex_count)[i];
    }

    igraph_vector_int_destroy(&giant_size);
    igraph_vector_int_destroy(&vertex_count);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


void test_bond(void) {
    // Test with normal graph and provided edge list
    igraph_t k_3, c_4, karate, random, null_graph, singleton;

    igraph_empty(&null_graph, 0, false);
    igraph_empty(&singleton, 1, false);
    igraph_full(&k_3, 3, false, false);
    igraph_small(&c_4, 4, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_famous(&karate, "Zachary");
    igraph_erdos_renyi_game_gnp(&random, 100, 0.01, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    printf("# Bond percolation test suite\n");
    printf("Null graph, no provided edge order.\n");
    percolate_bond(&null_graph, NULL, true);

    printf("Singleton graph, no provided edge order.\n");
    percolate_bond(&singleton, NULL, true);

    printf("K_3 graph percolation curve, no provided edge sequence.\n");
    percolate_bond(&k_3, NULL, true); // sequence is random, but since it should be consistent it doesn't matter

    igraph_vector_int_t edge_ids;
    igraph_vector_int_init_int(&edge_ids, 4, 0, 2, 1, 3);

    printf("C_4 graph with edge sequence ( 0 2 1 3 ).\n");
    percolate_bond(&c_4, &edge_ids, true);

    igraph_vector_int_destroy(&edge_ids);
    printf("Zachary karate graph, no edge list given.\n");
    // Karate graph, no edge list given
    percolate_bond(&karate, NULL, false);
    // Generated disconnected graph, 100 vertices, p=0.01

    igraph_vector_int_t storage_order;
    igraph_vector_int_init_range(&storage_order, 0, igraph_ecount(&karate));
    printf("Zachary karate graph, edges in storage order.\n");
    percolate_bond(&karate, &storage_order, true);
    igraph_vector_int_destroy(&storage_order);

    percolate_bond(&random, NULL, false); // sanity check

    printf("Error on duplicates\n");
    igraph_vector_int_t duplicates;
    igraph_vector_int_init_int(&duplicates, 5, 0,1,2,3,3);
    CHECK_ERROR(percolate_bond(&c_4, &duplicates, false), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&duplicates);
    printf("Null outputs\n");
    igraph_bond_percolation(&karate, NULL, NULL, NULL);

    igraph_destroy(&singleton);
    igraph_destroy(&null_graph);
    igraph_destroy(&k_3);
    igraph_destroy(&c_4);
    igraph_destroy(&karate);
    igraph_destroy(&random);

    VERIFY_FINALLY_STACK();
}

igraph_error_t percolate_site(igraph_t *graph, igraph_vector_int_t *vert_indices, igraph_bool_t printing) {
    igraph_vector_int_t outputs, edge_counts;
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t number_percolated;

    if (vert_indices != NULL) {
        number_percolated = igraph_vector_int_size(vert_indices);
    }
    else {
        number_percolated = vcount;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edge_counts, 0);

    IGRAPH_CHECK(igraph_site_percolation(graph, &outputs, &edge_counts, vert_indices));

    if (printing) {
        print_vector_int(&outputs);
        print_vector_int(&edge_counts);
    }

    igraph_vector_int_t component_sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&component_sizes, 0);

    IGRAPH_CHECK(igraph_connected_components(graph, NULL, &component_sizes, NULL, IGRAPH_WEAK));

    IGRAPH_ASSERT(igraph_vector_int_size(&outputs) == number_percolated);
    if (number_percolated > 1 && number_percolated == vcount) {
        // it is only guaranteed to have the same component vcount if all vertices are included
        IGRAPH_ASSERT(igraph_vector_int_max(&outputs) == igraph_vector_int_max(&component_sizes));
    }
    igraph_vector_int_destroy(&component_sizes);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_int_t prev = 0;
    for (igraph_int_t i = 0; i < number_percolated; i++) {
        IGRAPH_ASSERT(VECTOR(outputs)[i] <= vcount); //vcount cannot be greater than number of vertices
        IGRAPH_ASSERT(VECTOR(outputs)[i] > 0);      // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(outputs)[i] >= prev);   // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(outputs)[i] <= i + 1); // Largest component cannot be bigger than a tree with the same number of vertices.
        prev = VECTOR(outputs)[i];
    }
    prev = 0;
    for (igraph_int_t i = 0; i < number_percolated; i++) {
        IGRAPH_ASSERT(VECTOR(edge_counts)[i] >= 0);      // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(edge_counts)[i] >= prev);   // Size of largest component must be nondecreasing.
        prev = VECTOR(edge_counts)[i];
    }

    igraph_vector_int_destroy(&edge_counts);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_vector_int_destroy(&outputs);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}


void test_site(void) {

    // Test with normal graph and provided edge list
    igraph_t k_5, c_4, karate, random, null_graph, singleton;

    igraph_empty(&null_graph, 0, false);
    igraph_empty(&singleton, 1, false);
    igraph_full(&k_5, 5, false, false);
    igraph_small(&c_4, 4, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_famous(&karate, "Zachary");
    igraph_erdos_renyi_game_gnp(&random, 100, 0.01, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    printf("# Site percolation test suite\n");
    printf("Null graph, no provided vertex order.\n");
    percolate_site(&null_graph, NULL, true);

    printf("Singleton graph, no provided vertex order.\n");
    percolate_site(&singleton, NULL, true);

    printf("K_5 graph percolation curve, no provided vertex sequence.\n");
    percolate_site(&k_5, NULL, true); // sequence is random, but since it should be consistent it doesn't matter

    igraph_vector_int_t vertex_ids;
    igraph_vector_int_init_int(&vertex_ids, 4, 0, 2, 1, 3);

    printf("C_4 graph with vertex sequence ( 0 2 1 3 ).\n");
    percolate_site(&c_4, &vertex_ids, true);

    igraph_vector_int_destroy(&vertex_ids);
    printf("Zachary karate graph, no vertex list given.\n");
    // Karate graph, no vertex list given
    percolate_site(&karate, NULL, false);
    // Generated disconnected graph, 100 vertices, p=0.01

    igraph_vector_int_t storage_order;
    igraph_vector_int_init_range(&storage_order, 0, igraph_vcount(&karate));
    printf("Zachary karate graph, vertices in storage order.\n");
    percolate_site(&karate, &storage_order, true);
    igraph_vector_int_destroy(&storage_order);

    percolate_site(&random, NULL, false);

    igraph_vector_int_t bad_vert_list_repeat, bad_vert_list_too_big, vert_list_missing_vertex;
    igraph_vector_int_init_int(&bad_vert_list_too_big, 6, 0, 1, 2, 3, 4, 5);
    igraph_vector_int_init_int(&vert_list_missing_vertex, 3, 0, 1, 2);
    igraph_vector_int_init_int(&bad_vert_list_repeat,  5, 0, 0, 0, 0, 0);

    // should error due to vertex being too big
    printf("K_5 with too big vertex index\n");
    CHECK_ERROR(percolate_site(&k_5, &bad_vert_list_too_big, false), IGRAPH_EINVVID);

    // missing vertices are allowed
    printf("K_5 with missing vertices\n");
    percolate_site(&k_5, &vert_list_missing_vertex, true);

    // should error due to repeated vertices
    printf("K_5 with repeated vertices\n");
    CHECK_ERROR(percolate_site(&k_5, &bad_vert_list_repeat, false),  IGRAPH_EINVAL);

    printf("Null outputs\n");
    igraph_bond_percolation(&karate, NULL, NULL, NULL);

    igraph_destroy(&singleton);
    igraph_destroy(&null_graph);
    igraph_destroy(&k_5);
    igraph_destroy(&c_4);
    igraph_destroy(&karate);
    igraph_destroy(&random);

    igraph_vector_int_destroy(&vert_list_missing_vertex);
    igraph_vector_int_destroy(&bad_vert_list_repeat);
    igraph_vector_int_destroy(&bad_vert_list_too_big);

    VERIFY_FINALLY_STACK();
}

igraph_error_t percolate_edgelist(igraph_vector_int_t *edges, igraph_bool_t printing) {
    igraph_vector_int_t giant_size, vertex_count;
    igraph_t graph;
    igraph_int_t ecount;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&giant_size, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_count, 0);

    IGRAPH_CHECK(igraph_edgelist_percolation(edges, &giant_size, &vertex_count));

    if (printing) {
        print_vector_int(&giant_size);
        print_vector_int(&vertex_count);
    }

    igraph_vector_int_t component_sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&component_sizes, 0);

    IGRAPH_CHECK(igraph_create(&graph, edges, 0, false));
    ecount = igraph_ecount(&graph);
    IGRAPH_CHECK(igraph_connected_components(&graph,NULL, &component_sizes,NULL,IGRAPH_WEAK));
    igraph_destroy(&graph);

    IGRAPH_ASSERT(igraph_vector_int_size(&giant_size) == ecount);
    if (ecount > 1) {
        IGRAPH_ASSERT(igraph_vector_int_max(&giant_size) == igraph_vector_int_max(&component_sizes));
    }

    igraph_int_t prev = 0;
    for (igraph_int_t i = 0; i < ecount; i++) {
        // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(giant_size)[i] > 0);
        // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(giant_size)[i] >= prev);
        // Largest component cannot be bigger than a tree with the same number of edges.
        IGRAPH_ASSERT(VECTOR(giant_size)[i] <= i + 2);
        prev = VECTOR(giant_size)[i];
    }

    prev = 0;
    for (igraph_int_t i = 0; i < ecount; i++) {
        // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(vertex_count)[i] > 0);
        // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(vertex_count)[i] >= prev);
        // Largest component cannot be bigger than a tree with the same number of edges.
        IGRAPH_ASSERT(VECTOR(vertex_count)[i] <= 2*i + 2);
        prev = VECTOR(vertex_count)[i];
    }

    igraph_vector_int_destroy(&giant_size);
    igraph_vector_int_destroy(&component_sizes);
    igraph_vector_int_destroy(&vertex_count);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

void test_edgelist_percolation(void) {
    // Edge list percolation is already called from bond percolation,
    // so this mostly tests for expected errors that cannot occur from generated edge lists.
    igraph_vector_int_t odd, negative, loopy;
    igraph_vector_int_init_int(&odd, 5, 0, 1, 1, 2, 1);
    igraph_vector_int_init_int(&negative, 6, -1, 1, 0, 0, 1, -1);
    igraph_vector_int_init_int(&loopy, 8, 0, 0, 0, 1, 1, 2, 2, 0);

    printf("# Edge list percolation\n");
    printf("Percolation with ( 0 1 1 2 1 ), odd number of entries\n");
    CHECK_ERROR(percolate_edgelist(&odd, false), IGRAPH_EINVAL);
    printf("Percolation with ( -1 1 0 0 1 -1 ), negative numbers\n");
    CHECK_ERROR(percolate_edgelist(&negative, false), IGRAPH_EINVVID);

    printf("Percolation with ( 0 0 0 1 1 2 2 0 ), k_3 with loop\n");
    percolate_edgelist(&loopy, true);

    printf("Null outputs\n");
    igraph_edgelist_percolation(&loopy, NULL, NULL);

    igraph_vector_int_destroy(&odd);
    igraph_vector_int_destroy(&negative);
    igraph_vector_int_destroy(&loopy);
    VERIFY_FINALLY_STACK();
}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 30062025);

    test_bond();
    test_site();
    test_edgelist_percolation();

    return 0;
}
