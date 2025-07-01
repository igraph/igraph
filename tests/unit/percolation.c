/*
    IGraph library.
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

igraph_error_t percolate_b(igraph_t *graph, igraph_vector_int_t *edge_indices, igraph_bool_t printing) {
    igraph_vector_int_t outputs;
    igraph_integer_t size = igraph_ecount(graph);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);

    IGRAPH_CHECK(igraph_bond_percolation(graph, &outputs, edge_indices));

    if (printing) print_vector_int(&outputs);

    igraph_vector_int_t components;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&components, 0);

    IGRAPH_CHECK(igraph_connected_components(graph, NULL, &components, NULL, IGRAPH_WEAK));

    IGRAPH_ASSERT(igraph_vector_int_size(&outputs) == size);
    if (size > 1) {
        IGRAPH_ASSERT(igraph_vector_int_max(&outputs) == igraph_vector_int_max(&components));
    }
    igraph_vector_int_destroy(&components);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_integer_t prev = 0;
    for (igraph_integer_t i = 0; i < size; i++) {
        IGRAPH_ASSERT(VECTOR(outputs)[i] > 0);      // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(outputs)[i] >= prev);   // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(outputs)[i] <= i + 2); // Largest component cannot be bigger than a tree with the same number of edges.
        prev = VECTOR(outputs)[i];
    }

    igraph_vector_int_destroy(&outputs);
    IGRAPH_FINALLY_CLEAN(1);

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
    igraph_erdos_renyi_game_gnp(&random, 100, 0.01, false, false);
    printf("# Bond percolation test suite\n");
    printf("Null graph, no provided edge order.\n");
    percolate_b(&null_graph, NULL, true);

    printf("Singleton graph, no provided edge order.\n");
    percolate_b(&singleton, NULL, true);

    printf("K_3 graph percolation curve, no provided edge sequence.\n");
    percolate_b(&k_3, NULL, true); // sequence is random, but since it should be consistent it doesn't matter

    igraph_vector_int_t edge_ids;
    igraph_vector_int_init_int(&edge_ids, 4, 0, 2, 1, 3);

    printf("C_4 graph with edge sequence ( 0 2 1 3 ).\n");
    percolate_b(&c_4, &edge_ids, true);

    igraph_vector_int_destroy(&edge_ids);
    printf("Zachary karate graph, no edge list given.\n");
    // Karate graph, no edge list given
    percolate_b(&karate, NULL, false);
    // Generated disconnected graph, 100 vertices, p=0.01

    igraph_vector_int_t storage_order;
    igraph_vector_int_init_range(&storage_order, 0, igraph_ecount(&karate));
    printf("Zachary karate graph, edges in storage order.\n");
    percolate_b(&karate, &storage_order, true);
    igraph_vector_int_destroy(&storage_order);

    percolate_b(&random, NULL, false);
    igraph_destroy(&singleton);
    igraph_destroy(&null_graph);
    igraph_destroy(&k_3);
    igraph_destroy(&c_4);
    igraph_destroy(&karate);
    igraph_destroy(&random);
    VERIFY_FINALLY_STACK();

}

igraph_error_t percolate_s(igraph_t *graph, igraph_vector_int_t *vert_indices, igraph_bool_t printing) {
    igraph_vector_int_t outputs;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);
    igraph_integer_t size = igraph_vcount(graph);
    IGRAPH_CHECK(igraph_site_percolation(graph, &outputs, vert_indices));

    if (printing) print_vector_int(&outputs);

    igraph_vector_int_t components;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&components, 0);

    IGRAPH_CHECK(igraph_connected_components(graph, NULL, &components, NULL, IGRAPH_WEAK));

    IGRAPH_ASSERT(igraph_vector_int_size(&outputs) == size);
    if (size > 1) {
        IGRAPH_ASSERT(igraph_vector_int_max(&outputs) == igraph_vector_int_max(&components));
    }
    igraph_vector_int_destroy(&components);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_integer_t prev = 0;
    for (igraph_integer_t i = 0; i < size; i++) {
        IGRAPH_ASSERT(VECTOR(outputs)[i] > 0);      // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(outputs)[i] >= prev);   // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(outputs)[i] <= i + 1); // Largest component cannot be bigger than a tree with the same number of vertices.
        prev = VECTOR(outputs)[i];
    }


    igraph_vector_int_destroy(&outputs);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

igraph_error_t largest_component_s(igraph_t *graph, igraph_integer_t *size) {
    igraph_vector_int_t outputs;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);

    IGRAPH_CHECK(igraph_site_percolation(graph, &outputs, NULL));

    *size = VECTOR(outputs)[igraph_vector_int_size(&outputs) -1];

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
    igraph_erdos_renyi_game_gnp(&random, 100, 0.01, false, false);
    printf("# Site percolation test suite\n");
    printf("Null graph, no provided vertex order.\n");
    percolate_s(&null_graph, NULL, true);

    printf("Singleton graph, no provided vertex order.\n");
    percolate_s(&singleton, NULL, true);

    printf("K_5 graph percolation curve, no provided vertex sequence.\n");
    percolate_s(&k_5, NULL, true); // sequence is random, but since it should be consistent it doesn't matter

    igraph_vector_int_t vertex_ids;
    igraph_vector_int_init_int(&vertex_ids, 4, 0, 2, 1, 3);

    printf("C_4 graph with vertex sequence ( 0 2 1 3 ).\n");
    percolate_s(&c_4, &vertex_ids, true);

    igraph_vector_int_destroy(&vertex_ids);
    printf("Zachary karate graph, no vertex list given.\n");
    // Karate graph, no vertex list given
    percolate_s(&karate, NULL, false);
    // Generated disconnected graph, 100 vertices, p=0.01

    igraph_vector_int_t storage_order;
    igraph_vector_int_init_range(&storage_order, 0, igraph_vcount(&karate));
    printf("Zachary karate graph, vertices in storage order.\n");
    percolate_s(&karate, &storage_order, true);
    igraph_vector_int_destroy(&storage_order);

    percolate_s(&random, NULL, false);

    igraph_vector_int_t bad_vert_list_repeat, bad_vert_list_too_big, bad_vert_list_missing;
    igraph_vector_int_init_int(&bad_vert_list_too_big, 6, 0, 1, 2, 3, 4, 5);
    igraph_vector_int_init_int(&bad_vert_list_missing, 3, 0, 1, 2);
    igraph_vector_int_init_int(&bad_vert_list_repeat,  5, 0, 0, 0, 0, 0);
    // should error due to being too big
    printf("K_5 with too big vertex list\n");
    CHECK_ERROR(percolate_s(&k_5, &bad_vert_list_too_big, false), IGRAPH_EINVAL);
    // should error due to being too small
    printf("K_5 with too small vertex list\n");
    CHECK_ERROR(percolate_s(&k_5, &bad_vert_list_missing, false), IGRAPH_EINVAL);
    // should error due to repeated vertices
    printf("K_5 with repeated vertices\n");
    CHECK_ERROR(percolate_s(&k_5, &bad_vert_list_repeat, false),  IGRAPH_EINVAL);


    igraph_destroy(&singleton);
    igraph_destroy(&null_graph);
    igraph_destroy(&k_5);
    igraph_destroy(&c_4);
    igraph_destroy(&karate);
    igraph_destroy(&random);

    igraph_vector_int_destroy(&bad_vert_list_missing);
    igraph_vector_int_destroy(&bad_vert_list_repeat);
    igraph_vector_int_destroy(&bad_vert_list_too_big);

    VERIFY_FINALLY_STACK();
}


igraph_error_t el_percolate(igraph_vector_int_t * edge_list, igraph_bool_t printing) {
    igraph_vector_int_t outputs;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);
    IGRAPH_CHECK(igraph_edge_list_percolation(edge_list, &outputs));

    if (printing) print_vector_int(&outputs);

    igraph_vector_int_t components;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&components, 0);

    igraph_t graph;


    IGRAPH_CHECK(igraph_create(&graph, edge_list, 0, false));

    igraph_integer_t size = igraph_ecount(&graph);

    IGRAPH_CHECK(igraph_connected_components(&graph, NULL, &components, NULL, IGRAPH_WEAK));
    igraph_destroy(&graph);
    IGRAPH_ASSERT(igraph_vector_int_size(&outputs) == size);
    if (size > 1) {
        IGRAPH_ASSERT(igraph_vector_int_max(&outputs) == igraph_vector_int_max(&components));
    }
    igraph_vector_int_destroy(&components);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_integer_t prev = 0;
    for (igraph_integer_t i = 0; i < size; i++) {
        IGRAPH_ASSERT(VECTOR(outputs)[i] > 0);      // Sizes cannot be negative.
        IGRAPH_ASSERT(VECTOR(outputs)[i] >= prev);   // Size of largest component must be nondecreasing.
        IGRAPH_ASSERT(VECTOR(outputs)[i] <= i + 2); // Largest component cannot be bigger than a tree with the same number of edges.
        prev = VECTOR(outputs)[i];
    }

    igraph_vector_int_destroy(&outputs);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

void test_edge_list_percolation(void) {
    // Edge list percolation is already called from bond percolation,
    // so this mostly tests for expected errors that cannot occur from generated edge lists.
    igraph_vector_int_t odd, negative;
    igraph_vector_int_init_int(&odd, 5, 0,1,1,2,1);
    igraph_vector_int_init_int(&negative, 6, -1,1,0,0,1,-1);

    printf("# Edge list percolation\n");
    printf("Percolation with ( 0 1 1 2 1 ), odd number of entries\n");
    CHECK_ERROR(el_percolate(&odd, false), IGRAPH_EINVAL);
    printf("Percolation with ( -1 1 0 0 1 -1 ), negative numbers\n");
    CHECK_ERROR(el_percolate(&negative, false), IGRAPH_EINVVID);

    igraph_vector_int_destroy(&odd);
    igraph_vector_int_destroy(&negative);
    VERIFY_FINALLY_STACK();

}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 30062025);
    test_bond();
    test_site();
    test_edge_list_percolation();
    return 0;
}
