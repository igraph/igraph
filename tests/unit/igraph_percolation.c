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

igraph_error_t percolate_b(igraph_t *graph, igraph_vector_int_t *edge_indices) {
    igraph_vector_int_t outputs;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);

    IGRAPH_CHECK(igraph_bond_percolation(graph, &outputs, edge_indices));

    print_vector_int(&outputs);

    igraph_vector_int_destroy(&outputs);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

igraph_error_t largest_component_b(igraph_t *graph, igraph_integer_t *size) {
    igraph_vector_int_t outputs;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);

    IGRAPH_CHECK(igraph_bond_percolation(graph, &outputs, NULL));

    *size = VECTOR(outputs)[igraph_vector_int_size(&outputs) -1];

    igraph_vector_int_destroy(&outputs);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;

}

int test_bond(void) {
<<<<<<< HEAD
    // test with normal graph and provided edge list
=======
    // Test with normal graph and provided edge list
>>>>>>> dff3c4c5e749f6cdf5b02b3ae5200717d6d3b5c7

    igraph_t k_3, c_4, karate, random;
    igraph_integer_t size = 0;
    igraph_full(&k_3, 3, false, false);

    igraph_small(&c_4, 4, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 0, -1);

    printf("K_3 graph percolation curve, no provided edge sequence:\n");
    IGRAPH_CHECK(percolate_b(&k_3, NULL));

    igraph_vector_int_t edge_ids;
    IGRAPH_CHECK(igraph_vector_int_init_int(&edge_ids, 4, 0, 2, 1, 3));

    printf("C_4 graph with edge sequence 0, 2\n");
    IGRAPH_CHECK(percolate_b(&c_4, &edge_ids));

    igraph_vector_int_destroy(&edge_ids);
    igraph_destroy(&k_3);
    igraph_destroy(&c_4);

    igraph_famous(&karate, "Zachary");
    // Karate graph biggest component size, no edge list given
    largest_component_b(&karate, &size);
    IGRAPH_ASSERT(size == 34);
    // Generated disconnected graph, 100 vertices, p=0.01

    igraph_destroy(&karate);

    IGRAPH_CHECK(igraph_erdos_renyi_game_gnp(&random, 100, 0.01, false, false));

    igraph_vector_int_t components;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&components, 0);

    IGRAPH_CHECK(igraph_connected_components(&random, NULL, &components, NULL, IGRAPH_WEAK));

<<<<<<< HEAD

=======
>>>>>>> dff3c4c5e749f6cdf5b02b3ae5200717d6d3b5c7
    largest_component_b(&random, &size);
    IGRAPH_ASSERT(size == igraph_vector_int_max(&components));
    igraph_vector_int_destroy(&components);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_destroy(&random);
    VERIFY_FINALLY_STACK();
<<<<<<< HEAD
    return 0;
}
=======

    return 0;
}

>>>>>>> dff3c4c5e749f6cdf5b02b3ae5200717d6d3b5c7
igraph_error_t percolate_s(igraph_t *graph, igraph_vector_int_t *vert_indices) {
    igraph_vector_int_t outputs;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outputs, 0);

    IGRAPH_CHECK(igraph_site_percolation(graph, &outputs, vert_indices));

    print_vector_int(&outputs);

    igraph_vector_int_destroy(&outputs);
    IGRAPH_FINALLY_CLEAN(1);
<<<<<<< HEAD
=======

>>>>>>> dff3c4c5e749f6cdf5b02b3ae5200717d6d3b5c7
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

int test_site(void) {
    // test with normal graph and provided edge list

    igraph_t k_5, c_4, karate, random;
    igraph_integer_t size = 0;
    igraph_full(&k_5, 5, false, false);

    igraph_small(&c_4, 4, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 0, -1);

    printf("K_5 graph percolation curve, no provided vertex sequence:\n");
    IGRAPH_CHECK(percolate_s(&k_5, NULL));

    igraph_vector_int_t edge_ids;
    IGRAPH_CHECK(igraph_vector_int_init_int(&edge_ids, 4, 0, 2, 1, 3));

    printf("C_4 graph with vertex sequence 0, 2\n");
    IGRAPH_CHECK(percolate_s(&c_4, &edge_ids));

    igraph_vector_int_destroy(&edge_ids);

    igraph_destroy(&c_4);

    igraph_famous(&karate, "Zachary");

    //Karate graph biggest component size, no vertex list given: ");
    largest_component_s(&karate, &size);
    IGRAPH_ASSERT(size == 34);

    // Generated disconnected graph, 100 vertices, p=0.01

    igraph_destroy(&karate);

    IGRAPH_CHECK(igraph_erdos_renyi_game_gnp(&random, 100, 0.01, false, false));

    igraph_vector_int_t components;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&components, 0);

    IGRAPH_CHECK(igraph_connected_components(&random, NULL, &components, NULL, IGRAPH_WEAK));


    largest_component_s(&random, &size);
    IGRAPH_ASSERT(size == igraph_vector_int_max(&components));

    igraph_vector_int_destroy(&components);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vector_int_t bad_vert_list_repeat, bad_vert_list_too_big, bad_vert_list_missing;
    igraph_vector_int_init_int(&bad_vert_list_too_big, 6, 0, 1, 2, 3, 4, 5);
    igraph_vector_int_init_int(&bad_vert_list_missing, 3, 0, 1, 2);
    igraph_vector_int_init_int(&bad_vert_list_repeat,  5, 0, 0, 0, 0, 0);
    // should error due to being too big
    CHECK_ERROR(percolate_s(&k_5, &bad_vert_list_too_big), IGRAPH_EINVAL);
    // should error due to being too small
    CHECK_ERROR(percolate_s(&k_5, &bad_vert_list_missing), IGRAPH_EINVAL);
    // should error due to repeated vertices
    CHECK_ERROR(percolate_s(&k_5, &bad_vert_list_repeat),  IGRAPH_EINVAL);


    igraph_destroy(&k_5);
    igraph_destroy(&random);
    igraph_vector_int_destroy(&bad_vert_list_missing);
    igraph_vector_int_destroy(&bad_vert_list_repeat);
    igraph_vector_int_destroy(&bad_vert_list_too_big);

    VERIFY_FINALLY_STACK();
    return 0;
}


int main(void) {
    test_bond();
    test_site();
    return 0;
}
