/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "test_utilities.inc"

int test_from_prufer_back_to_prufer() {
    igraph_t graph;
    igraph_integer_t prufer[] = {2, 3, 2, 3};

    igraph_vector_int_t expected_prufer, output_prufer;

    igraph_bool_t success = 0;

    igraph_vector_int_view(&expected_prufer, prufer, 4);
    IGRAPH_CHECK(igraph_from_prufer(&graph, &expected_prufer));

    IGRAPH_CHECK(igraph_vector_int_init(&output_prufer, 4));
    IGRAPH_CHECK(igraph_to_prufer(&graph, &output_prufer));

    success = igraph_vector_int_all_e(&expected_prufer, &output_prufer);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&output_prufer);

    return success;
}

int test_from_prufer_back_to_prufer_with_resize() {
    igraph_t graph;
    igraph_integer_t prufer[] = {0, 2, 4, 1, 1, 0};

    igraph_vector_int_t expected_prufer, output_prufer;

    igraph_bool_t success;

    igraph_vector_int_view(&expected_prufer, prufer, 6);
    IGRAPH_CHECK(igraph_from_prufer(&graph, &expected_prufer));

    IGRAPH_CHECK(igraph_vector_int_init(&output_prufer, 0));
    IGRAPH_CHECK(igraph_to_prufer(&graph, &output_prufer));

    success = igraph_vector_int_all_e(&expected_prufer, &output_prufer);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&output_prufer);

    return success;
}

int test_from_prufer_back_to_prufer_with_resize2() {
    igraph_t graph;
    igraph_integer_t prufer[] = {2, 4, 5, 1, 3};

    igraph_vector_int_t expected_prufer, output_prufer;

    igraph_bool_t success;

    igraph_vector_int_view(&expected_prufer, prufer, 5);
    IGRAPH_CHECK(igraph_from_prufer(&graph, &expected_prufer));

    IGRAPH_CHECK(igraph_vector_int_init(&output_prufer, 0));
    IGRAPH_CHECK(igraph_to_prufer(&graph, &output_prufer));


    success = igraph_vector_int_all_e(&output_prufer, &expected_prufer);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&output_prufer);

    return success;
}

int random_tree(int size, igraph_t* tree, igraph_vector_int_t* prufer) {
    int i, j;
    int prufer_length;

    if (size < 0) {
        return IGRAPH_EINVAL;
    }

    if (size < 2) {
        return igraph_empty(tree, size, IGRAPH_UNDIRECTED);
    }

    prufer_length = size - 2;
    IGRAPH_CHECK(igraph_vector_int_resize(prufer, prufer_length));

    for (i = 0; i < prufer_length; ++i) {
        j = RNG_INTEGER(0, size - 1);
        VECTOR(*prufer)[i] = j;
    }

    IGRAPH_CHECK(igraph_from_prufer(tree, prufer));

    return IGRAPH_SUCCESS;
}

int test_from_random_prufer_back_to_prufer(int tree_size) {
    igraph_t graph;
    igraph_vector_int_t expected_prufer, output_prufer;

    igraph_bool_t success = 0;
    igraph_integer_t random_seed = 4096;

    IGRAPH_CHECK(igraph_vector_int_init(&output_prufer, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&expected_prufer, 0));

    igraph_rng_seed(igraph_rng_default(), random_seed);

    IGRAPH_CHECK(random_tree(tree_size, &graph, &expected_prufer));

    IGRAPH_CHECK(igraph_to_prufer(&graph, &output_prufer));

    success = igraph_vector_int_all_e(&output_prufer, &expected_prufer);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&expected_prufer);
    igraph_vector_int_destroy(&output_prufer);

    return success;
}

#undef RUN_TEST   /* from test_utilities.inc */

int test_num = 0;
#define RUN_TEST(TEST) \
    test_num++; \
    if(!(TEST)) { \
        return test_num; \
    }

int main() {
    RUN_TEST(test_from_prufer_back_to_prufer());
    RUN_TEST(test_from_prufer_back_to_prufer_with_resize());
    RUN_TEST(test_from_prufer_back_to_prufer_with_resize2());
    RUN_TEST(test_from_random_prufer_back_to_prufer(10));
    RUN_TEST(test_from_random_prufer_back_to_prufer(100));
    RUN_TEST(test_from_random_prufer_back_to_prufer(1000));
    RUN_TEST(test_from_random_prufer_back_to_prufer(10000));

    VERIFY_FINALLY_STACK();

    return 0;
}
