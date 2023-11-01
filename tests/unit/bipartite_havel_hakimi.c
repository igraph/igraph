//
// Created by Lara Holm on 1.11.2023.
//
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.h"
#include "math/safe_intop.h"

int main(void) {
    //igraph_t g;
    igraph_vector_int_t ds1, ds2;
    igraph_vector_int_t edges;
    igraph_integer_t ds1_sum;
    igraph_integer_t ds2_sum;

    igraph_integer_t deg1[] = {2, 3, 2, 1};
    igraph_integer_t deg2[] = {3, 1, 2, 1, 1};

    igraph_vector_int_init_array(&ds1, deg1, 4);
    igraph_vector_int_init_array(&ds2, deg2, 5);

    igraph_i_safe_vector_int_sum(&ds1, &ds1_sum);
    igraph_i_safe_vector_int_sum(&ds2, &ds2_sum);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ds1_sum+ds2_sum);

    igraph_vector_int_print(&ds1);
    igraph_vector_int_print(&ds2);


    igraph_realize_bipartite_degree_sequence(&ds1, &ds2, &edges, false);

    igraph_vector_int_print(&edges);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_vector_int_destroy(&edges);

    IGRAPH_FINALLY_CLEAN(1);
    VERIFY_FINALLY_STACK();

    //ds1 = [2, 3, 2, 1]
    //ds2 = [3, 1, 2, 1, 1]
    return 0;
}
