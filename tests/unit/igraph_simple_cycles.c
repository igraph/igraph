/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

void check_cycles_max(const igraph_t *graph, igraph_neimode_t mode, igraph_int_t expected, igraph_int_t max_cycle_length) {
    igraph_vector_int_list_t results_v;
    igraph_vector_int_list_t results_e;

    igraph_vector_int_list_init(&results_v, 0);
    igraph_vector_int_list_init(&results_e, 0);

    igraph_simple_cycles(graph, &results_v, &results_e, mode, 0, max_cycle_length, IGRAPH_UNLIMITED);

    printf("Finished search, found %" IGRAPH_PRId
           " cycles, expected %" IGRAPH_PRId " cycles."
           " Max cycle length was %" IGRAPH_PRId " vertices.\n\n",
           igraph_vector_int_list_size(&results_v), expected, max_cycle_length);

    if (igraph_vcount(graph) < 100) {
        printf("Vertex IDs in cycles:\n");
        print_vector_int_list(&results_v);

        printf("Edge IDs in cycles:\n");
        print_vector_int_list(&results_e);
    }

    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_v) == expected);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_e) == expected);

    for (igraph_int_t i = 0; i < expected; i++) {
        igraph_vector_int_t *vertices, *edges;

        vertices = igraph_vector_int_list_get_ptr(&results_v, i);
        edges = igraph_vector_int_list_get_ptr(&results_e, i);
        IGRAPH_ASSERT(igraph_vector_int_size(vertices) == igraph_vector_int_size(edges));
        if (max_cycle_length >= 0) {
            IGRAPH_ASSERT(igraph_vector_int_size(vertices) <= max_cycle_length);
        }
    }

    igraph_vector_int_list_destroy(&results_v);
    igraph_vector_int_list_destroy(&results_e);
}

void check_cycles(const igraph_t *graph, igraph_neimode_t mode, igraph_int_t expected) {
    check_cycles_max(graph, mode, expected, -1);
}

int main(void) {
    igraph_t g;
    igraph_t g_ring_undirected, g_star_undirected;

    printf("Testing null graph\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    check_cycles(&g, IGRAPH_OUT, 0);
    igraph_destroy(&g);

    printf("\nTesting edgeless graph\n");
    igraph_empty(&g, 5, IGRAPH_UNDIRECTED);
    check_cycles(&g, IGRAPH_OUT, 0);
    igraph_destroy(&g);

    printf("\nTesting directed cycle graph\n");
    igraph_ring(&g, 10, IGRAPH_DIRECTED, /*mutual=*/ false, /*circular=*/ true);
    check_cycles(&g, IGRAPH_OUT, 1);
    igraph_destroy(&g);

    // printf("\nTesting large directed cycle graph\n");
    // igraph_ring(&g, 10000, IGRAPH_DIRECTED, /*mutual=*/ false, /*circular=*/ true);
    // check_cycles(&g, 1);
    // igraph_destroy(&g);

    printf("\nTesting directed star\n");
    igraph_star(&g, 7, IGRAPH_STAR_OUT, 1);
    check_cycles(&g, IGRAPH_OUT, 0);
    igraph_destroy(&g);

    printf("\nTesting directed wheel\n");
    igraph_wheel(&g, 8, IGRAPH_WHEEL_OUT, 0);
    check_cycles(&g, IGRAPH_OUT, 1);
    check_cycles(&g, IGRAPH_IN, 1);
    check_cycles(&g, IGRAPH_ALL, 43);
    igraph_destroy(&g);

    printf("\nTesting undirected ring\n");
    igraph_ring(&g_ring_undirected, 10, IGRAPH_UNDIRECTED, /*mutual=*/ false, /*circular=*/ true);
    check_cycles(&g_ring_undirected, IGRAPH_OUT, 1);

    igraph_star(&g_star_undirected, 7, IGRAPH_STAR_UNDIRECTED, 1);
    printf("\nTesting undirected star\n");
    check_cycles(&g_star_undirected, IGRAPH_OUT, 0);

    igraph_disjoint_union(&g, &g_ring_undirected, &g_star_undirected);
    printf("\nTesting union of undirected wheel and star\n");
    check_cycles(&g, IGRAPH_OUT, 1);

    printf("\nTesting union of undirected wheel, star and a single edge\n");
    igraph_add_edge(&g, 7, 13); // add a random edge between the two structures to make them connected
    check_cycles(&g, IGRAPH_OUT, 1);
    igraph_destroy(&g);
    igraph_destroy(&g_star_undirected);
    igraph_destroy(&g_ring_undirected);

    printf("\nTesting a tree\n");
    igraph_kary_tree(&g, 20, 3, IGRAPH_TREE_OUT);
    check_cycles(&g, IGRAPH_OUT, 0);
    check_cycles(&g, IGRAPH_ALL, 0);
    igraph_destroy(&g);

    printf("\nTesting a complete DAG\n");
    igraph_full_citation(&g, 5, IGRAPH_DIRECTED);
    check_cycles(&g, IGRAPH_OUT, 0);
    check_cycles(&g, IGRAPH_ALL, 37);
    igraph_destroy(&g);

    igraph_t g_wheel_undirected_2;
    igraph_wheel(&g_wheel_undirected_2, 10, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("\nCreated undirected wheel\n");
    // call cycles finder, expect 73 cycles to be found (
    // 9 cycle of 10 nodes,
    // 10 cycles of 9 nodes,
    // 9 cycles of 8 nodes,
    // 9 cycles of 7 nodes,
    // 9 cycles of 6 nodes,
    // 9 cycles of 5 nodes,
    // 9 cycles of 4 nodes
    // 9 cycles of 3 nodes,
    // )
    check_cycles(&g_wheel_undirected_2, IGRAPH_OUT, 73);
    // test the max_cycle_length parameter
    check_cycles_max(&g_wheel_undirected_2, IGRAPH_OUT, 9, 3);
    check_cycles_max(&g_wheel_undirected_2, IGRAPH_OUT, 18, 4);
    // clean up
    igraph_destroy(&g_wheel_undirected_2);


    // Tests as requested in https://github.com/igraph/igraph/pull/2181#issuecomment-1326064152
    /*
     * This graph looks like:
     *
     * 1--2\
     * |  | |
     * 4--3/
     *
     */
    printf("\nTesting directed graph with a cycle of length 4 and a multi-edge\n");
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 1, 2,
                 2, 3,
                 2, 3,
                 3, 4,
                 4, 1,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 2);
    check_cycles_max(&g, IGRAPH_OUT, 2, 4);
    igraph_destroy(&g);

    // same, but undirected
    printf("\nTesting undirected graph with a cycle of length 4 and a multi-edge\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 1, 2,
                 2, 3,
                 2, 3,
                 3, 4,
                 4, 1,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 3);
    check_cycles_max(&g, IGRAPH_OUT, 3, 4);
    igraph_destroy(&g);

    // check that self-loops are handled
    printf("\nTesting directed graph with single self-loop\n");
    igraph_small(&g, 1, IGRAPH_DIRECTED, 0, 0, -1);
    check_cycles(&g, IGRAPH_OUT, 1);
    check_cycles(&g, IGRAPH_ALL, 1);
    igraph_destroy(&g);

    // Tests as requested in https://github.com/igraph/igraph/pull/2181#issuecomment-2243053770
    printf("\nTesting undirected graph with single length-2-loop\n");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 1, 0, 1, -1);
    check_cycles(&g, IGRAPH_OUT, 1);
    igraph_destroy(&g);

    printf("\nTesting undirected graph with 3 length-2-loops\n");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 1, 0, 1, 0, 1, -1);
    check_cycles(&g, IGRAPH_OUT, 3);
    igraph_destroy(&g);

    // and in https://github.com/igraph/igraph/pull/2181#issuecomment-2243060608
    printf("\nTesting undirected graph with single length-2-loop and a self-loop\n");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 1, 0, 1, 0, 0, -1);
    check_cycles(&g, IGRAPH_OUT, 2);
    igraph_destroy(&g);


    // Tests as requested in https://github.com/igraph/igraph/pull/2181#issuecomment-2243942240
    printf("\nTesting undirected graph of type 'envelope'\n");
    // expect:
    // 4 of length 3 vertices
    // 4 of length 5 vertices
    // 5 of length 4 vertices
    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 0, 1,
                 0, 3,
                 0, 4,
                 1, 2,
                 1, 3,
                 2, 3,
                 2, 4,
                 3, 4,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 13);
    igraph_destroy(&g);

    printf("\nTesting undirected graph of type 'boat'\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 0, 2,
                 0, 4,
                 1, 2,
                 1, 3,
                 1, 4,
                 2, 3,
                 2, 4,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 6);
    igraph_destroy(&g);


    // Tests as requested in https://github.com/igraph/igraph/pull/2181#issuecomment-2249751754
    printf("\nTesting undirected graph of type 'house'\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 0, 3,
                 0, 4,
                 1, 2,
                 1, 3,
                 1, 4,
                 2, 3,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 3);
    igraph_destroy(&g);

    printf("\nTesting undirected graph of type 'prism'\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0, 1,
                 0, 3,
                 0, 5,
                 1, 4,
                 1, 5,
                 2, 3,
                 2, 4,
                 2, 5,
                 3, 4,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 14);
    igraph_destroy(&g);


    // Tests as requested in https://github.com/igraph/igraph/pull/2181#issuecomment-2251350410
    printf("\nTesting undirected graph of type '7 vertices'\n");
    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                 0, 1,
                 0, 4,
                 0, 5,
                 0, 6,
                 1, 2,
                 1, 3,
                 1, 5,
                 1, 6,
                 2, 3,
                 2, 6,
                 3, 4,
                 3, 5,
                 3, 6,
                 4, 5,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 89);
    igraph_destroy(&g);

    printf("\nTesting undirected graph of type 'shooting star'\n");
    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                 0, 1,
                 0, 2,
                 0, 3,
                 0, 4,
                 0, 6,
                 1, 2,
                 1, 3,
                 1, 4,
                 1, 5,
                 1, 6,
                 2, 3,
                 2, 4,
                 2, 6,
                 3, 4,
                 3, 5,
                 4, 6,
                 5, 6,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 18 + 43 + 78 + 96 + 60);
    igraph_destroy(&g);


    // Tests as requested in https://github.com/igraph/igraph/pull/2181#issuecomment-2428987492
    igraph_small(&g, 5, IGRAPH_DIRECTED,
        0, 1,
        1, 2,
        2, 3,
        3, 0,
        1, 4,
        4, 2,
        -1);
    check_cycles(&g, IGRAPH_OUT, 2);
    check_cycles_max(&g, IGRAPH_OUT, 0, 3);
    check_cycles_max(&g, IGRAPH_OUT, 1, 4);
    check_cycles_max(&g, IGRAPH_OUT, 2, 5);
    igraph_destroy(&g);


    printf("\nTesting directed graph of type 'stable boat'\n");
    igraph_small(&g, 5, IGRAPH_DIRECTED,
         0, 2,
         0, 3,
         0, 4,
         1, 0,
         2, 3,
         3, 4,
         4, 1,
         4, 3,
        -1);
    check_cycles(&g, IGRAPH_OUT, 4);
    check_cycles_max(&g, IGRAPH_OUT, 1, 2);
    check_cycles_max(&g, IGRAPH_OUT, 2, 3);
    check_cycles_max(&g, IGRAPH_OUT, 3, 4);
    igraph_destroy(&g);

    printf("\nTesting directed graph of type 'stable letter'\n");
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 2, 0, 3, 1, 2, 1, 4, 2, 3, 2, 4, 3, 1, 4, 0, 4, 2,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 7);
    check_cycles_max(&g, IGRAPH_OUT, 0, 1);
    check_cycles_max(&g, IGRAPH_OUT, 1, 2);
    check_cycles_max(&g, IGRAPH_OUT, 3, 3);
    check_cycles_max(&g, IGRAPH_OUT, 5, 4);
    check_cycles_max(&g, IGRAPH_OUT, 7, 5);
    igraph_destroy(&g);

    printf("\nTesting directed graph of type 'double square'\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 2, 0, 4, 1, 3, 2, 5, 3, 0, 4, 1, 5, 4,
                 -1);
    check_cycles(&g, IGRAPH_OUT, 2);
    check_cycles_max(&g, IGRAPH_OUT, 0, 1);
    check_cycles_max(&g, IGRAPH_OUT, 0, 2);
    check_cycles_max(&g, IGRAPH_OUT, 0, 3);
    check_cycles_max(&g, IGRAPH_OUT, 1, 4);
    check_cycles_max(&g, IGRAPH_OUT, 1, 5);
    check_cycles_max(&g, IGRAPH_OUT, 2, 6);
    igraph_destroy(&g);

    printf("\nTesting undirected graph of type 'Mickey'\n");
    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                     0,1,
                     1,2,
                     2,0,
                     0,0,
                     0,3,
                     3,4,
                     4,5,
                     5,0,
                     -1);
    check_cycles(&g, IGRAPH_ALL, 3);
    igraph_destroy(&g);

    printf("\nTesting undirected graph of type 'Mickey2'\n");
    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                     0,1,
                     1,2,
                     2,0,
                     1,1,
                     0,3,
                     3,4,
                     4,5,
                     5,0,
                     -1);
    check_cycles(&g, IGRAPH_ALL, 3);
    igraph_destroy(&g);
    igraph_small(&g, 7, IGRAPH_DIRECTED,
                     0,1,
                     1,2,
                     2,0,
                     1,1,
                     0,3,
                     3,4,
                     4,5,
                     5,0,
                     -1);
    check_cycles(&g, IGRAPH_ALL, 3);
    igraph_destroy(&g);

    // as requested in https://github.com/igraph/igraph/issues/2692#issuecomment-2457627378
    printf("\nTesting undirected graph of type 'Mickey3'\n");
    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                     0, 1,
                     1, 2,
                     2, 0,
                     0, 3,
                     3, 4,
                     4, 5,
                     5, 0,
                     1, 1,
                     5, 6,
                     6, 5,
                     5, 5,
                     5, 5,
                     -1);
    check_cycles(&g, IGRAPH_ALL, 6);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    {
        igraph_vector_int_list_t v;

        igraph_vector_int_list_init(&v, 0);

        igraph_small(&g, 7, IGRAPH_DIRECTED,
                     0,1, 1,2, 2,0,
                     0,0,
                     0,3, 3,4, 4,5, 5,0,
                     -1);

        // Check that passing NULL vector lists doesn't crash.
        igraph_simple_cycles(&g, NULL, NULL, IGRAPH_ALL, -1, -1, IGRAPH_UNLIMITED);

        // Test limit on minimum cycle size.
        igraph_simple_cycles(&g, &v, NULL, IGRAPH_ALL, -1, -1, IGRAPH_UNLIMITED);
        IGRAPH_ASSERT(igraph_vector_int_list_size(&v) == 3);

        igraph_simple_cycles(&g, NULL, &v, IGRAPH_ALL, 1, -1, IGRAPH_UNLIMITED);
        IGRAPH_ASSERT(igraph_vector_int_list_size(&v) == 3);

        igraph_simple_cycles(&g, &v, NULL, IGRAPH_ALL, 2, -1, IGRAPH_UNLIMITED);
        IGRAPH_ASSERT(igraph_vector_int_list_size(&v) == 2);

        igraph_simple_cycles(&g, NULL, &v, IGRAPH_ALL, 2, 3, IGRAPH_UNLIMITED);
        IGRAPH_ASSERT(igraph_vector_int_list_size(&v) == 1);

        igraph_destroy(&g);
        igraph_vector_int_list_destroy(&v);
    }

    return 0;
}
