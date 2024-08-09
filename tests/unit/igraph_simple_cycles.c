/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "igraph_cycles.h"
#include "test_utilities.h"
#include <stdlib.h>

void check_cycles(const igraph_t *graph, igraph_integer_t expected) {
    igraph_vector_int_list_t results_v;
    igraph_vector_int_list_t results_e;
    igraph_integer_t i;

    igraph_vector_int_list_init(&results_v, 0);
    igraph_vector_int_list_init(&results_e, 0);

    igraph_simple_cycles_search_all(graph, &results_v, &results_e);

    printf("Finished search, found %" IGRAPH_PRId " cycles.\n\n", igraph_vector_int_list_size(&results_v));

    if (igraph_vcount(graph) < 100) {
        printf("Vertex IDs in cycles:\n");
        print_vector_int_list(&results_v);

        printf("Edge IDs in cycles:\n");
        print_vector_int_list(&results_e);
    }

    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_v) == expected);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_e) == expected);

    for (i = 0; i < expected; i++) {
        igraph_vector_int_t *vertices, *edges;

        vertices = igraph_vector_int_list_get_ptr(&results_v, i);
        edges = igraph_vector_int_list_get_ptr(&results_e, i);
        IGRAPH_ASSERT(igraph_vector_int_size(vertices) == igraph_vector_int_size(edges));
    }

    igraph_vector_int_list_destroy(&results_v);
    igraph_vector_int_list_destroy(&results_e);
}

int main(void) {
    igraph_t g, g_ring_undirected, g_star_undirected;

    printf("Testing null graph\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED, -1);
    check_cycles(&g, 0);
    igraph_destroy(&g);

    printf("\nTesting empty graph\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED, -1);
    check_cycles(&g, 0);
    igraph_destroy(&g);

    printf("\nTesting ring\n");
    igraph_ring(&g, 10, /*directed=*/ 1, /*mutual=*/ 0, /*circular=*/ 1);
    check_cycles(&g, 1);
    igraph_destroy(&g);

    printf("\nTesting large ring\n");
    igraph_ring(&g, 10000, /*directed=*/ 1, /*mutual=*/ 0, /*circular=*/ 1);
    check_cycles(&g, 1);
    igraph_destroy(&g);

    printf("\nTesting star\n");
    igraph_star(&g, 7, IGRAPH_STAR_OUT, 1);
    check_cycles(&g, 0);
    igraph_destroy(&g);

    printf("\nTesting directed wheel\n");
    igraph_wheel(&g, 10, IGRAPH_WHEEL_OUT, 0);
    check_cycles(&g, 1);
    igraph_destroy(&g);

    printf("\nTesting undirected ring\n");
    igraph_ring(&g_ring_undirected, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
    check_cycles(&g_ring_undirected, 1);

    igraph_star(&g_star_undirected, 7, IGRAPH_STAR_UNDIRECTED, 1);
    printf("\nTesting undirected star\n");
    check_cycles(&g_star_undirected, 0);

    igraph_disjoint_union(&g, &g_ring_undirected, &g_star_undirected);
    printf("\nTesting union of undirected wheel and star\n");
    check_cycles(&g, 1);

    printf("\nTesting union of undirected wheel, star and a single edge\n");
    igraph_add_edge(&g, 7, 13); // add a random edge between the two structures to make them connected
    check_cycles(&g, 1);
    igraph_destroy(&g);
    igraph_destroy(&g_star_undirected);
    igraph_destroy(&g_ring_undirected);

    /*
     * This test is temporarily disabled until we sort out the issues with the
     * smaller test cases
     *
    printf("\nTesting undirected wheel\n");
    igraph_wheel(&g, 10, IGRAPH_WHEEL_UNDIRECTED, 0);
    // call cycles finder, expect 65 cycles to be found (
    // 1 cycle of 10 nodes,
    // 10 cycles of 9 nodes,
    // 9 cycles of 8 nodes,
    // 9 cycles of 7 nodes,
    // 9 cycles of 6 nodes,
    // 9 cycles of 5 nodes,
    // 9 cycles of 4 nodes
    // 9 cycles of 3 nodes,
    // )
    check_cycles(&g, 65);
    igraph_destroy(&g);
    */

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
    check_cycles(&g_wheel_undirected_2, 73);
    // clean up
    igraph_destroy(&g_wheel_undirected_2);

    ////////////////////////////////
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
    igraph_small(&g, 5, IGRAPH_DIRECTED, 1, 2, 2, 3, 2, 3, 3, 4, 4, 1, -1);
    check_cycles(&g, 2);
    igraph_destroy(&g);

    // same, but undirected
    printf("\nTesting undirected graph with a cycle of length 4 and a multi-edge\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED, 1, 2, 2, 3, 2, 3, 3, 4, 4, 1, -1);
    check_cycles(&g, 2);
    igraph_destroy(&g);

    // check that self-loops are handled
    printf("\nTesting graph with single self-loop\n");
    igraph_small(&g, 1, IGRAPH_DIRECTED, 0, 0, -1);
    check_cycles(&g, 1);
    igraph_destroy(&g);

    ////////////////////////////////

    // clean up test
    VERIFY_FINALLY_STACK();

    return 0;
}
