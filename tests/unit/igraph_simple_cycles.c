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

void check_cycle_count_in_graph(igraph_t *graph, const igraph_integer_t expectedNrOfCycles)
{
    igraph_vector_int_list_t results_v;
    igraph_vector_int_list_t results_e;
    igraph_vector_int_list_init(&results_v, 0);
    igraph_vector_int_list_init(&results_e, 0);
    igraph_simple_cycles_search_all(graph, &results_v, &results_e);
    printf("Finished search, found %" IGRAPH_PRId " cycles.\n\n", igraph_vector_int_list_size(&results_v));
    if (igraph_vector_int_list_size(&results_v) != expectedNrOfCycles)
    {
        printf("Unexpectedly found %" IGRAPH_PRId " cycles instead of %" IGRAPH_PRId ": ", igraph_vector_int_list_size(&results_v), expectedNrOfCycles);
        for (int i = 0; i < igraph_vector_int_list_size(&results_v); i++)
        {
            print_vector_int(igraph_vector_int_list_get_ptr(&results_v, i));
        }
    }

    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_v) == expectedNrOfCycles);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_e) == igraph_vector_int_list_size(&results_v));
    igraph_vector_int_list_destroy(&results_v);
    igraph_vector_int_list_destroy(&results_e);
}

int main(void)
{
    igraph_t g_ring;
    // test: one cycle, expect to find 1
    igraph_ring(&g_ring, 10, /*directed=*/1, /*mutual=*/0, /*circular=*/1);
    printf("\nCreated ring\n");
    // call cycles finder, expect 1 cycle to be found
    check_cycle_count_in_graph(&g_ring, 1);
    // clean up
    igraph_destroy(&g_ring);

    igraph_t g_star;
    // test: star graph, does not have a cycle
    igraph_star(&g_star, 7, IGRAPH_STAR_OUT, 1);
    printf("\nCreated star\n");
    // call cycles finder, expect 0 cycle to be found
    check_cycle_count_in_graph(&g_star, 0);
    // clean up
    igraph_destroy(&g_star);

    igraph_t g_wheel;
    igraph_wheel(&g_wheel, 10, IGRAPH_WHEEL_OUT, 0);
    printf("\nCreated directed wheel\n");
    // call cycles finder, expect 1 cycle to be found
    check_cycle_count_in_graph(&g_wheel, 1);
    // clean up
    igraph_destroy(&g_wheel);

    igraph_t g_ring_undirected;
    // test: one cycle, expect to find 1
    igraph_ring(&g_ring_undirected, 10, /*directed=*/0, /*mutual=*/0, /*circular=*/1);
    printf("\nCreated undirected ring\n");
    // call cycles finder, expect 1 cycle to be found
    check_cycle_count_in_graph(&g_ring_undirected, 1);


    igraph_t g_star_undirected;
    // test: star graph, does not have a cycle
    igraph_star(&g_star_undirected, 7, IGRAPH_STAR_UNDIRECTED, 1);
    printf("\nCreated undirected star\n");
    // call cycles finder, expect 0 cycle to be found
    check_cycle_count_in_graph(&g_star_undirected, 0);

    igraph_t g_ring_plus_star_undirected;
    igraph_disjoint_union(&g_ring_plus_star_undirected, &g_ring_undirected, &g_star_undirected);
    printf("\nCreated union of undirected wheel and star\n");
    // call cycles finder, expect 1 cycle to be found
    check_cycle_count_in_graph(&g_ring_plus_star_undirected, 1);
    igraph_add_edge(&g_ring_plus_star_undirected, 7, 13); // add a random edge between the two structures to make them connected
    check_cycle_count_in_graph(&g_ring_plus_star_undirected, 1);
    // clean up
    igraph_destroy(&g_ring_plus_star_undirected);
    igraph_destroy(&g_star_undirected);
    igraph_destroy(&g_ring_undirected);

    igraph_t g_wheel_undirected;
    igraph_wheel(&g_wheel_undirected, 10, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("\nCreated undirected wheel\n");
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
    check_cycle_count_in_graph(&g_wheel_undirected, 65);
    // clean up
    igraph_destroy(&g_wheel_undirected);

    ////////////////////////////////
    // Tests as requested in https://github.com/igraph/igraph/pull/2181#issuecomment-1326064152
    igraph_t directed_multiedge;
    igraph_vector_int_t directed_multiedge_edges;
    /**
     * This graph looks like:
     *
     *  /\
     * 1--2
     * |  |
     * 3--4
     *
     */
    igraph_vector_int_init(&directed_multiedge_edges, 10);
    igraph_vector_int_set(&directed_multiedge_edges, 0, 1);
    igraph_vector_int_set(&directed_multiedge_edges, 1, 2);
    igraph_vector_int_set(&directed_multiedge_edges, 2, 2);
    igraph_vector_int_set(&directed_multiedge_edges, 3, 3);
    igraph_vector_int_set(&directed_multiedge_edges, 4, 2);
    igraph_vector_int_set(&directed_multiedge_edges, 5, 3);
    igraph_vector_int_set(&directed_multiedge_edges, 6, 3);
    igraph_vector_int_set(&directed_multiedge_edges, 7, 4);
    igraph_vector_int_set(&directed_multiedge_edges, 8, 4);
    igraph_vector_int_set(&directed_multiedge_edges, 9, 1);
    igraph_create(&directed_multiedge, &directed_multiedge_edges, 5, true);
    check_cycle_count_in_graph(&directed_multiedge, 2);
    igraph_destroy(&directed_multiedge);

    // same, but undirected
    igraph_t undirected_multiedge;
    igraph_create(&undirected_multiedge, &directed_multiedge_edges, 5, false);
    check_cycle_count_in_graph(&undirected_multiedge, 1);
    igraph_destroy(&undirected_multiedge);
    igraph_vector_int_destroy(&directed_multiedge_edges);

    // check that self-loops are handled
    igraph_t self_loop;
    igraph_vector_int_t self_loop_edges;
    igraph_vector_int_init(&self_loop_edges, 2);
    igraph_vector_int_set(&self_loop_edges, 0, 0);
    igraph_vector_int_set(&self_loop_edges, 1, 0);
    igraph_create(&self_loop, &self_loop_edges, 1, true);
    check_cycle_count_in_graph(&self_loop, 1);
    igraph_destroy(&self_loop);
    igraph_vector_int_destroy(&self_loop_edges);

    ////////////////////////////////

    // clean up test
    VERIFY_FINALLY_STACK();
    return 0;
}
