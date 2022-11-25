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

void checkGraphForNrOfCycles(igraph_t *graph, const igraph_integer_t expectedNrOfCycles, igraph_simple_cycle_search_mode_t search_mode)
{
    igraph_vector_int_list_t resultsT1;
    igraph_vector_int_list_init(&resultsT1, 0);
    igraph_simple_cycles_search_all(graph, &resultsT1, search_mode);
    printf("Finished search, found %" IGRAPH_PRId " cycles.\n\n", igraph_vector_int_list_size(&resultsT1));
    if (igraph_vector_int_list_size(&resultsT1) != expectedNrOfCycles)
    {
        printf("Unexpectedly found %" IGRAPH_PRId " cycles instead of %" IGRAPH_PRId ": ", igraph_vector_int_list_size(&resultsT1), expectedNrOfCycles);
        for (int i = 0; i < igraph_vector_int_list_size(&resultsT1); i++)
        {
            print_vector_int(igraph_vector_int_list_get_ptr(&resultsT1, i));
        }
    }

    IGRAPH_ASSERT(igraph_vector_int_list_size(&resultsT1) == expectedNrOfCycles);
    igraph_vector_int_list_destroy(&resultsT1);
}

int main(void)
{
    igraph_t g_ring;
    // test: one cycle, expect to find 1
    igraph_ring(&g_ring, 10, /*directed=*/1, /*mutual=*/0, /*circular=*/1);
    printf("\nCreated ring\n");
    // call cycles finder, expect 1 cycle to be found
    checkGraphForNrOfCycles(&g_ring, 1, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
    // clean up
    igraph_destroy(&g_ring);

    igraph_t g_star;
    // test: star graph, does not have a cycle
    igraph_star(&g_star, 7, IGRAPH_STAR_OUT, 1);
    printf("\nCreated star\n");
    // call cycles finder, expect 0 cycle to be found
    checkGraphForNrOfCycles(&g_star, 0, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
    // clean up
    igraph_destroy(&g_star);

    igraph_t g_wheel;
    igraph_wheel(&g_wheel, 10, IGRAPH_WHEEL_OUT, 0);
    printf("\nCreated directed wheel\n");
    // call cycles finder, expect 1 cycle to be found
    checkGraphForNrOfCycles(&g_wheel, 1, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
    // clean up
    igraph_destroy(&g_wheel);

    igraph_t g_ring_undirected;
    // test: one cycle, expect to find 1
    igraph_ring(&g_ring_undirected, 10, /*directed=*/0, /*mutual=*/0, /*circular=*/1);
    printf("\nCreated undirected ring\n");
    // call cycles finder, expect 1 cycle to be found
    checkGraphForNrOfCycles(&g_ring_undirected, 1, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
    checkGraphForNrOfCycles(&g_ring_undirected, 2, IGRAPH_UNDIRECTED_CYCLE_SEARCH_BOTH);


    igraph_t g_star_undirected;
    // test: star graph, does not have a cycle
    igraph_star(&g_star_undirected, 7, IGRAPH_STAR_UNDIRECTED, 1);
    printf("\nCreated undirected star\n");
    // call cycles finder, expect 0 cycle to be found
    checkGraphForNrOfCycles(&g_star_undirected, 0, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);

    igraph_t g_ring_plus_star_undirected;
    igraph_disjoint_union(&g_ring_plus_star_undirected, &g_ring_undirected, &g_star_undirected);
    printf("\nCreated union of undirected wheel and star\n");
    // call cycles finder, expect 1 cycle to be found
    checkGraphForNrOfCycles(&g_ring_plus_star_undirected, 1, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
    igraph_add_edge(&g_ring_plus_star_undirected, 7, 13); // add a random edge between the two structures to make them connected
    checkGraphForNrOfCycles(&g_ring_plus_star_undirected, 1, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
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
    checkGraphForNrOfCycles(&g_wheel_undirected, 65, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
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
    checkGraphForNrOfCycles(&directed_multiedge, 2, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
    igraph_destroy(&directed_multiedge);

    // same, but undirected
    igraph_t undirected_multiedge;
    igraph_create(&undirected_multiedge, &directed_multiedge_edges, 5, false);
    // NOTE: here, if we use IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE,
    // it currently only looks at the vertices -> would not find the duplicate with different edges
    checkGraphForNrOfCycles(&undirected_multiedge, 1, IGRAPH_UNDIRECTED_CYCLE_SEARCH_ONE);
    igraph_destroy(&undirected_multiedge);
    igraph_vector_int_destroy(&directed_multiedge_edges);

    // check that self-loops are handled
    igraph_t self_loop;
    igraph_vector_int_t self_loop_edges;
    igraph_vector_int_init(&self_loop_edges, 2);
    igraph_vector_int_set(&self_loop_edges, 0, 0);
    igraph_vector_int_set(&self_loop_edges, 1, 0);
    igraph_create(&self_loop, &self_loop_edges, 1, true);
    checkGraphForNrOfCycles(&self_loop, 1, IGRAPH_UNDIRECTED_CYCLE_SEARCH_BOTH);
    igraph_destroy(&self_loop);
    igraph_vector_int_destroy(&self_loop_edges);

    ////////////////////////////////

    // clean up test
    VERIFY_FINALLY_STACK();
    return 0;
}
