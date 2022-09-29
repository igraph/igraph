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

void checkGraphForNrOfCycles(igraph_t *graph, const igraph_integer_t expectedNrOfCycles)
{
  igraph_vector_int_list_t resultsT1;
  igraph_vector_int_list_init(&resultsT1, 0);
  igraph_simple_cycles_search_all(graph, &resultsT1);
  printf("Finished search, found %" IGRAPH_PRId " cylces.\n\n", igraph_vector_int_list_size(&resultsT1));
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
  igraph_t g_cycle;
  // test: one cycle, expect to find 1
  igraph_ring(&g_cycle, 10, /*directed=*/1, /*mutual=*/0, /*circular=*/1);
  printf("\nCreated ring\n");
  // call cycles finder, expect 1 cycle to be found
  checkGraphForNrOfCycles(&g_cycle, 1);
  // clean up
  igraph_destroy(&g_cycle);

  igraph_t g_no_cycle;
  // test: star graph, does not have a cycle
  igraph_star(&g_no_cycle, 7, IGRAPH_STAR_OUT, 1);
  printf("\nCreated star\n");
  // call cycles finder, expect 0 cycle to be found
  checkGraphForNrOfCycles(&g_no_cycle, 0);
  // clean up
  igraph_destroy(&g_no_cycle);

  igraph_t g_wheel;
  igraph_wheel(&g_wheel, 10, IGRAPH_WHEEL_OUT, 0);
  printf("\nCreated directed wheel\n");
  // call cycles finder, expect 1 cycle to be found
  checkGraphForNrOfCycles(&g_wheel, 1);
  // clean up
  igraph_destroy(&g_wheel);

  igraph_t g_ring_undirected;
  // test: one cycle, expect to find 1
  igraph_ring(&g_ring_undirected, 10, /*directed=*/0, /*mutual=*/0, /*circular=*/1);
  printf("\nCreated undirected ring\n");
  // call cycles finder, expect 1 cycle to be found
  checkGraphForNrOfCycles(&g_ring_undirected, 1);

  igraph_t g_star_undirected;
  // test: star graph, does not have a cycle
  igraph_star(&g_star_undirected, 7, IGRAPH_STAR_UNDIRECTED, 1);
  printf("\nCreated undirected star\n");
  // call cycles finder, expect 0 cycle to be found
  checkGraphForNrOfCycles(&g_star_undirected, 0);

  igraph_t g_ring_plus_star_undirected;
  igraph_disjoint_union(&g_ring_plus_star_undirected, &g_ring_undirected, &g_star_undirected);
  printf("\nCreated union of undirected wheel and star\n");
  // call cycles finder, expect 1 cycle to be found
  checkGraphForNrOfCycles(&g_ring_plus_star_undirected, 1);
  igraph_add_edge(&g_ring_plus_star_undirected, 7, 13); // add a random edge between the two structures to make them connected
  checkGraphForNrOfCycles(&g_ring_plus_star_undirected, 1);
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
  checkGraphForNrOfCycles(&g_wheel_undirected, 65);
  // clean up
  igraph_destroy(&g_wheel_undirected);

  // clean up test
  VERIFY_FINALLY_STACK();
  return 0;
}
