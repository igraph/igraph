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
#include "test_utilities.h"

int main()
{
  igraph_t g_cycle, g_no_cycle;

  // first test: star graph, does not have a cycle
  igraph_vector_int_list_t resultsT1;
  igraph_vector_int_list_t_init(&resultsT1, 0);
  igraph_star(&g_no_cycle, 7, IGRAPH_STAR_UNDIRECTED, 1);
  // call cycles finder, expect 0 cycle to be found
  igraph_simple_cylces_search_all(&g_no_cycle, &resultsT1);
  IGRAPH_ASSERT(igraph_list_size(resultsT1) == 0);
  // clean up
  igraph_vector_int_list_t_destroy(&resultsT1);
  igraph_destroy(&g_no_cycle);

  igraph_vector_int_list_t resultsT2;
  igraph_vector_int_list_t_init(&resultsT2, 0);
  igraph_ring(&g_cycle, 10, false, true, true);
  // call cycles finder, expect 1 cycle to be found
  igraph_simple_cylces_search_all(&g_no_cycle, &resultsT2);
  IGRAPH_ASSERT(igraph_list_size(resultsT2) == 1);
  // clean up
  igraph_vector_int_list_t_destroy(&resultsT2);
  igraph_destroy(&g_cycle);
}
