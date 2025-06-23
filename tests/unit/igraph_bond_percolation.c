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

#include "igraph_constructors.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "test_utilities.h"
#include <igraph.h>
#include "igraph_components.h"

igraph_error_t percolate(igraph_t *graph, igraph_vector_int_t *edge_indices) {
  igraph_vector_int_t outputs;
  IGRAPH_CHECK(igraph_vector_int_init(&outputs, 0));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &outputs);
  
  IGRAPH_CHECK(igraph_bond_percolation(graph, &outputs, edge_indices));

  print_vector_int(&outputs);
  
  igraph_vector_int_destroy(&outputs);
  IGRAPH_FINALLY_CLEAN(1);
  return IGRAPH_SUCCESS;
}

int main(void) {
  // test with normal graph and provided edge list

  igraph_t k_3, c_4, karate;
  igraph_full(&k_3, 3, false, false);
  
  igraph_small(&c_4, 4, IGRAPH_UNDIRECTED, 0, 1,1,2,2,3,3,0,-1);

  printf("K_3 graph percolation curve, no provided edge sequence:\n");
  IGRAPH_CHECK(percolate(&k_3, NULL));

  igraph_vector_int_t edge_ids;
  IGRAPH_CHECK(igraph_vector_int_init_int(&edge_ids, 4, 0,2,1,3));
  VECTOR(edge_ids);
  
  printf("C_4 graph with edge sequece 0, 2");
  IGRAPH_CHECK(percolate(&c_4, &edge_ids));
  
  igraph_vector_int_destroy(&edge_ids);
  igraph_destroy(&k_3);
  igraph_destroy(&c_4);

  igraph_famous(&karate, "Zachary");
  // TODO: verify that it counts the vertices correctly.

  // TODO: generate decently sized disconnected graph, verify that it doesn't still create one single component.
  igraph_destroy(&karate);

  
  VERIFY_FINALLY_STACK();
  return 0;
}
