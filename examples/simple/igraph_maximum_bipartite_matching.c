/* -*- mode: C -*-  */
/* vim:set ts=2 sw=2 sts=2 et: */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

int test_graph_from_leda_tutorial() {
  /* Test graph from the LEDA tutorial:
   * http://www.leda-tutorial.org/en/unofficial/ch05s03s05.html
   */
  igraph_t graph;
  igraph_vector_bool_t types;
  igraph_vector_long_t matching;
  igraph_integer_t matching_size;
  igraph_real_t matching_weight;
  igraph_bool_t is_matching;
  int i;

  igraph_small(&graph, 0, 0,
      0, 8, 0, 12, 0, 14,
      1, 9, 1, 10, 1, 13,
      2, 8, 2, 9,
      3, 10, 3, 11, 3, 13,
      4, 9, 4, 14,
      5, 14,
      6, 9, 6, 14,
      7, 8, 7, 12, 7, 14
      , -1);
  igraph_vector_bool_init(&types, 15);
  for (i = 0; i < 15; i++)
    VECTOR(types)[i] = (i >= 8);
  igraph_vector_long_init(&matching, 0);

  igraph_maximum_bipartite_matching(&graph, &types, &matching_size,
      &matching_weight, &matching, 0);
  if (matching_size != 6) {
    printf("matching_size is %ld, expected: 6\n", (long)matching_size);
    return 1;
  }
  if (matching_weight != 6) {
    printf("matching_weight is %ld, expected: 6\n", (long)matching_weight);
    return 2;
  }
  igraph_is_maximal_matching(&graph, &types, &matching, &is_matching);
  if (!is_matching) {
    printf("not a matching: ");
    igraph_vector_long_print(&matching);
    return 3;
  }

  igraph_vector_long_destroy(&matching);
  igraph_vector_bool_destroy(&types);
  igraph_destroy(&graph);

  return 0;
}

int main() {
  if (test_graph_from_leda_tutorial())
    return 1;

  if (!IGRAPH_FINALLY_STACK_EMPTY) {
    printf("Finally stack still has %d elements.\n", IGRAPH_FINALLY_STACK_SIZE());
    return 2;
  }

  return 0;
}
