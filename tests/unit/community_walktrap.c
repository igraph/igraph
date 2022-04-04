/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int main() {
  igraph_t graph;

  igraph_matrix_t merges;
  igraph_vector_t modularity;
  igraph_vector_t membership;

  /* Set default seed to get reproducible results */
  igraph_rng_seed(igraph_rng_default(), 42);

  /* Simple unweighted graph */
  igraph_small(&graph, 3, IGRAPH_UNDIRECTED,
      0,1, 1,2, 2,0, -1);

  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_init(&membership, 0);

  igraph_community_walktrap(&graph, NULL, 4, &merges, &modularity, &membership);
  printf("Merges:\n");
  igraph_matrix_print(&merges);

  printf("Modularity: ");
  igraph_vector_print(&modularity);

  printf("Membership: ");
  igraph_vector_print(&membership);

  igraph_vector_destroy(&membership);
  igraph_vector_destroy(&modularity);
  igraph_matrix_destroy(&merges);

  /* Test the case when modularity=0 and membership=0 as this caused a crash in
   * the R interface, see https://github.com/igraph/rigraph/issues/289 */

  igraph_matrix_init(&merges, 0, 0);

  igraph_community_walktrap(&graph, NULL, 4, &merges, NULL, NULL);
  printf("Merges:\n");
  igraph_matrix_print(&merges);
  igraph_matrix_destroy(&merges);

  igraph_destroy(&graph);

  VERIFY_FINALLY_STACK();

  return 0;
}
