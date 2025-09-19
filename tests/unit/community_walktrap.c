/*
   igraph library.
   Copyright (C) 2006-2025  The igraph development team <igraph@igraph.org>

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


int main(void) {
  igraph_t graph;

  igraph_matrix_int_t merges;
  igraph_vector_t modularity;
  igraph_vector_int_t membership;
  igraph_vector_t weights;

  /* Set default seed to get reproducible results */
  igraph_rng_seed(igraph_rng_default(), 42);

  printf("Basic test\n\n");

  /* Simple unweighted graph */
  igraph_small(&graph, 3, IGRAPH_UNDIRECTED,
      0,1, 1,2, 2,0, -1);

  igraph_matrix_int_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_int_init(&membership, 0);

  igraph_community_walktrap(&graph, NULL, 4, &merges, &modularity, &membership);
  printf("Merges:\n");
  print_matrix_int(&merges);

  /* We replace modularity values which are very close to zero
   * by exact zeros, so that we can have consistent test outputs
   * across platforms. */
  printf("Modularity: ");
  igraph_vector_zapsmall(&modularity, 1e-15);
  print_vector(&modularity);

  printf("Membership: ");
  print_vector_int(&membership);

  igraph_vector_int_destroy(&membership);
  igraph_vector_destroy(&modularity);
  igraph_matrix_int_destroy(&merges);

  /* Test the case when modularity=NULL and membership=NULL as this caused a crash in
   * the R interface, see https://github.com/igraph/rigraph/issues/289 */

  printf("\nWithout modularity and membership calculation\n\n");

  igraph_matrix_int_init(&merges, 0, 0);

  igraph_community_walktrap(&graph, NULL, 4, &merges, NULL, NULL);
  printf("Merges:\n");
  print_matrix_int(&merges);
  igraph_matrix_int_destroy(&merges);

  /* Test the case when merges=NULL but modularity is requested,
   * as the result was previously invalid due to not incrementing
   * the "merge index" during computation. */

  printf("\nWithout merges matrix calculation\n\n");

  igraph_vector_init(&modularity, 0);
  igraph_community_walktrap(&graph, NULL, 4, NULL, &modularity, NULL);
  printf("Modularity:\n");
  igraph_vector_zapsmall(&modularity, 1e-15);
  print_vector(&modularity);
  igraph_vector_destroy(&modularity);

  igraph_destroy(&graph);

  /* Test for bug https://github.com/igraph/igraph/issues/2042 */

  printf("\nBug 2042\n\n");

  igraph_small(&graph, 3, IGRAPH_UNDIRECTED,
               0,1, -1);

  igraph_matrix_int_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_int_init(&membership, 0);
  igraph_vector_init_real(&weights, 1, 0.2);

  igraph_community_walktrap(&graph, &weights, 4, &merges, &modularity, &membership);
  printf("Merges:\n");
  print_matrix_int(&merges);

  printf("Modularity: ");
  igraph_vector_zapsmall(&modularity, 1e-15);
  print_vector(&modularity);

  printf("Membership: ");
  print_vector_int(&membership);

  igraph_vector_destroy(&weights);
  igraph_vector_int_destroy(&membership);
  igraph_vector_destroy(&modularity);
  igraph_matrix_int_destroy(&merges);

  igraph_destroy(&graph);

  printf("\nSmall weighted graph\n\n");

  igraph_ring(&graph, 6, IGRAPH_UNDIRECTED, 0, 1);

  igraph_matrix_int_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_int_init(&membership, 0);
  igraph_vector_init_real(&weights, 6,
                          1.0, 0.5, 0.25, 0.75, 1.25, 1.5);

  igraph_community_walktrap(&graph, &weights, 4, &merges, &modularity, &membership);
  printf("Merges:\n");
  print_matrix_int(&merges);

  printf("Modularity: ");
  igraph_vector_zapsmall(&modularity, 1e-15);
  print_vector(&modularity);

  printf("Membership: ");
  print_vector_int(&membership);

  /* Negative weights are not allowed */
  VECTOR(weights)[0] = -0.1;
  CHECK_ERROR(igraph_community_walktrap(&graph, &weights, 4, &merges, &modularity, &membership), IGRAPH_EINVAL);

  /* Invalid weight vector size */
  VECTOR(weights)[0] = 0.1;
  igraph_vector_pop_back(&weights);
  CHECK_ERROR(igraph_community_walktrap(&graph, &weights, 4, &merges, &modularity, &membership), IGRAPH_EINVAL);

  igraph_vector_destroy(&weights);
  igraph_vector_int_destroy(&membership);
  igraph_vector_destroy(&modularity);
  igraph_matrix_int_destroy(&merges);

  igraph_destroy(&graph);

  printf("\nIsolated vertices\n\n");

  igraph_empty(&graph, 5, IGRAPH_UNDIRECTED);

  igraph_matrix_int_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_int_init(&membership, 0);

  igraph_community_walktrap(&graph, NULL, 4, &merges, &modularity, &membership);
  printf("Merges:\n");
  print_matrix_int(&merges);

  printf("Modularity: ");
  igraph_vector_zapsmall(&modularity, 1e-15);
  print_vector(&modularity);

  printf("Membership: ");
  print_vector_int(&membership);

  igraph_vector_destroy(&weights);
  igraph_vector_int_destroy(&membership);
  igraph_vector_destroy(&modularity);
  igraph_matrix_int_destroy(&merges);

  igraph_destroy(&graph);


  VERIFY_FINALLY_STACK();

  return 0;
}
