/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>

void print_matrix(igraph_matrix_t* m) {
  long int nr=igraph_matrix_nrow(m);
  long int nc=igraph_matrix_ncol(m);
  long int i, j;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      if (j!=0) { putchar(' '); }
      printf("%d", (int)MATRIX(*m, i, j));
    }
    printf("\n");
  }
}

igraph_bool_t check_laplacian(igraph_t* graph, igraph_matrix_t* matrix) {
  igraph_vector_t vec, res;
  long int i, j;

  igraph_vector_init(&vec, 0);
  igraph_vector_init(&res, igraph_vcount(graph));

  igraph_degree(graph, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS);

  for (i = 0; i < igraph_vcount(graph); i++) {
    VECTOR(vec)[i] = sqrt(VECTOR(vec)[i]);
  }

  for (i = 0; i < igraph_vcount(graph); i++) {
    for (j = 0; j < igraph_vcount(graph); j++) {
      VECTOR(res)[i] += MATRIX(*matrix, i, j) * VECTOR(vec)[j];
    }
  }

  if (igraph_vector_min(&res) > 1e-7) {
    printf("Invalid Laplacian matrix:\n");
    print_matrix(matrix);
    return 0;
  }

  igraph_vector_destroy(&vec);
  igraph_vector_destroy(&res);

  return 1;
}

int test_unnormalized_laplacian(igraph_bool_t dir) {
  igraph_t g;
  igraph_matrix_t m;
  igraph_vector_t vec;
  igraph_matrix_init(&m, 1, 1);

  /* No loop or multiple edges */
  igraph_ring(&g, 5, dir, 0, 1);
  igraph_laplacian(&g, &m, 0);
  print_matrix(&m);
  printf("===\n");

  /* Add some loop edges */
  igraph_vector_init_real(&vec, 4, 1.0, 1.0, 2.0, 2.0);
  igraph_add_edges(&g, &vec, 0);
  igraph_vector_destroy(&vec);

  igraph_laplacian(&g, &m, 0);
  print_matrix(&m);
  printf("===\n");

  /* Duplicate some edges */
  igraph_vector_init_real(&vec, 4, 1.0, 2.0, 3.0, 4.0);
  igraph_add_edges(&g, &vec, 0);
  igraph_vector_destroy(&vec);

  igraph_laplacian(&g, &m, 0);
  print_matrix(&m);

  igraph_destroy(&g);

  igraph_matrix_destroy(&m);

  return 0;
}

int test_normalized_laplacian(igraph_bool_t dir) {
  igraph_t g;
  igraph_matrix_t m;
  igraph_vector_t vec;
  igraph_matrix_init(&m, 1, 1);
  igraph_bool_t ok = 1;

  /* Undirected graph, no loop or multiple edges */
  igraph_ring(&g, 5, dir, 0, 1);
  igraph_laplacian(&g, &m, 1);
  ok = ok && check_laplacian(&g, &m);

  /* Add some loop edges */
  igraph_vector_init_real(&vec, 4, 1.0, 1.0, 2.0, 2.0);
  igraph_add_edges(&g, &vec, 0);
  igraph_vector_destroy(&vec);

  igraph_laplacian(&g, &m, 1);
  ok = ok && check_laplacian(&g, &m);

  /* Duplicate some edges */
  igraph_vector_init_real(&vec, 4, 1.0, 2.0, 3.0, 4.0);
  igraph_add_edges(&g, &vec, 0);
  igraph_vector_destroy(&vec);

  igraph_laplacian(&g, &m, 1);
  ok = ok && check_laplacian(&g, &m);

  igraph_destroy(&g);

  igraph_matrix_destroy(&m);

  if (ok)
    printf("OK\n");

  return !ok;
}

int main() {
  int res;
  int i;

  for (i = 0; i < 4; i++) {
    igraph_bool_t is_normalized = i / 2;
    igraph_bool_t dir = (i % 2 ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED);

    printf("=== %sormalized, unweighted, %sdirected\n",
      (is_normalized ? "N" : "Unn"),
      (dir == IGRAPH_DIRECTED ? "" : "un")
    );

    if (is_normalized)
      res = test_normalized_laplacian(dir);
    else
      res = test_unnormalized_laplacian(dir);

    if (res)
      return i+1;
  }

  return 0;
}
