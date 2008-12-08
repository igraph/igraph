/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include <stdlib.h>
#include <stdio.h>

void print_matrix(igraph_matrix_t *m, FILE *f) {
  long int i, j;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    for (j=0; j<igraph_matrix_ncol(m); j++) {
      fprintf(f, " %li", (long int)MATRIX(*m, i, j));
    }
    fprintf(f, "\n");
  }  
}

#define RNG() (long int)(((double)rand()) / RAND_MAX * 10);

void random_matrix(igraph_matrix_t *m) {
  long int i, j;

  for (i=0; i<igraph_matrix_nrow(m); i++) {
    for (j=0; j<igraph_matrix_ncol(m); j++) {
      MATRIX(*m, i, j) = RNG();
    }
  }  
}

void identity_matrix(igraph_matrix_t *m) {
  long int i, j;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    MATRIX(*m, i, i) = 1.0;
  }
}

int main() {
  
  igraph_matrix_t m1, m2, m3;

  igraph_matrix_init(&m1, 10, 10);
  igraph_matrix_init(&m2, 10, 5);
  igraph_matrix_init(&m3, 0, 0);

  identity_matrix(&m1);

  srand(13);
  random_matrix(&m2);

  /* Identity matrix from the left */
  print_matrix(&m2, stdout);
  igraph_matrix_mprod(&m1, &m2, &m3);
  print_matrix(&m3, stdout);

  igraph_matrix_crossprod(&m1, &m2, &m3);
  print_matrix(&m3, stdout);

  /* General matrices, all combination */

  printf("--\n");

  igraph_matrix_resize(&m1, 10, 5);
  igraph_matrix_resize(&m2, 5, 8);
  random_matrix(&m1);
  random_matrix(&m2);
  igraph_matrix_mprod(&m1, &m2, &m3);
  print_matrix(&m1, stdout);
  print_matrix(&m2, stdout);
  print_matrix(&m3, stdout);

  printf("--\n");
  
  igraph_matrix_resize(&m1, 5, 10);
  igraph_matrix_resize(&m2, 5, 8);
  random_matrix(&m1); 
  random_matrix(&m2);
  igraph_matrix_crossprod(&m1, &m2, &m3);
  print_matrix(&m1, stdout);
  print_matrix(&m2, stdout);
  print_matrix(&m3, stdout);

  printf("--\n");

  igraph_matrix_resize(&m1, 10, 5);
  igraph_matrix_resize(&m2, 8, 5);
  random_matrix(&m1);
  random_matrix(&m2);
  igraph_matrix_tcrossprod(&m1, &m2, &m3);
  print_matrix(&m1, stdout);
  print_matrix(&m2, stdout);
  print_matrix(&m3, stdout);

  igraph_matrix_destroy(&m1);
  igraph_matrix_destroy(&m2);
  igraph_matrix_destroy(&m3);

  return 0;
}
