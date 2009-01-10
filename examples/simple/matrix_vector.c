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

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n; i++) {
    fprintf(f, " %li", (long int)VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

void print_matrix(igraph_matrix_t *m, FILE *f) {
  long int i, j;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    for (j=0; j<igraph_matrix_ncol(m); j++) {
      fprintf(f, " %li", (long int)MATRIX(*m, i, j));
    }
    fprintf(f, "\n");
  }  
}

void identity_matrix(igraph_matrix_t *m) {
  long int i;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    MATRIX(*m, i, i) = 1.0;
  }
}

#define RNG() (long int)(((double)rand()) / RAND_MAX * 10);

void random_vector(igraph_vector_t *v) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n; i++) {
    VECTOR(*v)[i]=RNG();
  }
}

void random_matrix(igraph_matrix_t *m) {
  long int i, j;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    for (j=0; j<igraph_matrix_ncol(m); j++) {
      MATRIX(*m, i, j) = RNG();
    }
  }  
}

int main() {
  
  igraph_matrix_t m;
  igraph_vector_t v, res;

  srand(23);

  /***/
  
  igraph_matrix_init(&m, 10, 10);
  identity_matrix(&m);
  igraph_vector_init(&v, 10);
  random_vector(&v);
  
  igraph_vector_init(&res,0);
  igraph_matrix_vector_prod(&m, &v, &res);
  if (!igraph_vector_is_equal(&v, &res)) {
    return 0;
  }

  /***/

  igraph_matrix_resize(&m, 5, 10);
  random_matrix(&m);
  random_vector(&v);
  igraph_matrix_vector_prod(&m, &v, &res);

  print_matrix(&m, stdout);
  print_vector(&v, stdout);
  print_vector(&res, stdout);

  /***/
  
  igraph_vector_destroy(&res);
  igraph_vector_destroy(&res);
  igraph_matrix_destroy(&m);
  
  return 0;
}
