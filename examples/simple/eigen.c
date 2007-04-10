/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include <math.h>

int print_vector(const igraph_vector_t *v) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n; i++) {
    double out=VECTOR(*v)[i];
    printf("%5.2f", out);
    if (i!=n-1) { printf(" "); }
  }
  printf("\n");
  return 0;
}

int print_matrix(const igraph_matrix_t *m) {
  long int i, j, r=igraph_matrix_nrow(m), c=igraph_matrix_ncol(m);
  for (i=0; i<r; i++) {
    for (j=0; j<c; j++) {
      double out=MATRIX(*m, i, j);
    out=round(out*100)/100;
    if (out==0) { out=0; }
    printf("%5.2f", out);
    if (j!=c-1) { printf(" "); }
    }
    printf("\n");
  }
  return 0;
}

int main() {
  
  igraph_t g;
  igraph_matrix_t adj;
  int nodes, err, true=1;
  igraph_vector_t v;
  igraph_matrix_t ev;
  
  igraph_small(&g, 0, IGRAPH_UNDIRECTED,
	       0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
	       0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
	       0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
	       0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
	       1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
	       2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
	       2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
	       4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
	       6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
	       13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
	       18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
	       22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
	       23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
	       25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
	       28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
	       31, 32, 31, 33, 32, 33,
	       -1);

  igraph_matrix_init(&adj, 0, 0);
  igraph_get_adjacency(&g, &adj, IGRAPH_GET_ADJACENCY_BOTH);
  igraph_vector_init(&v, 0);
  igraph_matrix_init(&ev, 0, 0);
  igraph_eigen_rs(&adj, &v, &ev);

  print_vector(&v);
  print_matrix(&ev);

  igraph_eigen_rs(&adj, &v, 0);
  print_vector(&v);

  igraph_matrix_destroy(&ev);
  igraph_vector_destroy(&v);
  igraph_matrix_destroy(&adj);
  igraph_destroy(&g);  

  return 0;
}
