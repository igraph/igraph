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

int print_vector(const igraph_vector_t *v) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n; i++) {
    printf("%g", (double)VECTOR(*v)[i]);
    if (i!=n-1) { printf(" "); }
  }
  printf("\n");
}

int main() {

  igraph_t g;
  igraph_vector_t merges, membership;
  long int i, j;
  igraph_vector_t idx;
  
  /* Zachary Karate club */
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
 
  /* Make one step with both methods */
  igraph_vector_init(&merges, 0);
  igraph_vector_init(&membership, 0);
  igraph_community_leading_eigenvector_naive(&g, &merges, &membership, 1);

  /* Rewrite the membership to start with community 1 if needed */
  if (VECTOR(membership)[0] != 1) {
    for (i=0; i<igraph_vector_size(&membership); i++) {
      VECTOR(membership)[i] = 1-VECTOR(membership)[i];
    }
  }

  print_vector(&merges);
  print_vector(&membership);

  igraph_community_leading_eigenvector(&g, &merges, &membership, 1);

  /* Rewrite the membership to start with community 1 if needed */
  if (VECTOR(membership)[0] != 1) {
    for (i=0; i<igraph_vector_size(&membership); i++) {
      VECTOR(membership)[i] = 1-VECTOR(membership)[i];
    }
  }

  print_vector(&merges);
  print_vector(&membership);

  printf("\n");

  /* Make all the steps */
  igraph_community_leading_eigenvector(&g, &merges, &membership, igraph_vcount(&g));

  /* Rewrite the membership vector into some canonical form */
  /* Note that after this the merges and membership vectors do
     not match. For this is not a problem for this simple check. */
  igraph_vector_init(&idx, igraph_vector_size(&membership));
  for (i=0, j=0; i<igraph_vector_size(&membership); i++) {
    long int cl=VECTOR(membership)[i];
    if (VECTOR(idx)[cl] != 0) {
      VECTOR(membership)[i]=VECTOR(idx)[cl]-1;
    } else {
      VECTOR(idx)[cl]=++j;
      VECTOR(membership)[i]=j-1;
    }
  }
  igraph_vector_destroy(&idx);

  print_vector(&merges);
  print_vector(&membership);
					 
  igraph_vector_destroy(&membership);
  igraph_vector_destroy(&merges);
  igraph_destroy(&g);
  
  return 0;
}
