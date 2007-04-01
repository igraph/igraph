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
  long int i, size=igraph_vector_size(v);
  for (i=0; i<size; i++) {
    printf("%g ", (double) VECTOR(*v)[i]);
  }
  printf("\n");
}

int main() {
  
  igraph_t g;
  igraph_vector_t weights, partition, partition2, cut;
  igraph_integer_t value;
  long int i;
  
  igraph_small(&g, 0, IGRAPH_UNDIRECTED, 
	       0,1, 0,4, 1,2, 1,4, 1,5, 2,3, 2,6, 3,6, 3,7, 4,5, 5,6, 6,7,
	       -1);
  
  igraph_vector_init_int_end(&weights, -1, 2,3,3,2,2, 4,2,2,2,3, 1,3, -1);
  
  igraph_vector_init(&partition, 0);
  igraph_vector_init(&partition2, 0);
  igraph_vector_init(&cut, 0);
  
  igraph_mincut(&g, &value, &partition, &partition2, &cut, &weights);
  
  printf("mincut value: %g\n", (double) value);
  
  printf("first partition:  ");
  igraph_vector_sort(&partition);
  print_vector(&partition);
  printf("second partition: ");
  igraph_vector_sort(&partition2);
  print_vector(&partition2);
  printf("edges in the cut: ");
  for (i=0; i<igraph_vector_size(&cut); i++) {
    igraph_integer_t from, to;
    igraph_edge(&g, VECTOR(cut)[i], &from, &to);
    if (from > to) {
      igraph_integer_t tmp=from;
      from=to;
      to=tmp;
    }
    printf("%li-%li ", (long int)from, (long int)to);
  }
  printf("\n");

  igraph_vector_destroy(&cut);
  igraph_vector_destroy(&partition2);
  igraph_vector_destroy(&partition);
  igraph_vector_destroy(&weights);
  igraph_destroy(&g);
 
  return 0;
}
