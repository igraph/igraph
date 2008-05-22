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

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

int main() {
  igraph_i_cutheap_t ch;
  long int i;
  
  igraph_i_cutheap_init(&ch, 10);

  for (i=0; i<10; i++) {
    igraph_i_cutheap_update(&ch, i, i);
  }
/*   print_vector(&ch.heap, stdout); */
/*   print_vector(&ch.index, stdout); */
/*   print_vector(&ch.hptr, stdout); */
  while (!igraph_i_cutheap_empty(&ch)) {
    long int idx=igraph_i_cutheap_popmax(&ch);
    printf("%li ", idx);
/*     print_vector(&ch.heap, stdout); */
/*     print_vector(&ch.index, stdout); */
/*     print_vector(&ch.hptr, stdout); */
/*     printf("------------\n"); */
  }
  printf("\n");

  igraph_i_cutheap_destroy(&ch);

  if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;
  
  return 0;
}
