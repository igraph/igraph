/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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

int main() {
  
  igraph_t graph;
  igraph_vector_ptr_t separators;
  long int i, n;
  
  igraph_famous(&graph, "zachary");
  igraph_vector_ptr_init(&separators, 0);
  igraph_all_minimal_ab_separators(&graph, &separators);

  n=igraph_vector_ptr_size(&separators);
  for (i=0; i<n; i++) {
    igraph_bool_t res;
    igraph_vector_long_t *sep=VECTOR(separators)[i];
    igraph_is_separator(&graph, sep, &res);
    if (!res) { 
      printf("Vertex set %li is not a separator!\n", i);
      igraph_vector_long_print(sep);
      return 1;
    }
  }

  igraph_destroy(&graph);
  for (i=0; i<n; i++) {
    igraph_vector_long_t *v=VECTOR(separators)[i];
    igraph_vector_long_destroy(v);
    igraph_Free(v);
  }
  igraph_vector_ptr_destroy(&separators);
  
  return 0;
}
