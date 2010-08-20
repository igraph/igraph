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

int print_and_destroy(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(*ptr)[i];
    igraph_vector_print(v);
    igraph_vector_destroy(v);
    igraph_free(v);
  }
  igraph_vector_ptr_destroy(ptr);
  return 0;
}

int main() {
  igraph_t g;
  igraph_vector_ptr_t sep;
  
  igraph_small(&g, 7, IGRAPH_UNDIRECTED,
	       1,0, 2,0, 3,0, 4,0, 5,0, 6,0,
	       -1);
  igraph_vector_ptr_init(&sep, 0);
  igraph_minimum_size_separators(&g, &sep);
  print_and_destroy(&sep);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 5, IGRAPH_UNDIRECTED,
	       0,3, 1,3, 2,3,
	       0,4, 1,4, 2,4,
	       -1);
  igraph_vector_ptr_init(&sep, 0);
  igraph_minimum_size_separators(&g, &sep);
  print_and_destroy(&sep);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 5, IGRAPH_UNDIRECTED,
	       2,0, 3,0, 4,0,
	       2,1, 3,1, 4,1,
	       -1);
  igraph_vector_ptr_init(&sep, 0);
  igraph_minimum_size_separators(&g, &sep);
  print_and_destroy(&sep);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 10, IGRAPH_UNDIRECTED,
	       0,2, 0,3, 1,2, 1,3, 5,2, 5,3, 6,2, 6,3, 
	       7,2, 7,3, 8,2, 8,3, 9,2, 9,3,
	       2,4, 4,3,
	       -1);
  igraph_vector_ptr_init(&sep, 0);
  igraph_minimum_size_separators(&g, &sep);
  print_and_destroy(&sep);
  igraph_destroy(&g);  

  return 0;
}
