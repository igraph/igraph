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

int print_and_destroy(igraph_real_t value, 
		      igraph_vector_ptr_t *partitions) {
  long int i, n=igraph_vector_ptr_size(partitions);
  printf("Value: %g\n", value);
  for (i=0; i<n; i++) {
    igraph_vector_t *vec=VECTOR(*partitions)[i];
    igraph_vector_print(vec);
    igraph_vector_destroy(vec);
    igraph_free(vec);
  }  
  igraph_vector_ptr_destroy(partitions);
  
  return 0;
}

int main() {

  igraph_t g;
  igraph_vector_ptr_t partitions;
  igraph_real_t value;
  
  igraph_small(&g, 5, IGRAPH_DIRECTED,
	       0,1, 1,2, 2,3, 3,4,
	       -1);
  
  igraph_vector_ptr_init(&partitions, 0);
  igraph_all_st_mincuts(&g, &value, /*cuts=*/ 0, &partitions, 
			/*source=*/ 0, /*target=*/ 4,
			/*capacity=*/ 0);
  
  print_and_destroy(value, &partitions);
  igraph_destroy(&g);  

  /* ---------------------------------------------------------------- */
  
  igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 1,2, 1,3, 2,4, 3,4, 4,5, -1);
  igraph_vector_ptr_init(&partitions, 0);
  igraph_all_st_mincuts(&g, &value, /*cuts=*/ 0, &partitions,
			/*source=*/ 0, /*target=*/ 5, /*capacity=*/ 0);

  print_and_destroy(value, &partitions);
  igraph_destroy(&g);  

  /* ---------------------------------------------------------------- */
  
  igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 1,2, 1,3, 2,4, 3,4, 4,5, -1);
  igraph_vector_ptr_init(&partitions, 0);
  igraph_all_st_mincuts(&g, &value, /*cuts=*/ 0, &partitions,
			/*source=*/ 0, /*target=*/ 4, /*capacity=*/ 0);

  print_and_destroy(value, &partitions);
  igraph_destroy(&g);  
 
  return 0;
}
