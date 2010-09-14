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

int main() {
  
  igraph_t g;
  igraph_vector_ptr_t blocks;
  igraph_vector_t cohesion;
  igraph_vector_t parent;
  igraph_t block_tree;
  long int i;

  /* The graph from the Moody-White paper */

  igraph_small(&g, 23, IGRAPH_UNDIRECTED,
	       0,1, 0,2, 0,3, 0,4, 0,5,
	       1,2, 1,3, 1,4, 1,6,
	       2,3, 2,5, 2,6,
	       3,4, 3,5, 3,6,
	       4,5, 4,6, 4,20,
	       5,6, 
	       6,7, 6,10, 6,13, 6,18,
	       7,8, 7,10, 7,13,
	       8,9,
	       9,11, 9,12,
	       10,11, 10,13,
	       11,15,
	       12,15,
	       13,14,
	       14,15,
	       16,17, 16,18, 16,19,
	       17,19, 17,20,
	       18,19, 18,21, 18,22,
	       19,20,
	       20,21, 20,22,
	       21,22,
	       -1);
  
  igraph_vector_ptr_init(&blocks, 0);
  igraph_vector_init(&cohesion, 0);
  igraph_vector_init(&parent, 0);
  igraph_cohesive_blocks(&g, &blocks, &cohesion, &parent, 
			 &block_tree);
  
  for (i=0; i<igraph_vector_ptr_size(&blocks); i++) {
    igraph_vector_t *sg=VECTOR(blocks)[i];
    igraph_vector_print(sg);
    igraph_vector_destroy(sg);
    igraph_free(sg);
  }
  igraph_vector_print(&cohesion);
  igraph_vector_print(&parent);
  igraph_write_graph_edgelist(&block_tree, stdout);

  igraph_vector_ptr_destroy(&blocks);
  igraph_vector_destroy(&cohesion);
  igraph_vector_destroy(&parent);
  igraph_destroy(&block_tree);

  igraph_destroy(&g);
  
  return 0;
}
