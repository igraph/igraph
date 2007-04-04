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

int main() {
  igraph_t g;
  igraph_matrix_t merges;
  igraph_vector_t modularity;
  long int no_of_nodes;
  long int i;
  
  igraph_small(&g, 5, IGRAPH_UNDIRECTED, 
	       0,1,0,2,0,3,0,4, 1,2,1,3,1,4, 2,3,2,4, 3,4,
	       5,6,5,7,5,8,5,9, 6,7,6,8,6,9, 7,8,7,9, 8,9, 0,5, -1);
  igraph_vector_init(&modularity, 0);
  igraph_matrix_init(&merges, 0, 0);
  
  igraph_community_walktrap(&g, 0 /* no weights */,
			    4 /* steps */,
			    &merges, &modularity);
  
  no_of_nodes=igraph_vcount(&g);
  printf("Merges:\n");
  for (i=0; i<igraph_matrix_nrow(&merges); i++) {
    printf("%2.1li + %2.li -> %2.li (modularity %4.2f)\n", 
	   (long int)MATRIX(merges, i, 0), 
	   (long int)MATRIX(merges, i, 1), 
	   no_of_nodes+i,
	   VECTOR(modularity)[i]);
  }
  
  igraph_matrix_destroy(&merges);
  igraph_vector_destroy(&modularity);
  igraph_destroy(&g);
  
  return 0;
}
