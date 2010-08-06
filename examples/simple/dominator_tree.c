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
  igraph_t g, domtree;
  igraph_vector_t dom;

  igraph_vector_init(&dom, 0);
  igraph_small(&g, 13, IGRAPH_DIRECTED,
	       0,1, 0,2, 0,3,
	       1,4,
	       2,1, 2,4, 2,5,
	       3,6, 3,7,
	       4,12,
	       5,8, 
	       6,9,
	       7,9, 7,10,
	       8,5, 8,11,
	       9,11,
	       10,9,
	       11,9, 11,0,
	       12,8,
	       -1);

  /* Check NULL vector arguments */
  igraph_dominator_tree(&g, /*root=*/ 0, /*dom=*/ 0, /*domtree=*/ 0,
			/*leftout=*/ 0, /*mode=*/ IGRAPH_OUT);

  /* Proper calculation */
  igraph_dominator_tree(&g, /*root=*/ 0, &dom, /*domtree=*/ 0,
			/*leftout=*/ 0, /*mode=*/ IGRAPH_OUT);
  igraph_vector_print(&dom);

  /* Tree calculation */
  igraph_dominator_tree(&g, /*root=*/ 0, /*dom=*/ 0, /*domtree=*/ &domtree,
			/*leftout=*/ 0, /*mode=*/ IGRAPH_OUT);
  igraph_write_graph_edgelist(&domtree, stdout);
  
  igraph_vector_destroy(&dom);
  igraph_destroy(&domtree);
  igraph_destroy(&g);

  return 0;
}
