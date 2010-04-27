/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, 1005 Lausanne, Switzerland
   
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
  igraph_vector_t weight;
  igraph_attribute_combination_t comb;
  
  igraph_i_set_attribute_table(&igraph_cattribute_table);
  
  igraph_small(&g, 4, IGRAPH_DIRECTED, 
	       0, 1, 0, 1, 0, 1,
	       1, 2, 2, 3, 
	       -1);
  
  igraph_vector_init_seq(&weight, 1, igraph_ecount(&g));
  SETEANV(&g, "weight", &weight);
  igraph_vector_destroy(&weight);

  igraph_attribute_combination(&comb, 
			       "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
			       "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE, 
			       0);
  igraph_simplify(&g, /*multiple=*/ 1, /*loops=*/ 1, &comb);
  igraph_attribute_combination_destroy(&comb);
  
  igraph_write_graph_graphml(&g, stdout);
  
  igraph_destroy(&g);

  return 0;
}
