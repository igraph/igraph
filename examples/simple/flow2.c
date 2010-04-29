/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>

int main() {
  
  igraph_t g;
  igraph_real_t flow_value;
  igraph_vector_t cut;
  igraph_vector_t capacity;
  igraph_vector_t partition, partition2;
  long int i, n;
  igraph_integer_t source, target;
  FILE *infile;
  igraph_real_t flow_value2=0.0;

  igraph_vector_init(&capacity, 0);

  /***************/
  infile=fopen("ak-4102.max", "r");
  igraph_read_graph_dimacs(&g, infile, 0, 0, &source, &target, &capacity,
			   IGRAPH_DIRECTED);
  fclose(infile);

  igraph_vector_init(&cut, 0);
  igraph_vector_init(&partition, 0);
  igraph_vector_init(&partition2, 0);

  igraph_maxflow(&g, &flow_value, /*flow=*/ 0, &cut, &partition,
		 &partition2, source, target, &capacity);

  if (flow_value != 8207) {
    return 1;
  }

  n=igraph_vector_size(&cut);
  for (i=0; i<n; i++) {
    long int e=VECTOR(cut)[i];
    flow_value2 += VECTOR(capacity)[e];
  }
  if (flow_value != flow_value2) {
    return 2;
  }

  igraph_destroy(&g);  
  igraph_vector_destroy(&capacity);
  igraph_vector_destroy(&cut);
  igraph_vector_destroy(&partition);
  igraph_vector_destroy(&partition2);

  return 0;
}
