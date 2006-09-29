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

int main() {
  
  igraph_t g;
  igraph_real_t flow;
  igraph_vector_t capacity;
  long int i;

  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100, 5/100, 
			  IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  
  igraph_vector_init(&capacity, igraph_ecount(&g));
  for (i=0; i<igraph_ecount(&g); i++) {
    VECTOR(capacity)[i] = 1.0;
  }
  
  for (i=0; i<10; i++) {
    igraph_maxflow(&g, &flow, 0, igraph_vcount(&g)-1, &capacity);
  }
  
  igraph_vector_destroy(&capacity);
  igraph_destroy(&g);
  
  return 0;
}
