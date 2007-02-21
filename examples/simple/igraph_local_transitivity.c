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
  igraph_vs_t vertices;
  igraph_vector_t result;
  long int i;
  
  igraph_vs_seq(&vertices, 1, 101);
  igraph_barabasi_game(&g, 100000, 3, 0, 0, IGRAPH_DIRECTED);
  igraph_vector_init(&result, 0);
  
  for (i=0; i<1; i++) {
    igraph_transitivity_local_undirected2(&g, &result, igraph_vss_all());
  }

  for (i=0; i<1; i++) {
    igraph_transitivity_local_undirected4(&g, &result, igraph_vss_all());
  }
  
/*   for (i=0; i<igraph_vector_size(&result); i++) { */
/*     printf("%f ", VECTOR(result)[i]); */
/*   } */
/*   printf("\n"); */
  
  igraph_vector_destroy(&result);
  igraph_vs_destroy(&vertices);
  igraph_destroy(&g);
  
  return 0;
}
