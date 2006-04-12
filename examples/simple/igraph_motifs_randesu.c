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

void print_vector(igraph_vector_t *v) {
  long int i;
  igraph_real_t sum=igraph_vector_sum(v);
  for (i=0; i<igraph_vector_size(v); i++) {
    printf("%2.2f ", VECTOR(*v)[i]/sum);
  }
  printf("\n");
}

int main() {

  igraph_t g;
  igraph_vector_t hist;
  igraph_vector_t cp;

  igraph_vector_init_real(&cp, 8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  igraph_ring(&g, 1000, IGRAPH_DIRECTED, 1, 1);
  igraph_vector_init(&hist, 0);
  igraph_motifs_randesu(&g, &hist, 3, &cp);
  print_vector(&hist);
  igraph_destroy(&g);
  igraph_vector_destroy(&hist);
  igraph_vector_destroy(&cp);
  
  return 0;
}
