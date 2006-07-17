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

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

int main() {
  igraph_t g;
  igraph_integer_t eid;
  igraph_vector_t hist;
  long int i;
  int ret;

  /* DIRECTED */

  igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);

  igraph_vector_init(&hist, 9);

  for (i=1; i<10; i++) {
    igraph_get_eid(&g, &eid, 0, i);
    VECTOR(hist)[ (long int) eid ] = 1;
  }
  print_vector(&hist, stdout);
  
  igraph_vector_destroy(&hist);
  igraph_destroy(&g);

  /* UNDIRECTED */

  igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
  
  igraph_vector_init(&hist, 9);

  for (i=1; i<10; i++) {
    igraph_get_eid(&g, &eid, 0, i);
    VECTOR(hist)[ (long int) eid ] += 1;
    igraph_get_eid(&g, &eid, i, 0);
    VECTOR(hist)[ (long int) eid ] += 1;    
  }
  print_vector(&hist, stdout);
  
  igraph_vector_destroy(&hist);
  igraph_destroy(&g);

  /* NON-EXISTANT EDGE */
  
  igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);

  igraph_set_error_handler(igraph_error_handler_ignore);
  
  ret=igraph_get_eid(&g, &eid, 5, 6);
  if (ret != IGRAPH_EINVAL) {
    return 1;
  }
  
  igraph_destroy(&g);

  return 0;
}
