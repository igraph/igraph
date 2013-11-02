/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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
#include <stdlib.h>

void print_vector(igraph_vector_t *v) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n; i++) {
    printf(" %li", (long int) VECTOR(*v)[i]);
  }
  printf("\n");
}

void warning_handler_ignore(const char* reason,const char* file,int line,int e) {
}

int main() {
  
  igraph_t g;
  igraph_vector_ptr_t result;
  igraph_es_t es;
  igraph_integer_t omega;
  long int i, j, n;
  const int params[] = {4, -1, 2, 2, 0, 0, -1, -1};
 
  igraph_set_warning_handler(warning_handler_ignore);

  igraph_vector_ptr_init(&result, 0);
  igraph_full(&g, 6, 0, 0);
  igraph_es_pairs_small(&es, 0, 0, 1, 0, 2, 3, 5, -1);
  igraph_delete_edges(&g, es);
  igraph_es_destroy(&es);
  
  for (j=0; j<sizeof(params)/(2*sizeof(params[0])); j++) {
    if (params[2*j+1] != 0) {
      igraph_cliques(&g, &result, params[2*j], params[2*j+1]);  
    } else {
      igraph_largest_cliques(&g, &result);
    }
    n = igraph_vector_ptr_size(&result);
    printf("%ld cliques found\n", (long)n);
    for (i=0; i<n; i++) {
      igraph_vector_t* v;
      v=igraph_vector_ptr_e(&result,i);
      print_vector((igraph_vector_t*)v);
      igraph_vector_destroy(v);
      free(v);
    }
  }
   
  igraph_clique_number(&g, &omega);
  printf("omega=%ld\n", (long)omega);

  igraph_destroy(&g);

  igraph_tree(&g, 5, 2, IGRAPH_TREE_OUT);
  igraph_cliques(&g, &result, 5, 5);
  if (igraph_vector_ptr_size(&result) != 0) return 1;

  igraph_destroy(&g);
  igraph_vector_ptr_destroy(&result);

  return 0;
}
