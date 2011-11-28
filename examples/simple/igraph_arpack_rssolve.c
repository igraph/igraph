/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA
   
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
#include <stdlib.h>

extern int igraph_rng_inited;

/* First two eigenvalues of the adjacency matrix of a random graph */

int arpack_callback(igraph_real_t *to, const igraph_real_t *from, 
		    long int n, void *extra) {
  
  igraph_adjlist_t *adjlist=extra;
  igraph_vector_t *neis;
  long int i, j, nlen;
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(adjlist, i);
    nlen=igraph_vector_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=VECTOR(*neis)[j];
      to[i] += from[nei];
    }
  }				      
  
  return 0;
}

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %g", VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

void print_matrix(igraph_matrix_t *m, FILE* f) {
  long int i, j;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    for (j=0; j<igraph_matrix_ncol(m); j++) {
      fprintf(f, " %g", MATRIX(*m, i, j));
    }
    fprintf(f, "\n");
  }  
}

int main() {

  igraph_t g;
  igraph_adjlist_t adjlist;
  igraph_vector_t values;
  igraph_matrix_t vectors;
  igraph_arpack_options_t options;
  
  srand(1122);
  igraph_rng_inited=1;

  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100, 10.0/100,
			  /*directed=*/ 0, /*loops=*/ 0);

  igraph_adjlist_init(&g, &adjlist, IGRAPH_ALL);

  igraph_vector_init(&values, 0);
  igraph_matrix_init(&vectors, 0, 0);

  igraph_arpack_options_init(&options);
  options.n=igraph_vcount(&g);
  options.start=0;		/* random start vector */
  options.nev=2;
  options.ncv=5;
  options.which[0]='L'; options.which[1]='M';  

  igraph_arpack_rssolve(arpack_callback, &adjlist, &options, 
			/*storage=*/ 0, &values, &vectors);

  /* print_vector(&values, stdout); */
  /* print_matrix(&vectors, stdout); */

  igraph_matrix_destroy(&vectors);
  igraph_vector_destroy(&values);

  igraph_adjlist_destroy(&adjlist);
  igraph_destroy(&g);

  return 0;
}
