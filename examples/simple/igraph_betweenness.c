/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <csardi@rmki.kfki.hu>
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
  igraph_vector_t bet, bet2, weights;

  /*******************************************************/

  igraph_barabasi_game(/* graph= */    &g,
		       /* n= */        1000,
		       /* m= */        3,
		       /* outseq= */   0,
		       /* outpref= */  0,
		       /* directed= */ 0);
  
  igraph_simplify(&g, /* multiple= */ 1, /* loops= */ 1);
  
  igraph_vector_init(&bet, 0);
  
  igraph_betweenness_estimate(/* graph=     */ &g,
			      /* res=       */ &bet,
			      /* vids=      */ igraph_vss_all(),
			      /* directed = */ 0,
			      /* cutoff=    */ 2,
			      /* weights=   */ 0);
  
  igraph_vector_destroy(&bet);
  igraph_destroy(&g);

  /*******************************************************/

  igraph_tree(&g, 20000, 10, IGRAPH_TREE_UNDIRECTED);
  
  igraph_vector_init(&bet, 0);
  
  igraph_betweenness_estimate(/* graph=     */ &g,
			      /* res=       */ &bet,
			      /* vids=      */ igraph_vss_all(),
			      /* directed = */ 0,
			      /* cutoff=    */ 3,
			      /* weights=   */ 0);

  igraph_vector_init(&bet2, 0);
  igraph_vector_init(&weights, igraph_ecount(&g));
  igraph_vector_fill(&weights, 1.0);
  
  igraph_betweenness_estimate(/* graph=     */ &g, 
			      /* res=       */ &bet2,
			      /* vids=      */ igraph_vss_all(),
			      /* directed = */ 0,
			      /* cutoff=    */ 3,
			      /* weights=   */ &weights);

  if (!igraph_vector_is_equal(&bet, &bet2)) {
    printf("gebasz"); 
    return 1;
  }

  igraph_vector_destroy(&bet);
  igraph_vector_destroy(&bet2);
  igraph_vector_destroy(&weights);
  igraph_destroy(&g);
  
  return 0;
}

