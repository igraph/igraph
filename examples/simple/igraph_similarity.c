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

void print_matrix(igraph_matrix_t *m, FILE *f) {
  long int i, j;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    for (j=0; j<igraph_matrix_ncol(m); j++) {
      fprintf(f, " %.2f", MATRIX(*m, i, j));
    }
    fprintf(f, "\n");
  }
  fprintf(f, "==========\n");
}

int main() {
  
  igraph_t g;
  igraph_matrix_t m;
  
  igraph_small(&g, 0, IGRAPH_DIRECTED, 
	       0,1, 2,1, 2,0, 3,0, 
	       -1);
  
  igraph_matrix_init(&m, 0, 0);

  igraph_similarity_jaccard(&g, &m, igraph_vss_all(), IGRAPH_ALL, 1);
  print_matrix(&m, stdout);
  
  igraph_similarity_jaccard(&g, &m, igraph_vss_seq(1, 2), IGRAPH_ALL, 0);
  print_matrix(&m, stdout);
  
  igraph_similarity_jaccard(&g, &m, igraph_vss_all(), IGRAPH_OUT, 1);
  print_matrix(&m, stdout);
  
  igraph_similarity_jaccard(&g, &m, igraph_vss_all(), IGRAPH_IN, 0);
  print_matrix(&m, stdout);
  
  igraph_similarity_dice(&g, &m, igraph_vss_all(), IGRAPH_ALL, 1);
  print_matrix(&m, stdout);

  igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_ALL);
  print_matrix(&m, stdout);
  
  igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_OUT);
  print_matrix(&m, stdout);
  
  igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_IN);
  print_matrix(&m, stdout);
  
  igraph_matrix_destroy(&m);
  igraph_destroy(&g);
  
  return 0;
}
