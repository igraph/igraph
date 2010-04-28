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

int vector_print(const igraph_vector_t *v) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n-1; i++) {
    printf("%g ", VECTOR(*v)[i]);
  }
  if (n>0) {
    printf("%g", VECTOR(*v)[i]);
  }
  printf("\n");
  return 0;
}

int main() {
  igraph_t g;
  igraph_vector_t outdeg, indeg, vec;
  igraph_bool_t is_simple;

  igraph_vector_init_real(&outdeg, 10, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0);
  igraph_vector_init(&indeg, 0);
  igraph_vector_init(&vec, 0);

  /* checking the simple method, undirected graphs */
  igraph_degree_sequence_game(&g, &outdeg, 0, IGRAPH_DEGSEQ_SIMPLE);
  if (igraph_is_directed(&g) || igraph_vcount(&g) != 10)
	return 1;
  if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_OUT, 1))
	return 2;
  vector_print(&vec);
  igraph_destroy(&g);

  /* checking the Viger-Latapy method, undirected graphs */
  igraph_degree_sequence_game(&g, &outdeg, 0, IGRAPH_DEGSEQ_VL);
  if (igraph_is_directed(&g) || igraph_vcount(&g) != 10)
	return 3;
  if (igraph_is_simple(&g, &is_simple) || !is_simple)
	return 4;
  if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_OUT, 0))
	return 5;
  vector_print(&vec);
  igraph_destroy(&g);

  igraph_vector_destroy(&vec);
  igraph_vector_destroy(&outdeg);
  igraph_vector_destroy(&indeg);

  return 0;
}

