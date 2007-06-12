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
  igraph_vector_t v, v2;
  int i, ret;
  
  igraph_barabasi_game(&g, 10, 2, 0, 0, 1);
  if (igraph_ecount(&g) != 18) {
    return 1;
  }
  if (igraph_vcount(&g) != 10) {
    return 2;
  }
  if (!igraph_is_directed(&g)) {
    return 3;
  }

  igraph_vector_init(&v, 0);
  igraph_get_edgelist(&g, &v, 0);
  for (i=0; i<igraph_ecount(&g); i++) {
    if (VECTOR(v)[2*i] <= VECTOR(v)[2*i+1]) {
      return 4;
    }
  }
  igraph_destroy(&g);
  
  /* out degree sequence */
  igraph_vector_resize(&v, 10);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=3; VECTOR(v)[3]=3;
  VECTOR(v)[4]=4; VECTOR(v)[5]=5;
  VECTOR(v)[6]=6; VECTOR(v)[7]=7;
  VECTOR(v)[8]=8; VECTOR(v)[9]=9;
  
  igraph_barabasi_game(&g, 10, 0, &v, 0, 1);
  if (igraph_ecount(&g) != igraph_vector_sum(&v)) {
    return 5;
  }
  igraph_vector_init(&v2, 0);
  igraph_degree(&g, &v2, igraph_vss_all(), IGRAPH_OUT, 1);
  for (i=0; i<igraph_vcount(&g); i++) {
    if (VECTOR(v)[i] != VECTOR(v2)[i]) {
      return 6;
    }
  }
  igraph_vector_destroy(&v);
  igraph_vector_destroy(&v2);
  igraph_destroy(&g);
  
  /* outpref, we cannot really test this quantitatively,
     would need to set random seed */
  igraph_barabasi_game(&g, 10, 2, 0, 1, 1);
  igraph_vector_init(&v, 0);
  igraph_get_edgelist(&g, &v, 0);
  for (i=0; i<igraph_ecount(&g); i++) {
    if (VECTOR(v)[2*i] <= VECTOR(v)[2*i+1]) {
      return 7;
    }
  }
  if (!igraph_is_directed(&g)) {
    return 8;
  }
  igraph_vector_destroy(&v);
  igraph_destroy(&g);

  /* Error tests */
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_barabasi_game(&g, -10, 1, 0, 0, 0);
  if (ret != IGRAPH_EINVAL) {
    return 9;
  }
  ret=igraph_barabasi_game(&g, 10, -2, 0, 0, 0);
  if (ret != IGRAPH_EINVAL) {
    return 10;
  }
  igraph_vector_init(&v, 9);
  ret=igraph_barabasi_game(&g, 10, 0, &v, 0, 0);
  if (ret != IGRAPH_EINVAL) {
    return 11;
  }
  igraph_vector_destroy(&v);
  
  return 0;
}
