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
  igraph_vector_t v;
  int ret;

  igraph_vector_init(&v, 8);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=1; VECTOR(v)[3]=2;
  VECTOR(v)[4]=2; VECTOR(v)[5]=3;
  VECTOR(v)[6]=2; VECTOR(v)[7]=2;
  igraph_create(&g, &v, 0, 0);
  
  /* resize vector */
  igraph_vector_resize(&v, 2);
  VECTOR(v)[0]=3; VECTOR(v)[1]=2;
  igraph_delete_edges(&g, &v);
  if (igraph_ecount(&g) != 3) {
    return 1;
  }

  /* error test, no such edge to delete */
  igraph_set_error_handler(igraph_error_handler_ignore);
  VECTOR(v)[0]=3; VECTOR(v)[1]=2;
  ret=igraph_delete_edges(&g, &v);
  if (ret != IGRAPH_EINVAL) {
    return 2;
  } 
  if (igraph_ecount(&g) != 3) {
    return 3;
  }

  /* error test, invalid vertex id */
  VECTOR(v)[0]=10;
  ret=igraph_delete_edges(&g, &v);
  if (ret != IGRAPH_EINVVID) {
    return 4;
  } 
  if (igraph_ecount(&g) != 3) {
    return 5;
  }
  
  /* error test, invalid vertex id */
  igraph_vector_resize(&v, 3);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=2;
  ret=igraph_delete_edges(&g, &v);
  if (ret != IGRAPH_EINVEVECTOR) {
    return 6;
  } 
  if (igraph_ecount(&g) != 3) {
    return 7;
  }  

  igraph_vector_destroy(&v);
  igraph_destroy(&g);
  
  return 0;
}
