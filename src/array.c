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

#include "types.h"

int igraph_array3_init(igraph_array3_t *a, long int n1, long int n2, 
		       long int n3) {
  int ret;
  ret=igraph_vector_init(&a->data, n1*n2*n3);
  a->n1=n1;
  a->n2=n2;
  a->n3=n3;
  a->n1n2=n1*n2;
  
  return ret;
}

void igraph_array3_destroy(igraph_array3_t *a) {
  igraph_vector_destroy(&a->data);
}

long int igraph_array3_size(const igraph_array3_t *a) {
  return (a->n1n2) * (a->n3);
}

long int igraph_array3_n(const igraph_array3_t *a, long int idx) {
  switch (idx) {
  case 1: return a->n1;
    break;
  case 2: return a->n2;
    break;
  case 3: return a->n3;
    break;
  }
  return 0;
}

int igraph_array3_resize(igraph_array3_t *a, long int n1, long int n2, 
			 long int n3) {
  int ret=igraph_vector_resize(&a->data, n1*n2*n3);
  a->n1=n1;
  a->n2=n2;
  a->n3=n3;
  a->n1n2=n1*n2;
  
  return ret;
}

void igraph_array3_null(igraph_array3_t *a) {
  igraph_vector_null(&a->data);
}
