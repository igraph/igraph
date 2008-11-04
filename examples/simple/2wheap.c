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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <time.h>
#include <stdlib.h>

int main() {
  
  igraph_vector_t elems;
  igraph_2wheap_t Q;
  long int i;
  igraph_real_t prev=IGRAPH_INFINITY;

  srand(time(0));
  
  igraph_vector_init(&elems, 100);
  for (i=0; i<igraph_vector_size(&elems); i++) {
    VECTOR(elems)[i] = rand()/(double)RAND_MAX;
  }
 
  igraph_2wheap_init(&Q, igraph_vector_size(&elems));
  for (i=0; i<igraph_vector_size(&elems); i++) {
    igraph_2wheap_push_with_index(&Q, i, VECTOR(elems)[i]);
  }

  /*****/

  for (i=0; i<igraph_vector_size(&elems); i++) {
    if (VECTOR(elems)[i] != igraph_2wheap_get(&Q, i)) { return 1; }
  }

  /*****/
  
  for (i=0; i<igraph_vector_size(&elems); i++) {
    long int j;
    igraph_real_t tmp=igraph_2wheap_max(&Q);
    if (tmp > prev) { return 2; }
    if (tmp != igraph_2wheap_delete_max_index(&Q, &j)) { return 3; }
    if (VECTOR(elems)[j] != tmp) { return 4; }
    prev=tmp;
  }

  /*****/

  for (i=0; i<igraph_vector_size(&elems); i++) {
    igraph_2wheap_push_with_index(&Q, i, VECTOR(elems)[i]);
  }
  if (igraph_2wheap_size(&Q) != igraph_vector_size(&elems)) { return 5; }
  for (i=0; i<igraph_vector_size(&elems); i++) {
    VECTOR(elems)[i] = rand()/(double)RAND_MAX;
    igraph_2wheap_modify(&Q, i, VECTOR(elems)[i]);
  }
  for (i=0; i<igraph_vector_size(&elems); i++) {
    if (VECTOR(elems)[i] != igraph_2wheap_get(&Q, i)) { return 6; }
  }
  prev=IGRAPH_INFINITY;
  for (i=0; i<igraph_vector_size(&elems); i++) {
    long int j;
    igraph_real_t tmp=igraph_2wheap_max(&Q);
    if (tmp > prev) { return 7; }
    if (tmp != igraph_2wheap_delete_max_index(&Q, &j)) { return 8; }
    if (VECTOR(elems)[j] != tmp) { return 9; }
    prev=tmp;
  }
  if (!igraph_2wheap_empty(&Q)) { return 10; }
  if (igraph_2wheap_size(&Q) != 0) { return 11; }
  
  igraph_2wheap_destroy(&Q);
  igraph_vector_destroy(&elems);
  
  return 0;
}
