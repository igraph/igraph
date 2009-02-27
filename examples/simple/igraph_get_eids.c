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
  igraph_vector_t vec;
  igraph_vector_t eids, eids2;
  int ret;
  long int i;

  igraph_real_t q1[] = { 0,1, 0,1 };
  igraph_real_t q2[] = { 0,1, 0,1, 0,1 };
  igraph_real_t q3[] = { 1,0, 3,4, 1,0, 0,1, 3,4, 0,1 };

  igraph_vector_init(&eids, 0);

  /*********************************/
  igraph_small(&g, /*n=*/ 10, /*directed=*/ 1, 
	       0,1, 0,1, 1,0, 1,2, 3,4, 3,4, 3,4, 3,5, 3,7, 
	       9,8,
	       -1);

  igraph_vector_view(&vec, q1, sizeof(q1) / sizeof(igraph_real_t));
  igraph_get_eids(&g, &eids, &vec, /*directed=*/ 1);
  igraph_vector_sort(&eids);
  print_vector(&eids, stdout);

  igraph_vector_view(&vec, q2, sizeof(q2) / sizeof(igraph_real_t));
  igraph_get_eids(&g, &eids, &vec, /*directed=*/ 0);
  igraph_vector_sort(&eids);
  print_vector(&eids, stdout);

  igraph_vector_view(&vec, q2, sizeof(q2) / sizeof(igraph_real_t));
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_get_eids(&g, &eids, &vec, /*directed=*/ 1);
  if (ret != IGRAPH_EINVAL) { return 1; } 
  igraph_set_error_handler(igraph_error_handler_abort);

  igraph_destroy(&g);
  /*********************************/

  /*********************************/
  igraph_small(&g, /*n=*/10, /*directed=*/0, 
	       0,1, 1,0, 0,1, 3,4, 3,4, 5,4, 9,8, 
	       -1);
  
  igraph_vector_view(&vec, q1, sizeof(q1) / sizeof(igraph_real_t));
  igraph_get_eids(&g, &eids, &vec, /*directed=*/1);
  igraph_vector_sort(&eids);
  print_vector(&eids, stdout);

  igraph_vector_view(&vec, q3, sizeof(q3) / sizeof(igraph_real_t));
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_get_eids(&g, &eids, &vec, /*directed=*/0);
  if (ret != IGRAPH_EINVAL) { return 2; }
  igraph_set_error_handler(igraph_error_handler_abort);
  
  igraph_destroy(&g);

  /*********************************/

  igraph_vector_destroy(&eids);

  /*********************************/
  /* Speed tests */

#define NODES 10000
  igraph_barabasi_game(&g, /*n=*/ NODES, /*m=*/ 3, 
		       /*outseq=*/ 0, /*outpref=*/ 0,
		       /*directed=*/ 1);
  igraph_simplify(&g, /*multiple=*/ 1, /*loops=*/ 0);

  igraph_vector_init(&eids, NODES/2);
  igraph_random_sample(&eids, 0, igraph_ecount(&g)-1, NODES/2);
  igraph_vector_init(&vec, NODES);
  for (i=0; i<NODES/2; i++) {
    VECTOR(vec)[2*i]   = IGRAPH_FROM(&g, VECTOR(eids)[i]);
    VECTOR(vec)[2*i+1] = IGRAPH_TO(&g, VECTOR(eids)[i]);
  }
  igraph_vector_init(&eids2, 0);
  igraph_get_eids(&g, &eids2, &vec, /*directed=*/ 1);
  if (!igraph_vector_is_equal(&eids, &eids2)) {
    return 3;
  }

  /**/

  for (i=0; i<NODES/2; i++) {
    VECTOR(vec)[2*i]   = IGRAPH_TO(&g, VECTOR(eids)[i]);
    VECTOR(vec)[2*i+1] = IGRAPH_FROM(&g, VECTOR(eids)[i]);
  }
  igraph_get_eids(&g, &eids2, &vec, /*directed=*/ 0);
  if (!igraph_vector_is_equal(&eids, &eids2)) {
    return 4;
  }

  igraph_vector_destroy(&eids);
  igraph_vector_destroy(&eids2);
  igraph_vector_destroy(&vec);
  igraph_destroy(&g);
		  
  /*********************************/

  return 0;
}
