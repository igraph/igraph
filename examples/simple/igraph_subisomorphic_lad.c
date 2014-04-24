/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>

int main() {
  igraph_t target, pattern;
  igraph_bool_t iso;
  igraph_vector_t map;
  igraph_vector_ptr_t maps;
  int i, n;
  int domainsvec[] = { 0, 2, 8, -1,
		       4, 5, 6, 7, -1, 
		       1, 3, 5, 6, 7, 8, -1,
		       0, 2, 8, -1,
		       1, 3, 7, 8, -1, -2 };
  igraph_vector_ptr_t domains;
  igraph_vector_t *v=0;

  igraph_small(&target, 9, IGRAPH_UNDIRECTED, 
	       0,1,0,4,0,6,
	       1,0,1,4,1,2,
	       2,1,2,3,
	       3,2,3,4,3,5,3,7,3,8,
	       4,0,4,1,4,3,4,5,4,6,
	       5,6,5,4,5,3,5,8,
	       6,0,6,4,6,5,
	       7,3,7,8,
	       8,5,8,3,8,7, 
	       -1);
  igraph_simplify(&target, /*multiple=*/ 1, /*loops=*/ 0, /*edge_comb=*/ 0);

  igraph_small(&pattern, 5, IGRAPH_UNDIRECTED,
	       0,1,0,4,
	       1,0,1,4,1,2,
	       2,1,2,3,
	       3,2,3,4,
	       4,3,4,1,4,0,
	       -1);
  igraph_simplify(&pattern, /*multiple=*/ 1, /*loops=*/ 0, /*edge_comb=*/ 0);
  
  igraph_vector_init(&map, 0);
  igraph_vector_ptr_init(&maps, 0);
  
  igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0, &iso, &map, 
			   &maps, /*induced=*/ 0, /*time_limit=*/ 0);

  if (!iso) { return 1; }
  igraph_vector_print(&map);
  n=igraph_vector_ptr_size(&maps);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(maps)[i];
    igraph_vector_print(v);
    igraph_vector_destroy(v);
    igraph_free(v);
  }
 
  printf("---------\n");

  igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0, &iso, &map, 
			   &maps, /*induced=*/ 1, /*time_limit=*/ 0);

  if (!iso) { return 2; }
  igraph_vector_print(&map);
  n=igraph_vector_ptr_size(&maps);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(maps)[i];
    igraph_vector_print(v);
    igraph_vector_destroy(v);
    igraph_free(v);
  }

  printf("---------\n");

  igraph_vector_ptr_init(&domains, 0);
  i=0; while (1) { 
    if (domainsvec[i] == -2) { 
      break; 
    } else if (domainsvec[i] == -1) {
      igraph_vector_ptr_push_back(&domains, v);
      v=0;
    } else {
      if (!v) { 
	v=(igraph_vector_t *) malloc(sizeof(igraph_vector_t));
	igraph_vector_init(v, 0);
      }
      igraph_vector_push_back(v, domainsvec[i]);
    }
    i++;
  }

  igraph_subisomorphic_lad(&pattern, &target, &domains, &iso, &map, &maps,
			   /*induced=*/ 0, /*time_limit=*/ 0);
  
  if (!iso) { return 3; }
  igraph_vector_print(&map);
  n=igraph_vector_ptr_size(&maps);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(maps)[i];
    igraph_vector_print(v);
    igraph_vector_destroy(v);
    igraph_free(v);
  }

  n=igraph_vector_ptr_size(&domains);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(domains)[i];
    igraph_vector_destroy(v);
    free(v);
  }
  
  igraph_vector_ptr_destroy(&domains);
  igraph_vector_destroy(&map);
  igraph_vector_ptr_destroy(&maps);

  igraph_destroy(&pattern);
  igraph_destroy(&target);
  
  return 0;
}
