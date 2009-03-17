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
#include <stdlib.h>

int check_projection(const igraph_t *graph,
		     const igraph_vector_bool_t *types,
		     const igraph_t *proj1,
		     const igraph_t *proj2) {
  igraph_integer_t vcount1, ecount1, vcount2, ecount2;
  igraph_bipartite_projection_size(graph, types, &vcount1, &ecount1,
				   &vcount2, &ecount2);
  if (proj1 && igraph_vcount(proj1) != vcount1) {
    exit(10);
  }
  if (proj1 && igraph_ecount(proj1) != ecount1) {
    exit(11);
  }
  if (proj2 && igraph_vcount(proj2) != vcount2) {
    exit(12);
  }
  if (proj2 && igraph_ecount(proj2) != ecount2) {
    exit(13);
  }
  return 0;
}

int main() {
  
  igraph_t g, p1, p2, full, ring;
  igraph_vector_bool_t types;
  igraph_bool_t iso;
  long int i;
  
  /*******************************************************/
  /* Full bipartite graph -> full graphs                 */
  /*******************************************************/

  igraph_vector_bool_init(&types, 0);
  igraph_full_bipartite(&g, &types, 5, 3, /*directed=*/ 0, 
			/*mode=*/ IGRAPH_ALL);

  /* Get both projections */
  igraph_bipartite_projection(&g, &types, &p1, &p2, /*probe1=*/ -1);
  check_projection(&g, &types, &p1, &p2);

  /* Check first projection */
  igraph_full(&full, igraph_vcount(&p1), /*directed=*/0, /*loops=*/0);
  igraph_isomorphic_bliss(&p1, &full, &iso, 0, 0, 
			  IGRAPH_BLISS_FM, IGRAPH_BLISS_FM, 0, 0);
  if (!iso) { return 1; }
  igraph_destroy(&full);
  
  /* Check second projection */
  igraph_full(&full, igraph_vcount(&p2), /*directed=*/0, /*loops=*/0);
  igraph_isomorphic_bliss(&p2, &full, &iso, 0, 0, 
			  IGRAPH_BLISS_FM, IGRAPH_BLISS_FM, 0, 0);
  if (!iso) { return 2; }
  igraph_destroy(&full);
  
  igraph_destroy(&p1);
  igraph_destroy(&p2);
  igraph_destroy(&g);
  igraph_vector_bool_destroy(&types);

  /*******************************************************/
  /* More sophisticated test
  /*******************************************************/
  
  igraph_ring(&g, 100, /*directed=*/ 1, /*mutual=*/ 1, 
	      /*circular=*/ 1);
  igraph_vector_bool_init(&types, igraph_vcount(&g));
  for (i=0; i<igraph_vector_bool_size(&types); i++) {
    VECTOR(types)[i] = i%2 ? 0 : 1;
  }
  
  /* Get both projections */
  igraph_bipartite_projection(&g, &types, &p1, &p2, /*probe1=*/ -1);
  check_projection(&g, &types, &p1, &p2);
  
  /* Check first projection */
  igraph_ring(&ring, igraph_vcount(&g)/2, /*directed=*/ 0, 
	      /*mutual=*/ 0, /*circular=*/ 1);
  igraph_isomorphic_bliss(&p1, &ring, &iso, 0, 0, 
			  IGRAPH_BLISS_FM, IGRAPH_BLISS_FM, 0, 0);
  if (!iso) { return 1; }
  
  /* Check second projection */
  igraph_isomorphic_bliss(&p2, &ring, &iso, 0, 0, 
			  IGRAPH_BLISS_FM, IGRAPH_BLISS_FM, 0, 0);
  if (!iso) { return 2; }
  igraph_destroy(&ring);
  
  igraph_destroy(&p1);
  igraph_destroy(&p2);
  igraph_destroy(&g);
  igraph_vector_bool_destroy(&types);

  return 0;
}
