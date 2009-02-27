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

#include "igraph.h"
#include "attributes.h"

int igraph_i_bipartite_projection(const igraph_t *graph,
				  const igraph_vector_bool_t *types,
				  igraph_t *proj,
				  int which) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i, remaining_nodes=0;
  igraph_vector_t vertex_perm, vertex_index;
  igraph_vector_t edges;
  igraph_vector_t neis;

  if (which < 0) { return 0; }
  
  IGRAPH_VECTOR_INIT_FINALLY(&vertex_perm, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&vertex_perm, no_of_nodes));
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vertex_index, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(*types)[i] == which) {
      VECTOR(vertex_index)[i] = ++remaining_nodes;
      igraph_vector_push_back(&vertex_perm, i);
    }
  }
  
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(*types)[i] != which) {
      long int j, k, n;
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_ALL));
      n=igraph_vector_size(&neis);
      for (j=0; j<n; j++) {
	long int from=VECTOR(neis)[j];
	long int from2=VECTOR(vertex_index)[from]-1;
	if (j>0 && from == VECTOR(neis)[j-1]) {
	  continue;
	}
	for (k=j+1; k<n; k++) {
	  long int to=VECTOR(neis)[k];
	  long int to2=VECTOR(vertex_index)[to]-1;
	  if (from != to) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, from2));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, to2));
	  }
	}
      }
    }
  }
  
  igraph_vector_destroy(&vertex_index);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(2);
  
  IGRAPH_CHECK(igraph_create(proj, &edges, remaining_nodes, 
			     /*directed=*/ 0));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_destroy, proj);
  
  IGRAPH_I_ATTRIBUTE_DESTROY(proj);
  IGRAPH_I_ATTRIBUTE_COPY(proj, graph, 1, 0, 0);
  /*  For this we need the new attribute handling interface first... */
  /*   IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph, proj, &vertex_perm)); */
  igraph_vector_destroy(&vertex_perm);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

/**
 * \function igraph_bipartite_projection
 * Create one or both projections of a bipartite (two-mode) network
 * 
 * TODO
 */

int igraph_bipartite_projection(const igraph_t *graph, 
				const igraph_vector_bool_t *types,
				igraph_t *proj1,
				igraph_t *proj2,
				igraph_integer_t probe1) {
  
  long int no_of_nodes=igraph_vcount(graph);

  /* t1 is -1 if proj1 is omitted, it is 0 if it belongs to type zero,
     it is 1 if it belongs to type one. The same for t2 */
  int t1, t2;
  
  if (igraph_vector_bool_size(types) != no_of_nodes) {
    IGRAPH_ERROR("Invalid bipartite type vector size", IGRAPH_EINVAL);
  }
  
  if (probe1 >= no_of_nodes) {
    IGRAPH_ERROR("No such vertex to probe", IGRAPH_EINVAL);
  }
  
  if (probe1 >= 0 && !proj1) {
    IGRAPH_ERROR("`probe1' given, but `proj1' is a null pointer", IGRAPH_EINVAL);
  }
  
  if (probe1 >=0) {
    t1=VECTOR(*types)[(long int)probe1];
    if (proj2) {
      t2=1-t1;
    } else {
      t2=-1;
    }
  } else {
    t1 = proj1 ? 0 : -1;
    t2 = proj2 ? 1 : -1;
  }
  
  IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj1, t1));
  IGRAPH_FINALLY(igraph_destroy, proj1);
  IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj2, t2));
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}
  
  
  
      
    
